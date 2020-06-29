/*
 * Copyright (C) 2014-2016 Brian L. Browning
 * Copyright (C) 2019 Altti I. Maarala
 *
 * This file is part of SparkBeagle
 *
 * SparkBeagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SparkBeagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.ngseq.sparkbeagle;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.readers.TabixReader;
import org.apache.commons.cli.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.compress.CompressionCodecFactory;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.sql.SparkSession;
import org.ngseq.sparkbeagle.beagleutil.ChromInterval;
import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.blbutil.*;
import org.ngseq.sparkbeagle.haplotype.BitHapPair;
import org.ngseq.sparkbeagle.haplotype.HapPairPhasedGT;
import org.ngseq.sparkbeagle.imp.ImpData;
import org.ngseq.sparkbeagle.imp.ImpLS;
import org.ngseq.sparkbeagle.imp.StateProbs;
import org.ngseq.sparkbeagle.vcf.*;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.function.Supplier;

/**
 * Class {@code Main} is the entry class for the Beagle program.
 * See {@code Par.usage()} and online program documentation for usage
 * instructions.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Main {

    private static SparkSession ss;
    private static JavaSparkContext sc;

    /**
     * The program name and version.
     */
    public static final String VERSION = "(version 5.0)";
    public static final String PROGRAM = "beagle.28Sep18.793.jar";
    public static final String COMMAND = "java -jar beagle.28Sep18.793.jar";

    /**
     * The copyright string.
     */
    public static final String COPYRIGHT = "Copyright (C) 2014-2018 Brian L. Browning";

    /**
     * The program name and a brief help message.
     */
    public static final String SHORT_HELP = Main.PROGRAM + " " + VERSION
            + Const.nl + Main.COPYRIGHT
            + Const.nl + "Enter \"java -jar beagle.28Sep18.793.jar\" to "
            + "list command line argument";

    private final Par par;
    private final GeneticMap genMap;
    private final Data data;
    private static RunStats runStats = null;
    private final WindowWriter windowWriter;

    /**
     * Entry point to Beagle program.  See {@code Parameters.usage()} and
     * online program documentation for usage instructions.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {

        ss = SparkSession.builder().appName("SparkBeagle").getOrCreate();
        sc = new JavaSparkContext(ss.sparkContext());

        Options options = new Options();
        options.addOption(new Option( "gt",true, "HDFS path to target data" ));
        options.addOption(new Option( "gl",true, "gl" ));
        options.addOption(new Option( "ref",true, "HDFS path to reference data" ));
        options.addOption(new Option( "map",true, "HDFS path to genetic map data" ));
        options.addOption(new Option( "out",true, "HDFS output folder" ));
        options.addOption(new Option( "nthreads",true, "Number of threads used per imputation interval" ));
        options.addOption(new Option( "window",true, "imputation window size, default is 10 cM" ));
        options.addOption(new Option( "overlap",true, "imputation window overlap size, default is 1 cM" ));

        options.addOption(new Option( "niterations",true, "Phasing iterations, default is 0" ));
        options.addOption(new Option( "rheader",true, "vcf reference header" ));
        options.addOption(new Option( "numparts",true, "Number of partitions to be used" ));
        options.addOption(new Option( "gprobs",true, "gprobs true or false" ));
        options.addOption(new Option( "phaseits",true, "Phasing iterations, 0 is default and phasing is skipped" ));
        options.addOption(new Option( "burnits",true, "Burnin iterations, 0 is default" ));
        options.addOption(new Option( "chrom",true, "chromosome to be processed" ));
        options.addOption(new Option( "hdfs",true,"HDFS url" ));
        options.addOption(new Option( "index",true,"tbi index file" ));
        options.addOption(new Option( "cluster",true,"cluster length in cM" ));


        //options.addOption(new Option( "mode",true, "gt, gl, gtgl or ref" ));
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp( "spark-submit <spark specific args>", options, true );

        CommandLineParser parser = new BasicParser();
        CommandLine cmd = null;

        try {
            // parse the command line arguments
            cmd = parser.parse( options, args );

        }
        catch( ParseException exp ) {
            // oops, something went wrong
            System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
        }

	    Locale.setDefault(Locale.US);

        Par par = parameters(cmd);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", String.valueOf(par.nthreads()));
        RunStats runStats = new RunStats(par);
        runStats.printStartInfo();

        phaseData(cmd, par);
    }

    private Main(Par par, Data data, WindowWriter windWriter, RunStats runStats) {
        assert par!=null;
        assert data!=null;
        assert windWriter!=null;
        assert runStats!=null;
        this.par = par;
        this.genMap = data.genMap();
        this.data = data;
        this.runStats = runStats;
        this.windowWriter = windWriter;
    }

    private static Data data(Par par, ByteArrayOutputStream split) {
        Filter<String> sFilter = null;
        Filter<Marker> mFilter = null;

        ChromInterval chromInterval = par.chromInt();
        if (par.ref()==null) {
            Supplier<SampleFileIt<GTRec>> targSupplier =
                    () -> targIt(par, mFilter, sFilter, chromInterval);
            return TargetData.targetData(par, targSupplier);
        }
        else {
            SampleFileIt<GTRec> targIt = targIt(par, mFilter, sFilter, chromInterval);
            Supplier<SampleFileIt<RefGTRec>> refSupplier = refSupplier(par,
                    mFilter, sFilter, chromInterval,split);
            return AllData.allData(refSupplier, targIt, par);
        }
    }


    /*
     * Phases the data, imputes ungenotyped markers, and performed IBD segment
     * detection.
     */
     private static void phaseData(CommandLine cmd, Par parameters) {

         String header = parameters.getHeader();
         Broadcast<String> rheaderBc = sc.broadcast(header);

         List<String> intervals = getIntervals(parameters);

         Broadcast<CommandLine> params = sc.broadcast(cmd);

         JavaRDD<String> intervalsRDD = sc.parallelize(intervals, intervals.size());

         intervalsRDD.foreach(i -> {

             String refheader = rheaderBc.getValue();
             Par par = parameters(params.getValue());

             String[] is = i.split(":");
             String[] intvals = is[1].split("-");
             Interval interval = new Interval(is[0],Integer.valueOf(intvals[0]), Integer.valueOf(intvals[1]));
             int prevEnd = Integer.valueOf(intvals[2]);

             TabixReader tabixReader = null;
             try {
                 tabixReader = new TabixReader(par.getRefPath(), par.getIndexPath());
             } catch (IOException e) {
                 e.printStackTrace();
             }
             TabixReader.Iterator tit = tabixReader.query(interval.getContig(), interval.getStart(), interval.getEnd());
             int overlapcnt = 0;

             ByteArrayOutputStream bos = new ByteArrayOutputStream(1 << 26);
             BlockCompressedOutputStream bao = new BlockCompressedOutputStream(bos, new File(""));

             try {

                 boolean hasnext = true;
                 bao.write(refheader.trim().getBytes());
                 bao.write("\n".getBytes());

                 while(hasnext) {
                     String tabix = tit.next();

                     if(tabix == null)
                         hasnext = false;
                     else{
                         int pos = Integer.valueOf(tabix.substring(0,30).split("\t")[1]);
                         if(pos <= prevEnd)
                             overlapcnt++;
                         bao.write(tabix.trim().getBytes());
                         bao.write("\n".getBytes());
                     }
                 }

             } catch (IOException e) {
                 e.printStackTrace();
             }
             bos.close();
             bao.close();

            Data data = data(par, bos);

            RunStats runStats = new RunStats(par);
            runStats.printSampleSummary(data);

            WindowWriter windowWriter = new WindowWriter(data.targetSamples(), par.hdfsout(), par.hdfs());

            MainHelper mh = new MainHelper(par, data.genMap(), runStats);
            GT overlapHaps = null;
            int window = 0;

            do {
                if (++window > 1) {
                    data.advanceWindowCm();
                }

                runStats.printWindowUpdate(data);

                CurrentData cd = new CurrentData(par, data.genMap(), data, overlapHaps);
                GT phasedTarg = mh.phase(cd);

                printOutput(cd, phasedTarg, windowWriter, par, data.genMap(), prevEnd, overlapcnt);

                overlapHaps = overlapHaps(cd, phasedTarg);
            } while (data.canAdvanceWindow());

        });
    }

    private static List<String> getIntervals(Par parameters) {

        List<String> intervals = new ArrayList<>();

        int chrom = Integer.valueOf(parameters.chromInt().chrom());

        PlinkGenMap pgenMap = PlinkGenMap.fromPlinkMapFile(parameters.map(), String.valueOf(chrom));
        int intvStart = pgenMap.basePos(chrom,0);
        int intvEnd = pgenMap.basePos(chrom, parameters.window());
        int lastGenPos = pgenMap.index2BasePos(chrom, pgenMap.nMapPositions(chrom)-1);
        int prevEnd = 0;

        intervals.add(chrom+":"+intvStart+"-"+intvEnd+"-0");

        while (intvEnd < lastGenPos){

            prevEnd = intvEnd;
            intvStart=pgenMap.basePos(chrom,pgenMap.genPos(chrom,prevEnd)-parameters.overlap()); //substracts overlap from previous end position
            intvEnd = pgenMap.basePos(chrom, pgenMap.genPos(chrom,intvStart)+parameters.window());
            intervals.add(chrom+":"+intvStart+"-"+intvEnd +"-"+prevEnd);

        }
        return intervals;
    }

    /*
     * Initialize GenotypeValues to have values 1 at known genotypes.
     */
    private static void initializeGV(GenotypeValues gv, GT gl) {
        assert gv.markers().equals(gl.markers());
        assert gv.samples().equals(gl.samples());
        int nMarkers = gl.nMarkers();
        int nSamples = gl.nSamples();
        for (int m=0; m<nMarkers; ++m) {
            for (int s=0; s<nSamples; ++s) {
                int a1 = gl.allele1(m, s);
                int a2 = gl.allele2(m, s);
                if (a1>=0 && a2>=0) {
                    int gt = VcfRecord.gtIndex(a1, a2);
                    gv.add(m, s, gt, 1.0);
                }
            }
        }
    }

    private static void printOutput(CurrentData cd, GT phasedTarg, WindowWriter windowWriter, Par par, GeneticMap genMap, int prevEnd, int overlap) {
        assert par.gt()!=null;

        int refStart = overlap;
        int refEnd = cd.nextSpliceStart();
        int nThreads = par.nthreads();
        if (cd.nMarkers()==cd.nTargMarkers() || par.impute() == false) {
            windowWriter.print(phasedTarg, refStart, refEnd, nThreads);
        }
        else {
            long t0 = System.nanoTime();
            ImpData impData = new ImpData(par, cd, phasedTarg, genMap);
            AtomicReferenceArray<StateProbs> stateProbs = ImpLS.stateProbs(impData);
            windowWriter.print(impData, stateProbs, refStart, refEnd, prevEnd, cd.markers().marker(0).pos());
            runStats.imputationNanos(System.nanoTime() - t0);
            runStats.printImputationUpdate();
        }
    }

    private static GT overlapHaps(CurrentData cd, GT phasedTarg) {
        assert phasedTarg.isPhased();
        int nextOverlap = cd.nextTargetOverlapStart();
        int nextSplice = cd.nextTargetSpliceStart();
        if (cd.nextOverlapStart() == cd.nextSpliceStart()) {    // xxx change 2nd to cd.nMarkers()?
            return null;
        }
        int nSamples = phasedTarg.nSamples();
        int nMarkers = nextSplice - nextOverlap;
        Markers markers = phasedTarg.markers().restrict(nextOverlap, nextSplice);
        Samples samples = phasedTarg.samples();
        List<BitHapPair> list = new ArrayList<>(nSamples);
        int[] a1 = new int[nMarkers];
        int[] a2 = new int[nMarkers];
        for (int s = 0; s < nSamples; ++s) {
            for (int m = 0; m < nMarkers; ++m) {
                a1[m] = phasedTarg.allele1(nextOverlap + m, s);
                a2[m] = phasedTarg.allele2(nextOverlap + m, s);
            }
            list.add(new BitHapPair(markers, samples.idIndex(s), a1, a2));
        }
        return new HapPairPhasedGT(phasedTarg.samples(), list);
    }

        private static SampleFileIt<GTRec> targIt(Par par,
                                              Filter<Marker> markerFilter, Filter<String> sampleFilter,
                                              ChromInterval chromInterval) {
        CompressionCodecFactory factory = new CompressionCodecFactory( new Configuration());

        FileIt<String> it = InputIt.fromGzipFile(par.gt(), factory.getCodec(new Path(par.gtpath())), (par.gtpath().endsWith(".gz") ||  par.gtpath().endsWith(".zip") ||  par.gtpath().endsWith(".tar")));

        SampleFileIt<GTRec> targIt = VcfIt.create(it, sampleFilter,
                markerFilter,  VcfIt.toBitSetGT);
        if (chromInterval!=null) {
            targIt = new IntervalVcfIt<>(targIt, chromInterval);
        }
        return targIt;
    }

    private static Supplier<SampleFileIt<RefGTRec>> refSupplier(Par par,
                                                                Filter<Marker> mFilter, Filter<String> sFilter,
                                                                ChromInterval chromInt, ByteArrayOutputStream split) {
        Filter<Marker> mFilter2 = updateFilter(par, mFilter, sFilter, chromInt);
        return () -> {
            SampleFileIt<RefGTRec> refIt;
            String filename = par.ref().toString();
            if (filename.endsWith(".bref")) {
                String s = Const.nl + "ERROR: bref format (.bref) is not supported"
                         + Const.nl + "       Reference files should be in bref3 format (.brer3)" ;
                Utilities.exit(s);
            }
            //TODO: implement support for Bref
            /*if (par.isBref()) {
                refIt = new Bref3It(par.ref(), mFilter2);
            }*/
            /*else {
                if (filename.endsWith(".vcf")==false
                        && filename.endsWith(".vcf.gz")==false) {
                    runStats.println(Const.nl
                            + "WARNING: unrecognized reference file type "
                            + "(expected \".bref3\", \".vcf\", or \".vcf.gz\")"
                            + Const.nl);
                }*/

            //FileIt<String> it = InputIt.fromGzipFile(par.ref());

            FileIt<String> it = InputIt.fromBGZFStream(split, false);

            refIt = RefIt.create(it, sFilter, mFilter2,
                        RefIt.MAX_EM_BUFFER_SIZE);
            //}
            if (chromInt!=null) {
                refIt = new IntervalVcfIt<>(refIt, chromInt);
            }
            return refIt;
        } ;
    }

    private static Filter<Marker> updateFilter(Par par, Filter<Marker> mFilter,
                                               Filter<String> sFilter, ChromInterval chromInt) {
        if (par.impute() && par.gt()!=null) {
            return mFilter;
        }
        else {
            Set<Marker> includedMarkers = new HashSet<>(50000);
            try (SampleFileIt<GTRec> vcfIt = targIt(par, mFilter, sFilter,
                    chromInt)) {
                while (vcfIt.hasNext()) {
                    includedMarkers.add(vcfIt.next().marker());
                }
            }
            return Filter.includeFilter(includedMarkers);
        }
    }

    /*
     * Checks that certain parameters are consistent, and prints error
     * message and exits if parameters are inconsistent.
     *
     * @param args the command line arguments.
     */
    private static Par parameters(CommandLine cmd) {
        // warnings are printed in RunStats.startInfo() method
        Par par = new Par(cmd);
        //checkOutputPrefix(par);
        if (1.1*par.overlap() >= par.window()) {
            String s = SHORT_HELP + Const.nl
                    + Const.nl + "ERROR: The \"window\" parameter must be at least "
                    + "1.1 times the \"overlap\" parameter"
                    + Const.nl + "Exiting program.";
            Utilities.exit(s);
        }
        return par;
    }

    private static void checkOutputPrefix(Par par) {
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String s = "ERROR: \"out\" parameter cannot be a directory: \""
                    + par.out() + "\"";
            Utilities.exit(Par.usage() + s);
        }

        File vcfOut = new File(par.out() + ".vcf.gz");
        if (vcfOut.equals(par.ref())) {
            String s = "ERROR: VCF output file equals input file: " + par.ref();
            Utilities.exit(Par.usage() + s);
        }
        if (vcfOut.equals(par.gt())) {
            String s = "ERROR: VCF output file equals input file: " + par.gt();
            Utilities.exit(Par.usage() + s);
        }
    }
}
