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

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.blbutil.BGZIPOutputStream;
import org.ngseq.sparkbeagle.blbutil.FileUtil;
import org.ngseq.sparkbeagle.blbutil.Utilities;
import org.ngseq.sparkbeagle.imp.ImpData;
import org.ngseq.sparkbeagle.imp.ImputedVcfWriter;
import org.ngseq.sparkbeagle.imp.StateProbs;
import org.ngseq.sparkbeagle.vcf.GT;
import org.ngseq.sparkbeagle.vcf.VcfWriter;

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReferenceArray;

/**
 * <p>Class {@code WindowWriter} writes VCF and IBD output data.
 * </p>
 * <p>Instances of class {@code WindowWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class WindowWriter implements Closeable {

    private final Samples samples;
    //private final String outPrefix;
    private final File vcfOutFile = null;
    private final String hdfs;
    private final String hdfsPath;
    private FileSystem fs;


    /**
     * Constructs a new {@code WindowWriter} object.
     * @param samples the sample whose data will be printed
     * @param outPrefix the output file prefix
     *
     * @throws IllegalArgumentException if {@code outPrefix.length() == 0}
     * @throws NullPointerException if
     * {@code samples == null || outPrefix == null}
     */
    public WindowWriter(Samples samples, String hdfsPath, String hdfs) {
        if (samples==null) {
            throw new NullPointerException("samples==null");
        }

        this.samples = samples;
        this.hdfs=hdfs;
        this.hdfsPath = hdfsPath;

        try {
            fs = FileSystem.get( new URI(this.hdfs), new Configuration() );
        } catch (URISyntaxException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Returns the samples whose data is written by {@code this}.
     * @return the samples whose data is written by {@code this}
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Prints VCF records with GT and GP format fields for markers with
     * index between {@code cd.lastSplice()} (inclusive) and
     * {@code cd.nextSplice()} (exclusive).
     *
     * @param cd the input data for the current marker window
     * @param gv scaled genotype probabilities for the target samples
     *
     * @throws NullPointerException if {@code cd == null || gv == null}
     */
    public void printGV(CurrentData cd, GenotypeValues gv) {
        boolean append = true;
        try (PrintWriter vcfOut = FileUtil.bgzipPrintWriter(vcfOutFile, append)) {
            VcfWriter.appendRecords(gv, cd.prevTargetSpliceStart(),
                    cd.nextTargetSpliceStart(), vcfOut);
        }
    }

    /**
     * Prints the data in {@code alProbs} for markers
     * with index between {@code refStart} (inclusive) and
     * {@code refEnd} (exclusive) to the output
     * VCF file: {@code this.outPrefix() + ".vcf.gz"}.
     *
     * @param impData the input data for genotype imputation
     * @param stateProbs the imputed state probabilities
     * @param refStart the starting ref marker index (inclusive)
     * @param refEnd the ending ref marker index (exclusive)
     *
     * @throws IllegalArgumentException if
     * {@code stateProbs.size() != impData.nTargHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code refStart < 0 || refEnd > impData.refGT().nMarkers()}
     * @throws NullPointerException if {@code impData==null || stateProbs==null}
     * @throws NullPointerException if any element of {@code stateProbs} is
     * {@code null}
     */
    public void print(ImpData impData, AtomicReferenceArray<StateProbs> stateProbs,
            int refStart, int refEnd) {
        if (stateProbs.length() != impData.nTargHaps()) {
            throw new IllegalArgumentException("inconsistent data:");
        }
        int nThreads = impData.par().nthreads();
        final AtomicInteger counter = new AtomicInteger(0);
        final AtomicReferenceArray<byte[]> output
                = new AtomicReferenceArray<>(impData.nClusters());

        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(new RunImputeOutput(impData, stateProbs, refStart, refEnd,
                    counter, output));
        }
        try {
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable e) {
            Utilities.exit("ERROR", e);
        }
        print(output, vcfOutFile);

    }

    public void print(ImpData impData, AtomicReferenceArray<StateProbs> stateProbs,
                      int refStart, int refEnd, int startpos, int endpos ) {
        if (stateProbs.length() != impData.nTargHaps()) {
            throw new IllegalArgumentException("inconsistent data:");
        }
        final AtomicInteger counter = new AtomicInteger(0);

        int nClusters = impData.nClusters();

        FSDataOutputStream fos = null;
        try {
            NumberFormat nf = new DecimalFormat("#00000000000000");
            fos = fs.create(new Path(hdfs+"/"+hdfsPath+"/"+ nf.format(startpos) +"_"+endpos+".gz"));
        } catch (IOException e) {
            e.printStackTrace();
        }

        int cluster = counter.getAndIncrement();
        while (cluster < nClusters) {
                ImputedVcfWriter ivw = new ImputedVcfWriter(impData, cluster, refStart, refEnd);
                ByteArrayOutputStream baos = new ByteArrayOutputStream();
                try (PrintWriter out=new PrintWriter(new BGZIPOutputStream(baos, false))) {
                    ivw.appendRecords(stateProbs, out);
                }
                try {
                    fos.write(baos.toByteArray());
                } catch (IOException e) {
                    e.printStackTrace();
                }
                cluster = counter.getAndIncrement();
        }
        try {
            fos.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Writes the data in phased genotypes for the specified markers
     * to the output VCF file: {@code this.outPrefix() + ".vcf.gz"}.
     *
     * @param phasedTarg the estimated target haplotypes
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param nThreads the number of parallel threads to use
     *
     * @throws IllegalArgumentException if
     * {@code isImputed.length != alProbs.nMarkers()}
     * @throws IllegalArgumentException if {@code phasedTarg.isPhased() == false}
     * @throws IllegalArgumentException if {@code nThreads < 1}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end > phasedTarg.nMarkers() || start > end}
     */
    public void print(GT phasedTarg, int start, int end, int nThreads) {
        int step = nMarkersPerStep(phasedTarg.nSamples());
        int nSteps = nSteps(end-start, step);
        final AtomicReferenceArray<byte[]> output
                = new AtomicReferenceArray<>(nSteps);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(new RunPhaseOutput(phasedTarg, start, end, step, nSteps,
                    output)) ;
        }
        try {
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable t) {
            Utilities.exit(t);
        }
        print(output, vcfOutFile);
    }

    private static void print(AtomicReferenceArray<byte[]> output,
            File outFile)  {
        boolean append = true;
        try {
            try (OutputStream fos = new BufferedOutputStream(
                    new FileOutputStream(outFile, append))) {
                for (int j=0, n=output.length(); j<n; ++j) {
                    fos.write(output.get(j));
                }
            }
        } catch (IOException e) {
            fileOutputError(outFile, e);
        }
    }

    private static int nMarkersPerStep(int nSamples) {
        int nBytesPerStep = 5*(1 << 20);
        int bytesPerSample = 4;
        int nMarkersPerStep = nBytesPerStep / (nSamples*bytesPerSample);
        return (nMarkersPerStep == 0) ? 1 : (nSamples*bytesPerSample);
    }

    private static int nSteps(int n, int step) {
        int nSteps = n / step;
        if (nSteps * step != n) {
            ++nSteps;
        }
        return nSteps;
    }

    @Override
    public void close() {
        boolean append = true;
        try {
            try (FileOutputStream fos = new FileOutputStream(vcfOutFile, append);
                    BufferedOutputStream bos = new BufferedOutputStream(fos);
                    BGZIPOutputStream bgzip = new BGZIPOutputStream(bos, true)) {
                // write empty BGZIP block to bgzip by closing bgzip
            }
        } catch (IOException e) {
            Utilities.exit("Error closing file: " + vcfOutFile, e);
        }
    }

    private static void fileOutputError(File file, Exception e) {
        Utilities.exit("Error writing to file: " + file, e);
    }

    private static class RunImputeOutput implements Runnable {

        private final ImpData impData;
        private final int refStart;
        private final int refEnd;
        private final int nClusters;
        private final AtomicInteger counter;
        private final AtomicReferenceArray<StateProbs> stateProbs;
        private final AtomicReferenceArray<byte[]> output;

        public RunImputeOutput(ImpData impData,
                AtomicReferenceArray<StateProbs> stateProbs, int refStart, int refEnd,
                AtomicInteger clusterIndex, AtomicReferenceArray<byte[]> output) {
            this.impData = impData;
            this.refStart = refStart;
            this.refEnd = refEnd;
            this.nClusters = impData.nClusters();
            this.counter = clusterIndex;
            this.stateProbs = stateProbs;
            this.output = output;

        }

        @Override
        public void run() {
            try {
                int cluster = counter.getAndIncrement();
                while (cluster < nClusters) {
                    ImputedVcfWriter ivw = new ImputedVcfWriter(
                            impData, cluster, refStart, refEnd);
                    ByteArrayOutputStream baos = new ByteArrayOutputStream();
                    try (PrintWriter out=new PrintWriter(
                            new BGZIPOutputStream(baos, false))) {
                        ivw.appendRecords(stateProbs, out);
                    }
                    output.set(cluster, baos.toByteArray());
                    cluster = counter.getAndIncrement();
                }
            }
            catch(Throwable t) {
                Utilities.exit(t);
            }
        }
    }

    private static class RunPhaseOutput implements Runnable {

        private final GT phasedTarg;
        private final int start;
        private final int end;
        private final AtomicInteger counter;
        private final int step;
        private final int nSteps;
        private final AtomicReferenceArray<byte[]> output;

        public RunPhaseOutput(GT phasedTarg, int start, int end,
                int step, int nSteps, AtomicReferenceArray<byte[]> output) {
            this.phasedTarg = phasedTarg;
            this.start = start;
            this.end = end;
            this.step = step;
            this.nSteps = nSteps;
            this.counter = new AtomicInteger(0);
            this.output = output;
        }

        @Override
        public void run() {
            try {
                boolean printDS = false;
                boolean printGP = false;
                int index = counter.getAndIncrement();
                while (index < nSteps) {
                    int segStart = start + step*index;
                    int segEnd = Math.min(segStart + step, end);
                    ByteArrayOutputStream baos = new ByteArrayOutputStream();
                    try (PrintWriter vcfOut=new PrintWriter(
                            new BGZIPOutputStream(baos, false))) {
                        VcfWriter.appendRecords(phasedTarg, segStart,
                                segEnd, printDS, printGP, vcfOut);
                    }
                    output.set(index, baos.toByteArray());
                    index = counter.getAndIncrement();
                }
            }
            catch(Throwable t) {
                Utilities.exit(t);
            }
        }
    }
}
