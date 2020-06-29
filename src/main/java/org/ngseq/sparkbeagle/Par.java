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

import org.apache.commons.cli.CommandLine;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.ngseq.sparkbeagle.beagleutil.ChromInterval;
import org.ngseq.sparkbeagle.blbutil.Const;
import org.ngseq.sparkbeagle.blbutil.Utilities;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * <p>Class {@code Parameters} represents the parameters for a Beagle analysis.
 * </p>
 * <p>Instances of class {@code Parameters} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Par {


    private String header;

    public boolean isBref() {
        return bref;
    }

    private String index;

    private boolean bref = false;
    private String hdfsout;
    private String hdfs;
    // data parameters
    private FileSystem fs = null;

    private  FSDataInputStream gt;
    private  FSDataInputStream gl;
    private  FSDataInputStream gtgl;
    private  FSDataInputStream rheader;
    private  String ref;
    private String gtpath;

    private String out;
    private FSDataInputStream map = null;
    private ChromInterval chromInt;

    // default phasing parameters

    private static final int D_BURNIN = 6;
    private static final int D_ITERATIONS = 12;
    private static final int D_PHASE_STATES = 280;
    private static final float D_PHASE_SEGMENT = 4.0f;

    // default imputation parameters
    private static final boolean D_IMPUTE = true;
    private static final int D_IMP_STATES = 1600;
    private static final float D_IMP_SEGMENT = 6.0f;
    private static final float D_CLUSTER = 0.005f;
    private static final boolean D_AP = false;
    private static final boolean D_GP = false;

    // default general parameters
    private static final int D_NE = 1_000_000;
    private static final float D_ERR = 0.0001f;
    private static final float D_WINDOW = 40.0f;
    private static final float D_OVERLAP = 4.0f;
    private static final int D_SEED = -99999;
    private static final int D_NTHREADS = Integer.MAX_VALUE;
    private static final float D_STEP = 0.1f;
    private static final int D_NSTEPS = 7;


    // phasing parameters
    private int burnin = D_BURNIN;
    private int iterations = D_ITERATIONS;
    private int phase_states = D_PHASE_STATES;
    private float phase_segment = D_PHASE_SEGMENT;

    // imputation parameters
    private boolean impute = true;
    private int imp_states = D_IMP_STATES;
    private float imp_segment = D_IMP_SEGMENT;
    private float cluster = D_CLUSTER;
    private boolean ap= D_AP;
    private boolean gp = D_GP;

    // general parameters
    private float ne = D_NE;
    private float err = D_ERR;
    private float window = D_WINDOW;
    private float overlap = D_OVERLAP;
    private long seed = D_SEED;
    private int nthreads;
    private float step = D_STEP;
    private int nsteps = D_NSTEPS;

    // undocumented parameters
    private FSDataInputStream truth = null;
    private  FSDataInputStream excludesamples = null;
    private  FSDataInputStream excludemarkers = null;


    private  int pol;


    /**
     * Constructs a new {@code Parameters} instance from the specified
     * command line arguments.
     * @param args the command line arguments
     * @throws IllegalArgumentException if a command line argument
     * is incorrectly specified
     * @throws NumberFormatException if a numeric value for a parameter
     * is incorrectly specified
     * @throws NullPointerException if {@code args ==  null}
     */
    public Par(CommandLine cmd) {
        int IMAX = Integer.MAX_VALUE;
        long LMIN = Long.MIN_VALUE;
        long LMAX = Long.MAX_VALUE;
        float FMIN = Float.MIN_VALUE;
        float FMAX = Float.MAX_VALUE;


        hdfs = (cmd.hasOption("hdfs")==true)? (cmd.getOptionValue("hdfs")):"";

        try {
            fs = FileSystem.get(new Configuration());
        } catch (IOException e) {
            e.printStackTrace();
        }

        //mode = cmd.getOptionValue("mode");
        try {
            //header = getRheader(hdfs+cmd.getOptionValue("rheader"));
            header = cmd.getOptionValue("rheader");
            gl = (cmd.hasOption("gl")==true)? fs.open(new Path(hdfs+cmd.getOptionValue("gl"))):null;

            ref = (cmd.hasOption("ref")==true)?  hdfs+cmd.getOptionValue("ref"):null;

            gt = (cmd.hasOption("gt")==true)? fs.open(new Path(hdfs+cmd.getOptionValue("gt"))):null;
            gtpath =  (cmd.hasOption("gt")==true)? hdfs+cmd.getOptionValue("gt"):null;

            gtgl = (cmd.hasOption("gtgl")==true)?  fs.open(new Path(hdfs+cmd.getOptionValue("gtgl"))):null;

            map = (cmd.hasOption("map")==true)? fs.open(new Path(hdfs+cmd.getOptionValue("map"))):null;

            index = (cmd.hasOption("index")==true)?  hdfs+cmd.getOptionValue("index"):null;


            //rheader = (cmd.hasOption("rheader")==true)? fs.open(new Path(cmd.getOptionValue("rheader"))):null;


        } catch (IOException e) {
            e.printStackTrace();
        }

        hdfsout = (cmd.hasOption("out")==true)? (cmd.getOptionValue("out")):"beagle";
        nthreads = (cmd.hasOption("nthreads")==true)? Integer.valueOf(cmd.getOptionValue("nthreads")):1;
        out = (cmd.hasOption("localout")==true)? (cmd.getOptionValue("localout")):"beagle";
        iterations = (cmd.hasOption("niterations")==true)? Integer.valueOf(cmd.getOptionValue("niterations")):0;
        //gprobs = (cmd.hasOption("gprobs")==true)? Boolean.getBoolean(cmd.getOptionValue("gprobs")):true;
        burnin = (cmd.hasOption("burnits")==true)? Integer.valueOf(cmd.getOptionValue("burnits")):0;
        phase_states = (cmd.hasOption("phaseits")==true)? Integer.valueOf(cmd.getOptionValue("phaseits")):0;
        window = (cmd.hasOption("window")==true)? Float.valueOf(cmd.getOptionValue("window")):D_WINDOW;
        overlap = (cmd.hasOption("overlap")==true)? Float.valueOf(cmd.getOptionValue("overlap")):D_OVERLAP;
        step = (cmd.hasOption("step")==true)? Float.valueOf(cmd.getOptionValue("step")):D_STEP;
        chromInt = (cmd.hasOption("chrom")==true)? parseChromInt(cmd.getOptionValue("chrom")):null;
        cluster = (cmd.hasOption("cluster")==true)? Float.valueOf(cmd.getOptionValue("cluster")):D_CLUSTER;
        //nsamples = (cmd.hasOption("nsamples")==true)? Integer.valueOf(cmd.getOptionValue("nsamples"))4;
        //nsamples=4; //Beagle
        pol = (cmd.hasOption("pol")==true)? Integer.valueOf(cmd.getOptionValue("pol")):5000;
        bref = (cmd.hasOption("bref")==true)? true:false;

    }



    /**
     * Returns a description of the Beagle command line arguments.
     * @return a description of the Beagle command line arguments
     */
    public static String usage() {
        String nl = Const.nl;
        return  "Usage: " + Main.COMMAND + " [arguments]" + nl
                + nl
                + "data parameters ..." + nl
                + "  gt=<VCF file: use GT field>                        (optional)" + nl
                + "  ref=<bref3 or VCF file with phased genotypes>      (optional)" + nl
                + "  out=<output file prefix>                           (required)" + nl
//                + "  ped=<linkage format pedigree file>                 (optional)" + nl
                + "  map=<PLINK map file with cM units>                 (optional)" + nl
                + "  chrom=<[chrom] or [chrom]:[start]-[end]>           (optional)" + nl
                + "  excludesamples=<file with 1 sample ID per line>    (optional)" + nl
                + "  excludemarkers=<file with 1 marker ID per line>    (optional)" + nl + nl

                + "phasing parameters ..." + nl
                + "  burnin=<number of burnin iterations>               (default=" + D_BURNIN + ")" + nl
                + "  iterations=<number of phasing iterations>          (default=" + D_ITERATIONS + ")" + nl
                + "  phase-states=<model states for phasing>            (default=" + D_PHASE_STATES + ")" + nl
                + "  phase-segment=<min haplotype segment length (cM)>  (default=" + D_PHASE_SEGMENT + ")" + nl + nl

                + "imputation parameters ..." + nl
                + "  impute=<impute ungenotyped markers (true/false)>   (default=" + D_IMPUTE + ")" + nl
                + "  imp-states=<model states for imputation>           (default=" + D_IMP_STATES + ")" + nl
                + "  imp-segment=<min haplotype segment length (cM)>    (default=" + D_IMP_SEGMENT + ")" + nl
                + "  imp-cluster=<max cM in a marker cluster>           (default=" + D_CLUSTER + ")" + nl
                + "  imp-ap=<print posterior allele probabilities>      (default=" + D_AP + ")" + nl
                + "  imp-gp=<print posterior genotype probabilities>    (default=" + D_GP + ")" + nl + nl

                + "general parameters ..." + nl
                + "  ne=<effective population size>                     (default=" + D_NE + ")" + nl
                + "  err=<allele mismatch rate>                         (default=" + D_ERR + ")" + nl
                + "  window=<window length in cM>                       (default=" + D_WINDOW + ")" + nl
                + "  overlap=<window overlap in cM>                     (default=" + D_OVERLAP + ")" + nl
                + "  seed=<random seed>                                 (default=" + D_SEED + ")" + nl
                + "  nthreads=<number of threads>                       (default: machine-dependent)" + nl
                + "  step=<IBS step length (cM)>                        (default=" + D_STEP + ")" + nl
                + "  nsteps=<number of IBS steps>                       (default=" + D_NSTEPS + ")" + nl + nl;
    }

    private static ChromInterval parseChromInt(String str) {
        ChromInterval ci = ChromInterval.parse(str);
        if (str!=null && str.length()>0 && ci==null) {
            throw new IllegalArgumentException("Invalid chrom parameter: " + str);
        }
        return ci;
    }

    /**
     * Returns the nthreads parameter, which is equal to
     * {@code Runtime.getRuntime().availableProcessors()} if
     * {@code nthreads == Integer.MAX_VALUE}.
     * @return the nthreads parameter
     */
    private static int modNthreads(int nthreads) {
        if (nthreads==Integer.MAX_VALUE) {
            return Runtime.getRuntime().availableProcessors();
        }
        else {
            return nthreads;
        }
    }


    /**
     * Returns the gt parameter or {@code null} if no gt parameter was
     * specified.
     * @return the gt parameter or {@code null} if no gt parameter was
     * specified
     */
    public FSDataInputStream gt() {
        return gt;
    }

    /**
     * Returns the ref parameter or {@code null} if no ref parameter was
     * specified.
     * @return the ref parameter or {@code null} if no ref parameter was
     * specified
     */
    public String ref() {
        return ref;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter
     */
    public String out() {
        return out;
    }

    /**
     * Returns the ped parameter or {@code null}
     * if no ped parameter was specified.
     *
     * @return the ped parameter or {@code null}
     * if no ped parameter was specified
     */
    public FSDataInputStream ped() {
        return null; // ignoring ped file at present
    }

    /**
     * Returns the map parameter.
     * @return the map parameter
     */
    public FSDataInputStream map() {
        return map;
    }

    /**
     * Returns the chromosome interval or {@code null} if no chrom
     * parameter was specified.
     *
     * @return the chromosome interval or {@code null} if no chrom
     * parameter was specified.
     */
    public ChromInterval chromInt() {
        return chromInt;
    }
    /**
     * Returns the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     *
     * @return the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified
     */
    public FSDataInputStream excludesamples() {
        return excludesamples;
    }

    /**
     * Returns the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified.
     *
     * @return the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified
     */
    public FSDataInputStream excludemarkers() {
        return excludemarkers;
    }

    // phasing parameters

    /**
     * Returns the burnin parameter.
     * @return the burnin parameter
     */
    public int burnin() {
        return burnin;
    }

    /**
     * Returns the iterations parameter.
     * @return the iterations parameter
     */
    public int iterations() {
        return iterations;
    }

    /**
     * Returns the phase-states parameter.
     * @return the phase-states parameter
     */
    public int phase_states() {
        return phase_states;
    }

    /**
     * Returns the phase-segment parameter.
     * @return the phase-segment parameter
     */
    public float phase_segment() {
        return phase_segment;
    }

    // imputation parameters

    /**
     * Returns the impute parameter.
     * @return the impute parameter
     */
    public boolean impute() {
        return impute;
    }

    /**
     * Returns the imp-states parameter.
     * @return the imp-states parameter
     */
    public int imp_states() {
        return imp_states;
    }

    /**
     * Returns the imp-segment parameter.
     * @return the imp-segment parameter
     */
    public float imp_segment() {
        return imp_segment;
    }

    /**
     * Returns the cluster parameter.
     * @return the cluster parameter
     */
    public float cluster() {
        return cluster;
    }

    /**
     * Returns the ap parameter.
     * @return the ap parameter
     */
    public boolean ap() {
        return ap;
    }

    /**
     * Returns the gp parameter.
     * @return the gp parameter
     */
    public boolean gp() {
        return gp;
    }

    // general parameters

    /**
     * Returns the ne parameter
     * @return the ne parameter
     */
    public float ne() {
        return ne;
    }

    /**
     * Returns the err parameter.
     * @return the err parameter
     */
    public float err() {
        return err;
    }

    /**
     * Returns the window parameter.
     * @return the window parameter
     */
    public float window() {
        return window;
    }

    /**
     * Return the overlap parameter.
     * @return the overlap parameter.
     */
    public float overlap() {
        return overlap;
    }

    /**
     * Returns the seed parameter.
     * @return the seed parameter
     */
    public long seed() {
        return seed;
    }

    /**
     * Returns the nthreads parameter.
     * @return the nthreads parameter
     */
    public int nthreads() {
        return nthreads;
    }

    /**
     * Returns the step parameter.
     * @return the step parameter
     */
    public float step() {
        return step;
    }

    /**
     * Returns the nsteps parameter.
     * @return the nsteps parameter
     */
    public int nsteps() {
        return nsteps;
    }

    // undocumented parameters

    /**
     * Returns the truth parameter
     * @return the truth
     */
    public FSDataInputStream truth() {
        return truth;
    }

    public String gtpath() {
        return gtpath;
    }

    public String hdfs() {
        return hdfs;
    }
    public String hdfsout() {
        return hdfsout;
    }

    public int getPartitionOverlap() {
        return pol;
    }

    public String getRefPath() {
        return ref;
    }

    public String getIndexPath() {
        return index;
    }
    public String getHeaderPath() {
        return hdfs+header;
    }
    public String getHeader() {
        String h ="";
        try {

            InputStreamReader isr = new InputStreamReader(fs.open(new Path(getHeaderPath())));
            BufferedReader br = new BufferedReader(isr);
            String line;
            while ((line = br.readLine()) != null) {
                h += line;
                h += "\n";
            }
            if(!h.startsWith("##")) //Somehow ## is lost from start of header
                h="##"+h;
            System.out.print(h);
        }
        catch(IOException e) {
            Utilities.exit("Error reading header", e);
        }
        return h;

    }
}
