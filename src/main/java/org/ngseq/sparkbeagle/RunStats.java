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

import org.ngseq.sparkbeagle.blbutil.Const;
import org.ngseq.sparkbeagle.blbutil.FileUtil;
import org.ngseq.sparkbeagle.blbutil.Utilities;
import org.ngseq.sparkbeagle.vcf.Data;
import org.ngseq.sparkbeagle.vcf.Marker;
import org.ngseq.sparkbeagle.vcf.Markers;

import java.io.File;
import java.io.PrintWriter;

/**
 * Class {@code RunStats} contains methods for storing and printing
 * statistics describing a Beagle analysis.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RunStats {

    private final Par par;
    private final PrintWriter log;
    private final long startNanos;

    private long buildNanos = 0;
    private long totalBuildNanos = 0;

    private long totalIterationNanos = 0;

    private long imputeNanos = 0;
    private long totalImputeNanos = 0;

    /**
     * Constructs a new {@code RunStats} instance.
     * @param par the analysis parameters
     * @throws NullPointerException if {@code par == null}
     */
    RunStats(Par par) {
        this.startNanos = System.nanoTime();
        this.par = par;
        this.log = FileUtil.printWriter(new File(par.out() + ".log"));
    }

    /**
     * Prints initial information about the analysis to a log
     * file and to standard output.
     */
    public void printStartInfo() {
        Utilities.duoPrint(log, Main.SHORT_HELP + Const.nl);
        Utilities.duoPrintln(log, "Start time: " + Utilities.timeStamp());
        //Utilities.duoPrint(log, Utilities.commandLine(Main.PROGRAM, par.args()));
        if (par.ped() != null) {
            String s = Const.nl + "WARNING: This version will not model"
                    + " duos or trios in the pedigree file";
            Utilities.duoPrintln(log, s);
        }
        if (par.map() == null) {
            String s = Const.nl + "No genetic map is specified: using 1 cM = 1 Mb";
            Utilities.duoPrintln(log, s);
        }
        log.flush();
    }

   /**
     * Prints information about the samples to a log
     * file and to standard output.
     * @param data the input data
     */
    public void printSampleSummary(Data data) {
        Pedigree fam = data.ped();
        Utilities.duoPrint(log, Const.nl);
        Utilities.duoPrint(log, String.format("Reference samples: %,11d%n",
                data.nRefSamples()));
        Utilities.duoPrint(log, String.format("Study samples:     %,11d%n",
                data.nTargetSamples()));
        if (par.ped() != null) {
            Utilities.duoPrint(log, "  ");
            Utilities.duoPrint(log, String.valueOf(fam.nSingles()));
            Utilities.duoPrintln(log, " singles");
            Utilities.duoPrint(log, "  ");
            Utilities.duoPrint(log, String.valueOf(fam.nDuos()));
            Utilities.duoPrintln(log, " duos");
            Utilities.duoPrint(log, "  ");
            Utilities.duoPrint(log, String.valueOf(fam.nTrios()));
            Utilities.duoPrintln(log, " trios");
        }
        log.flush();
    }

   /**
     * Prints information about the marker window to a log
     * file and to standard output.
     * @param data the input genotype data
     */
    public void printWindowUpdate(Data data) {
        Markers markers = data.markers();
        Marker first = markers.marker(0);
        Marker last = markers.marker(markers.nMarkers() - 1);
        StringBuilder sb = new StringBuilder(30);
        sb.append(Const.nl);
        sb.append("Window ");
        sb.append(data.windowIndex());
        sb.append(" (");
        String chr = first.chrom();
        if (chr.equals(Const.MISSING_DATA_STRING)==false) {
            sb.append(chr);
            sb.append(Const.colon);
        }
        sb.append(first.pos());
        sb.append(Const.hyphen);
        if (chr.equals(last.chrom())==false) {
            sb.append(last.chrom());
            sb.append(Const.colon);
        }
        sb.append(last.pos());
        sb.append(')');
        sb.append(Const.nl);
        if (data.nRefSamples()>0) {
            sb.append(String.format("Reference markers: %,11d%n", data.nMarkers()));
        }
        sb.append(String.format("Study markers:     %,11d%n", data.nTargetMarkers()));
        Utilities.duoPrint(log, sb.toString());
        log.flush();
    }

    /**
     * Prints information about the complete analysis to a log
     * file and to standard output, and closes the log file.
     * @param nTargetMarkers the total number of target markers analyzed
     * @param nMarkers the total number of markers analyzed
     */
    public void printSummaryAndClose(int nTargetMarkers, int nMarkers) {
        long totalTime = System.nanoTime() - startNanos;
        Utilities.duoPrint(log, Const.nl);
        Utilities.duoPrintln(log, "Cumulative Statistics:" + Const.nl);
        if (nTargetMarkers != nMarkers) {
            Utilities.duoPrint(log,
                    String.format("Reference markers: %,11d%n", nMarkers));
        }
        Utilities.duoPrint(log,
                String.format("Study markers:     %,11d%n%n", nTargetMarkers));

        if (totalBuildNanos > 0) {
            duoPrintNanos("Model building time:           ", totalBuildNanos);
        }
        if (totalIterationNanos > 1000) {
            duoPrintNanos("Haplotype phasing time:        ", totalIterationNanos);
        }
        if (totalImputeNanos > 0) {
            duoPrintNanos("Imputation time:               ", imputeNanos);
        }
        duoPrintNanos(    "Total time:                    ", totalTime);
        Utilities.duoPrintln(log, Const.nl + "End time: "
                + Utilities.timeStamp());
        Utilities.duoPrintln(log, Main.PROGRAM + " finished");
        log.close();
    }

    /**
     * Increases the cumulative time to build the DAG models by the
     * specified number of nanoseconds.
     * @param nanos the nanoseconds required to build an instance
     * of the DAG model
     */
    public void buildNanos(long nanos) {
        buildNanos = nanos;
        totalBuildNanos += nanos;
    }

    /**
     * Stores the time for updating the haplotype estimates and increases the
     * cumulative phasing time by the specified number of nanoseconds.
     * @param nanos the nanoseconds required to updating the haplotype
     * estimates
     */
    public void iterationNanos(long nanos) {
        totalIterationNanos += nanos;
    }

    /**
     * Stores the time for imputing ungenotyped marker and increases
     * the cumulative imputation time by the specified number
     * of nanoseconds.
     * @param nanos the nanoseconds required to impute ungenotyped
     * markers
     */
    public void imputationNanos(long nanos) {
        imputeNanos = nanos;
        totalImputeNanos += nanos;
    }

    /**
     * Prints run time for most recent imputation to a log file
     * and to standard output.
     */
    public void printImputationUpdate() {
        Utilities.duoPrint(log, Const.nl);
        duoPrintNanos("Imputation time:               ", imputeNanos);
        log.flush();
    }

    /**
     * Prints the specified string to the log file and to standard out.
     * @param msg the message to be printed
     */
    public void println(String msg) {
        Utilities.duoPrintln(log, msg);
        log.flush();
    }

    /**
     * Prints information about the specified iteration.
     * @param window the window
     * @param iter the iteration
     */
    public void printIterationUpdate(int window, int iter) {
        Utilities.duoPrint(log, Const.nl + "Window=" + window
                + " Iteration=" + iter + Const.nl);
        duoPrintNanos("Time for building model:         ", buildNanos);
        log.flush();
    }

    public void printPhasingIterationUpdate(int it, boolean burnin, long elapsedNanos) {
        iterationNanos(elapsedNanos);
        String msg;
        if (burnin) {
            if (it==1) {
               println("");
            }
            msg = "Burnin  iteration " + it + ":";
        }
        else {
            it -= par.burnin();
            if (it==1) {
                println("");
            }
            msg = "Phasing iteration " + it + ":";
        }
        duoPrintNanos(String.format("%1$-31s", msg), elapsedNanos);
    }

    /**
     * Print the specified message followed by the human
     * elapsed time as formatted by
     * {@code blbutil.Utilities.elapsedNanos(nanos)}
     * @param message the message to be printed
     * @param nanos the elapsed time in nanoseconds
     */
    public void duoPrintNanos(String message, long nanos) {
        Utilities.duoPrint(log, message);
        Utilities.duoPrintln(log, Utilities.elapsedNanos(nanos));
    }
}
