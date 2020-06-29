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
package org.ngseq.sparkbeagle.bref;

import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.blbutil.Const;
import org.ngseq.sparkbeagle.blbutil.FileIt;
import org.ngseq.sparkbeagle.blbutil.SampleFileIt;
import org.ngseq.sparkbeagle.blbutil.Utilities;
import org.ngseq.sparkbeagle.vcf.RefGTRec;
import org.ngseq.sparkbeagle.vcf.RefIt;

import java.io.File;

/**
 * <p>Class {@code Bref3} converts files in VCF format into 
 * bref version 3 format.
 * </p>
 * <p>Instances of class {@code Bref3} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Bref3 {

    private static final String program = "bref3.__REV__.jar";

    /**
     * The {@code main()} method is the entry point to the bref program.
     * See the usage() method for usage instructions.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args.length>2) {
            System.out.println(usage());
            System.exit(0);
        }
        if (args.length==1 && args[0].equalsIgnoreCase("help")) {
            System.out.println(usage());
            System.exit(0);
        }
        boolean useStdIn = useStdIn(args);
        int maxNSeq = -1;
        if (args.length==2 || (useStdIn && args.length==1)) {
            maxNSeq = maxNSeq(args[args.length-1]);
        }
        String inputFile = useStdIn ? null : args[0];
        writeBref(inputFile, maxNSeq);
    }

    private static boolean useStdIn(String[] sa) {
        if (sa.length==0) {
            return true;
        }
        else if (sa.length==1) {
            return (sa[0].endsWith(".vcf") || sa[0].endsWith(".vcf.gz"))==false;
        }
        else {
            return false;
        }
    }

    private static int maxNSeq(String arg) {
        int maxNSeq = -1;
        try {
            maxNSeq = Integer.parseInt(arg);
            if (maxNSeq < 1 || maxNSeq > Character.MAX_VALUE) {
                exit("Error: invalid <nSeq> " + arg);
            }
        }
        catch (NumberFormatException e) {
            exit("Error: <nSeq> is not a parsable integer: " + arg);
        }
        return maxNSeq;
    }

    private static void exit(String msg) {
        System.out.println(usage());
        Utilities.exit(msg);
    }

    private static void writeBref(String fileName, int maxNSeq) {
        try (SampleFileIt<RefGTRec> it = refIt(fileName);
                BrefWriter brefOut = brefOut(it.samples(), maxNSeq)) {
            while (it.hasNext()) {
                brefOut.write(it.next());
            }
        }
    }

    private static SampleFileIt<RefGTRec> refIt(String fileName) {
        FileIt<String> it = null;
        //TODO: modify to use FSDataInputStream
        /*if (fileName==null) {
            it = InputIt.fromStdIn();
        }
        else if (fileName.endsWith(".gz")) {
            it = InputIt.fromGzipFile(new File(fileName));
        }
        else {
            it = InputIt.fromTextFile(new File(fileName));
        }*/
        return RefIt.create(it);
    }

    private static BrefWriter brefOut(Samples samples, int maxNSeq) {
        File outFile = null;    // write to standard output
        if (maxNSeq<0) {
            return new AsIsBref3Writer(program, samples, outFile);
        }
        else {
            return new CompressBref3Writer(program, samples, maxNSeq, outFile);
        }
    }

    private static String usage() {
        StringBuilder sb = new StringBuilder(500);
        sb.append("usage:");
        sb.append(Const.nl);
        sb.append("  java -jar ");
        sb.append(program);
        sb.append(" help");
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("  java -jar ");
        sb.append(program);
        sb.append(" [vcf] <nseq>  > [bref3]");
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("  cat   [vcf]   | java -jar ");
        sb.append(program);
        sb.append(" <nseq>  > [bref3]");
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("where");
        sb.append(Const.nl);
        sb.append("  [bref3]  = the output bref3 file");
        sb.append(Const.nl);        
        sb.append("  [vcf]    = A VCF file with phased, non-missing genotype data.  If the");
        sb.append(Const.nl);
        sb.append("             file is gzip-compressed, its filename must end in \".gz\"");
        sb.append(Const.nl);
        sb.append("             and \"cat\" must be replaced with \"zcat\"");
        sb.append(Const.nl);
        sb.append("  <nseq>   = optional argument for maximum number of unique sequences");
        sb.append(Const.nl);
        sb.append("             in a bref3 block. If there are N reference samples,");
        sb.append(Const.nl);
        sb.append("             the default value is: <max-seq>=2^(2*log10(N) + 1)");
        sb.append(Const.nl);
        return sb.toString();
    }
}
