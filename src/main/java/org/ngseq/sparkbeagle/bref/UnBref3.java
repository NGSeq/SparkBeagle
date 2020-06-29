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

import org.ngseq.sparkbeagle.blbutil.Const;
import org.ngseq.sparkbeagle.blbutil.FileUtil;
import org.ngseq.sparkbeagle.blbutil.SampleFileIt;
import org.ngseq.sparkbeagle.blbutil.Utilities;
import org.ngseq.sparkbeagle.vcf.RefGTRec;
import org.ngseq.sparkbeagle.vcf.VcfWriter;

import java.io.PrintWriter;

/**
 * <p>Class {@code UnBref3} converts files in bref version 3 format
 * into VCF format.
 * </p>
 * <p>Instances of class {@code UnBref3} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class UnBref3 {

    private static final String program = "unbref3.__REV__.jar";

    /**
     * The {@code main()} method is the entry point to the bref program.
     * See the usage() method for usage instructions.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args.length>1) {
            exit(usage());
        }
        if (args.length==1 && args[0].equalsIgnoreCase("help")) {
            exit(usage());
        }
        String fileName = args.length==0 ? null : args[0];
        writeVcf(fileName);
    }

    private static void exit(String msg) {
        System.out.println(usage());
        Utilities.exit(msg);
    }

    private static void writeVcf(String fileName) {
        try (PrintWriter out = FileUtil.stdOutPrintWriter();
                SampleFileIt<RefGTRec> brefIt = bref3It(fileName)) {
            if (brefIt.hasNext()) {
                RefGTRec ve = brefIt.next();
                VcfWriter.writeMetaLinesGT(ve.samples().ids(), program, out);
                out.println(ve.toString());
            }
            while (brefIt.hasNext()) {
                out.println(brefIt.next().toString());
            }
        }
    }

    private static SampleFileIt<RefGTRec> bref3It(String fileName) {
        //TODO: Change to use FSDataInputStream..
        //File file = fileName==null ? null : new File(fileName);
        //return new Bref3It(file);
        return null;
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
        sb.append(" [bref3] > [vcf])");
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("  cat  [bref3]  | java -jar ");
        sb.append(program);
        sb.append(" > [vcf]");
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("where");
        sb.append(Const.nl);
        sb.append("  [bref3]  = the input bref3 file");
        sb.append(Const.nl);
        sb.append("  [vcf]    = the ouput VCF file");
        sb.append(Const.nl);
        return sb.toString();
    }
}
