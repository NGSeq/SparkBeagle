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
import org.ngseq.sparkbeagle.vcf.RefGTRec;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>Class {@code CompressBref3Writer} writes phased, non-missing genotypes
 * to a binary reference format v3 (bref) file.
 * The {@code close()} method must be called after the last invocation of
 * the {@code write()} method in order to ensure that all buffered
 * data is written to the output binary reference file.
 * </p>
 * <p>Instances of class {@code CompressBrerf3Writer} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class CompressBref3Writer implements BrefWriter {

    private final int maxNAlleles;
    private final List<RefGTRec> buffer;
    private final SeqCoder3 seqCoder;
    private final AsIsBref3Writer brefWriter;
    private final int nonMajThreshold;

    /**
     * Constructs a new {@code CompressBref3fWriter} for the specified data.
     * The Java virtual machine will exit with an error message if an I/O
     * error occurs.
     *
     * @param program the name of the program which is creating the
     * binary reference file.
     * @param samples the list of samples whose genotype data will
     * be written in binary reference format
     * @param maxNSeq the maximum number of distinct allele sequences
     * in a compressed block
     * @param brefFile name of the output binary reference file or
     * {@code null} if the output should be directed to standard output
     * @throws IllegalArgumentException
     * {@code maxNSeq < 0 || maxNSeq >= Chracter.MAX_VALUE}
     * @throws NullPointerException if {@code program == null || samples == null}
     */
    public CompressBref3Writer(String program, Samples samples, int maxNSeq,
            File brefFile) {
        this.maxNAlleles = SeqCoder3.MAX_NALLELES;
        this.buffer = new ArrayList<>(500);
        this.seqCoder = new SeqCoder3(samples, maxNSeq);
        this.brefWriter = new AsIsBref3Writer(program, samples, brefFile);
        this.nonMajThreshold = (maxNSeq/4) + 1;
    }

    @Override
    public Samples samples() {
        return brefWriter.samples();
    }

    @Override
    public void write(RefGTRec rec) {
        if (rec.isAlleleCoded()==false) {
            rec = RefGTRec.alleleCodedInstance(rec);
        }
        if (buffer.size()==Integer.MAX_VALUE) {
            flushBuffer();
        }
        if (convertToSeqCoding(rec)) {
            boolean success = seqCoder.add(rec);
            if (success == false) {
                flushBuffer();
                success = seqCoder.add(rec);
                assert success;
            }
            buffer.add(null);
        }
        else {
            buffer.add(rec);
        }
    }

    private boolean convertToSeqCoding(RefGTRec rec) {
        assert rec.isAlleleCoded();
        if (rec.marker().nAlleles() > maxNAlleles) {
            return false;
        }
        int majAllele = rec.majorAllele();
        int nonMajCnt = 0;
        for (int a=0, n=rec.nAlleles(); a<n; ++a) {
            if (a!=majAllele) {
                nonMajCnt += rec.alleleCount(a);
            }
        }
        return nonMajCnt >= nonMajThreshold;
    }

    private void flushBuffer() {
        List<RefGTRec> list = seqCoder.getCompressedList();
        int index = 0;
        for (int j=0, n=buffer.size(); j<n; ++j) {
            RefGTRec rec = buffer.get(j);
            brefWriter.write( rec==null ? list.get(index++) : rec );
        }
        assert index==list.size();
        buffer.clear();
    }

    @Override
    public void close() {
        flushBuffer();
        brefWriter.close();
    }
}
