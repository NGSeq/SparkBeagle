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
import org.ngseq.sparkbeagle.ints.IntArray;
import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.vcf.Marker;
import org.ngseq.sparkbeagle.vcf.RefGTRec;
import org.ngseq.sparkbeagle.vcf.SeqCodedRefGT;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * <p>Class {@code SeqCoder3} compresses a sequence of allele-coded
 * {@code RefGTRec} objects.  The class is designed for use with brev v3 format.
 * Compression is performed by storing the list of distinct allele sequences
 * and the allele sequence carried by each haplotype.
 * </p>
 * <p>Class {@code SeqCoder3} is not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SeqCoder3 {

    /**
     * The maximum number of alleles are that permitted in order
     * for the {@code add()} method to return {@code true}
     */
    public static final int MAX_NALLELES = 255;

    /**
     * The major allele frequency threshold for allele coding.
     * Sequence coding should only be applied if the major allele frequency
     * less than or equal to this threshold.
     */
    public static final float COMPRESS_FREQ_THRESHOLD = 0.995f;

    private final Samples samples;
    private final int maxNSeq;
    private final List<RefGTRec> recs;
    private final int[] hap2Seq;
    private final IntList seq2Cnt;  // for processiing allele-coded markers
    private final List<IntList> seq2AlleleMap;

    /**
     * Constructs a new {@code SeqCoder3} for the specified samples.
     * @param samples the list of samples whose data will be compressed
     * @throws NullPointerException if {@code samples == null}
     */
    public SeqCoder3(Samples samples) {
        this(samples, defaultMaxNSeq(samples.nSamples()));
    }

    /**
     * Constructs a new {@code SeqCoder3} for the specified samples.
     * @param samples the list of samples whose data will be compressed
     * @param maxNSeq the maximum number of distinct allele sequences
     * permitted if the {@code add()} method returns {@code true}
     * @throws NullPointerException if {@code samples == null}
     * @throws IllegalArgumentException
     * {@code maxNSeq < 0 || maxNSeq >= Chracter.MAX_VALUE}
     */
    public SeqCoder3(Samples samples, int maxNSeq) {
        if (maxNSeq < 0 || maxNSeq >= Character.MAX_VALUE) {
            throw new IllegalArgumentException(String.valueOf(maxNSeq));
        }
        this.samples = samples;
        this.maxNSeq = maxNSeq;
        this.recs = new ArrayList<>(100);
        this.hap2Seq = new int[2*samples.nSamples()];
        this.seq2Cnt = new IntList(3*maxNSeq/2 + 1);
        this.seq2AlleleMap = new ArrayList<>(maxNSeq + 1);
        initialize();
    }

    /**
     * Returns the default maximum number of sequences for the specified
     * number of samples.  The default value is equal to
     * {@code (int) Math.min((long) Math.pow(2, 2*Math.log10(nSamples) + 1),
     * Character.MAX_VALUE)}
     * @param nSamples the number of samples
     * @return the default maximum number of sequences for the specified
     * number of samples
     * @throws IllegalArgumentException if {@code nSamples < 1}
     */
    public static int defaultMaxNSeq(int nSamples) {
        if (nSamples < 1) {
            throw new IllegalArgumentException(String.valueOf(nSamples));
        }
        if (nSamples==1) {
            return 3;
        }
        else {
            double exponent = 2*Math.log10(nSamples) + 1;
            long maxNSeq = (long) Math.floor(Math.pow(2.0, exponent));
            if (maxNSeq > Character.MAX_VALUE) {
                return Character.MAX_VALUE;
            }
            else {
                return (int) maxNSeq;
            }
        }
    }

    /**
     * Returns the list of samples whose phased genotype data will be compressed.
     * @return the list of samples whose phased genotype data will be compressed
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns the number of compressed {@code RefGTRec} objects.
     * @return the number of compressed {@code RefGTRec} objects
     */
    public int nRecs() {
        return recs.size();
    }

    /**
     * Returns the maximum number of distinct allele sequences.
     * @return the maximum number of distinct allele sequences
     */
    public int maxNSeq() {
        return maxNSeq;
    }

    /**
     * Attempts to add the specified {@code RefGTRec} object to the list of
     * compressed {@code RefGTRec} objects, and returns {@code true}
     * if the {@code RefGTRec} object was added.
     *
     * @param rec reference genotypes for a marker
     * @return {@code true} if the specified {@code RefGTRec} object was
     * added to the list of compressed markers
     *
     * @throws IllegalArgumentException if
     * {@code em.samples().equals(this.samples()) == false}
     * @throws IllegalArgumentException if
     * {@code rec.isAlleleCoded() == false}
     * @throws NullPointerException if {@code em == null}
     */
    public boolean add(RefGTRec rec) {
        if (rec.samples().equals(samples)==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (rec.isAlleleCoded()==false) {
            throw new IllegalArgumentException(rec.getClass().toString());
        }
        boolean success = setAlleleMap(rec);
        if (success) {
            recs.add(rec);
            int majorAllele = rec.majorAllele();
            for (int a=0, n=rec.nAlleles(); a<n; ++a) {
                if (a!=majorAllele) {
                    int nCopies = rec.alleleCount(a);
                    for (int c=0; c<nCopies; ++c) {
                        int h = rec.hapIndex(a, c);
                        int oldSeq = hap2Seq[h];
                        IntList list = seq2AlleleMap.get(oldSeq);
                        int index=0;
                        while (index<list.size() && list.get(index)!=a) {
                            index+=2;
                        }
                        int newSeq = list.get(index+1);
                        if (newSeq != oldSeq) {
                            while (newSeq >= seq2Cnt.size()) {
                                seq2Cnt.add(0);
                            }
                            hap2Seq[h] = newSeq;
                            seq2Cnt.decrementAndGet(oldSeq);
                            seq2Cnt.incrementAndGet(newSeq);
                        }
                    }
                }
            }
        }
        assert seq2Cnt.size()==seq2AlleleMap.size();
        return success;
    }

    private void clearAlleleMap() {
        for (int j=0, n=seq2AlleleMap.size(); j<n; ++j) {
            seq2AlleleMap.get(j).clear();
        }
    }

    private boolean setAlleleMap(RefGTRec rec) {
        assert seq2Cnt.size()==seq2AlleleMap.size();
        int nStartSeq = seq2Cnt.size();
        int[] seq2NonMajorCnt = new int[nStartSeq];
        clearAlleleMap();
        int nAlleles = rec.nAlleles();
        int majorAllele = rec.majorAllele();
        for (int a=0; a<nAlleles; ++a) {
            if (a!=majorAllele) {
                int nCopies = rec.alleleCount(a);
                for (int c=0; c<nCopies; ++c) {
                    int h = rec.hapIndex(a, c);
                    int seq = hap2Seq[h];
                    ++seq2NonMajorCnt[seq];
                    IntList list = seq2AlleleMap.get(seq);
                    if (list.isEmpty()) {
                        list.add(a);
                        list.add(seq);
                    }
                    else {
                        int index=0;
                        while (index<list.size() && list.get(index)!=a) {
                            index+=2;
                        }
                        if (index==list.size()) {
                            list.add(a);
                            list.add(seq2AlleleMap.size());
                            seq2AlleleMap.add(new IntList(4));
                        }
                    }
                }
            }
        }
        addMajorAllele(seq2NonMajorCnt, majorAllele);
        if (seq2AlleleMap.size() >= maxNSeq) {
            seq2AlleleMap.subList(nStartSeq, seq2AlleleMap.size()).clear();
            return false;
        }
        else {
            return true;
        }
    }

    private void addMajorAllele(int[] seq2NonMajorCnt, int majorAllele) {
        for (int s=0; s<seq2NonMajorCnt.length; ++s) {
            if (seq2NonMajorCnt[s] < seq2Cnt.get(s)) {
                IntList list = seq2AlleleMap.get(s);
                if (list.isEmpty()) {
                    list.add(majorAllele);
                    list.add(s);
                }
                else {
                    // assign major allele the existing sequence index
                    list.add(list.get(0));
                    assert list.get(1)==s;
                    list.add(seq2AlleleMap.size());
                    list.set(0, majorAllele);
                    seq2AlleleMap.add(new IntList(4));
                }
            }
        }
    }

    /**
     * Returns and clears the stored list of compressed {@code RefGTRec}
     * objects.
     *
     * @return the list of compressed {@code RefGTRec} objects
     */
    public List<RefGTRec> getCompressedList() {
        if (recs.isEmpty()) {
            return Collections.emptyList();
        }
        List<RefGTRec> list = new ArrayList<>(recs.size());
        int[] seq2Hap = seq2Hap();  // can this be created when adding records?
        IntArray hap2seq =  IntArray.create(hap2Seq, seq2Hap.length);
        for (int j=0, n=recs.size(); j<n; ++j) {
            RefGTRec rec = recs.get(j);
            Marker m = rec.marker();
            IntArray seq2allele = seq2Allele(rec, seq2Hap);
            list.add(new SeqCodedRefGT(m, samples, hap2seq, seq2allele));
        }
        initialize();
        return list;
    }

    private int[] seq2Hap() {
        int[] seqToHap = new int[seq2AlleleMap.size()];
        Arrays.fill(seqToHap, -1);
        for (int h=0; h<hap2Seq.length; ++h) {
            int seq = hap2Seq[h];
            if (seqToHap[seq] == -1) {
                seqToHap[seq] = h;
            }
        }
        return seqToHap;
    }

    private IntArray seq2Allele(RefGTRec rec, int[] seq2Hap) {
        int[] seq2Allele = new int[seq2Hap.length];
        for (int j=0; j<seq2Allele.length; ++j) {
            seq2Allele[j] = rec.get(seq2Hap[j]);
        }
        return IntArray.create(seq2Allele, rec.nAlleles());
    }

    /**
     * Clears the list of compressed {@code RefGTRec} objects.
     */
    private void initialize() {
        recs.clear();
        seq2Cnt.clear();
        seq2AlleleMap.clear();

        // initialize with empty sequence (seq index is 0)
        Arrays.fill(hap2Seq, 0);
        seq2Cnt.add(hap2Seq.length);
        seq2AlleleMap.add(new IntList(4));
    }
}