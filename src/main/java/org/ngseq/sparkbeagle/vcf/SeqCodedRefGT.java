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
package org.ngseq.sparkbeagle.vcf;

import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.ints.IntArray;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code SeqCodedRefGT}  represents phased, non-missing
 * genotypes for a list of reference samples at a single marker.
 * Genotype emission probabilities are determined by the sample
 * genotypes.
 * </p>
 * <p>Instances of class {@code SeqCodedRefGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SeqCodedRefGT implements RefGTRec {

    private final Marker marker;
    private final Samples samples;
    private final IntArray hapToSeq;
    private final IntArray seqToAllele;

    /**
     * Creates a new {@code SeqCodedRefGT} instance with phased,
     * non-missing genotypes from the specified marker, samples,
     * and haplotype alleles.  The contract for the constructed object
     * is undefined if any element of {@code hapToSeq} is negative or
     * greater than or equal to {@code seqToAllele.size()} or if any element
     * of {@code seqToAllele} is negative or greater than or equal to
     * {@code marker.nAlleles()}.
     *
     * @param marker the marker
     * @param samples the samples
     * @param hapToSeq an array whose {@code j}-th element is the index
     * of the distinct allele sequence carried by the {@code j}-th haplotype
     * @param seqToAllele an array whose {@code j}-th element is the marker
     * allele carried by the {@code j}-th distinct allele sequence
     *
     * @throws IllegalArgumentException if
     * {@code hapToSeq.size() != 2*samples.nSamples()}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public SeqCodedRefGT(Marker marker, Samples samples, IntArray hapToSeq,
        IntArray seqToAllele) {
        if (hapToSeq.size() != 2*samples.nSamples()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        this.marker = marker;
        this.samples = samples;
        this.hapToSeq = hapToSeq;
        this.seqToAllele = seqToAllele;
    }

    @Override
    public boolean isPhased(int sample) {
        if (sample < 0 || sample >= this.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        return true;
    }

    /**
     * Returns {@code true}.
     * @return {@code true}
     */
    @Override
    public boolean isPhased() {
        return true;
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int size() {
        return hapToSeq.size();
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public boolean isGTData() {
        return true;
    }

    @Override
    public float gl(int sample, int allele1, int allele2) {
        boolean match = allele1 == allele1(sample)
                && allele2 == allele2(sample);
        return match ? 1f : 0f;
    }

    @Override
    public int[][] hapIndices() {
        int[] alCnts = new int[marker().nAlleles()];
        for (int h=0, n=size(); h<n; ++h) {
            ++alCnts[get(h)];
        }
        int majAllele = 0;
        for (int al=1; al<alCnts.length; ++al) {
            if (alCnts[al] > alCnts[majAllele]) {
                majAllele = al;
            }
        }
        int[][] hapIndices = new int[alCnts.length][];
        for (int al=0; al<alCnts.length; ++al) {
            if (al!=majAllele) {
                hapIndices[al] = new int[alCnts[al]];
            }
        }
        Arrays.fill(alCnts, 0);
        for (int h=0, n=size(); h<n; ++h) {
            int al = get(h);
            if (al!=majAllele) {
                hapIndices[al][alCnts[al]++] = h;
            }
        }
        return hapIndices;
    }

    @Override
    public boolean isAlleleCoded() {
        return false;
    }

    @Override
    public int majorAllele() {
        int[][] hapIndices = hapIndices();
        int majAllele = -1;
        for (int j=0; j<hapIndices.length && majAllele==-1; ++j) {
            if (hapIndices[j]==null) {
                majAllele = j;
            }
        }
        return majAllele;
    }

    @Override
    public int alleleCount(int allele) {
        int[][] hapIndices = hapIndices();
        if (hapIndices[allele]==null) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return hapIndices[allele].length;
        }
    }

    @Override
    public int allele1(int sample) {
        return seqToAllele.get(hapToSeq.get(2*sample));
    }

    @Override
    public int allele2(int sample) {
        return seqToAllele.get(hapToSeq.get(2*sample + 1));
    }

    @Override
    public int get(int hap) {
        return seqToAllele.get(hapToSeq.get(hap));
    }

    @Override
    public int[] alleles() {
        return IntStream.range(0, hapToSeq.size())
                .map(h -> get(h))
                .toArray();
    }

    @Override
    public int nAlleles() {
        return this.marker().nAlleles();
    }

    @Override
    public int hapIndex(int allele, int copy) {
        int[][] hapIndices = hapIndices();
        if (hapIndices[allele]==null) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return hapIndices[allele][copy];
        }
    }

    @Override
    public boolean isCarrier(int allele, int hap) {
        return get(hap)==allele;
    }

    /**
     * Returns the data represented by {@code this} as a VCF
     * record with a GT format field. The returned VCF record
     * will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return the data represented by {@code this} as a VCF
     * record with a GT format field
     */
    @Override
    public String toString() {
        return GTRec.toVcfRec(this);
    }

    @Override
    public int nMaps() {
        return 2;
    }

    @Override
    public IntArray[] maps() {
        return new IntArray[] {hapToSeq, seqToAllele};
    }

    @Override
    public IntArray map(int index) {
        if (index==0) {
            return hapToSeq;
        }
        else if (index==1) {
            return seqToAllele;
        }
        else {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
    }
}
