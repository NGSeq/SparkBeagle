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
package org.ngseq.sparkbeagle.haplotype;

import org.ngseq.sparkbeagle.vcf.Marker;
import org.ngseq.sparkbeagle.vcf.Markers;

/**
 * <p>Class {@code BitHapPair} represents a pair of haplotypes for a sample.
 * The class stores alleles as a bit array.
 * </p>
 * Instances of class {@code BitHapPair} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BitHapPair implements HapPair {

    private static final int LOG2_BITS_PER_WORD = 6;

    private final Markers markers;
    private final int idIndex;
    private final long[] bits1;
    private final long[] bits2;

    /**
     * Constructs a new {@code BitHapPair} instance.
     * @param markers the sequence of markers
     * @param idIndex the sample identifier index
     * @param alleles1 the sequence of allele indices for the first haplotype
     * @param alleles2 the sequence of alleles indices for the second haplotype
     *
     * @throws IllegalArgumentException if
     * {@code alleles1.length != markers.nMarkers()
     * || alleles2.length != markers.nMarkers()}
     * @throws IllegalArgumentException if {@code alleles1[k] < 0 ||
     * allele1[k] >= markers.marker(k).nAlleles()} for some {@code k} satisfying
     * {@code 0 <= k && k < markers.nMarkers()}
     * @throws IllegalArgumentException if {@code alleles2[k] < 0 ||
     * allele2[k] >= markers.marker(k).nAlleles()} for some {@code k} satisfying
     * {@code 0 <= k && k < markers.nMarkers()}
     * @throws IndexOutOfBoundsException if {@code idIndex < 0}
     * @throws NullPointerException if
     * {@code marker == null || alleles1 == null || allele2 == null}
     */
    public BitHapPair(Markers markers, int idIndex, int[] alleles1,
            int[] alleles2) {
        if (alleles1.length != markers.nMarkers()
                || alleles2.length != markers.nMarkers()) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (idIndex < 0) {
            throw new IndexOutOfBoundsException(String.valueOf(idIndex));
        }
        this.markers = markers;
        this.idIndex = idIndex;
        this.bits1 = toBitArray(markers, alleles1);
        this.bits2 = toBitArray(markers, alleles2);
    }

    private static long[] toBitArray(Markers markers, int[] alleles) {
        int nWords = (markers.sumHaplotypeBits() + (Long.SIZE-1)) >> LOG2_BITS_PER_WORD;
        long[] bits = new long[nWords];
        int bitIndex = 0;
        for (int k=0; k<alleles.length; ++k) {
            int allele = alleles[k];
            if (allele < 0 || allele >= markers.marker(k).nAlleles()) {
                String s = "allele \"" + allele + "\" out of bounds for marker: "
                        + markers.marker(k);
                throw new IllegalArgumentException(s);
            }
            int mask = 1;
            int nBits = markers.sumHaplotypeBits(k+1) - markers.sumHaplotypeBits(k);
            for (int l=0; l<nBits; ++l) {
                if ((allele & mask)==mask) {
                    int wordIndex =  bitIndex >> LOG2_BITS_PER_WORD;
                    bits[wordIndex] |= (1L << bitIndex);
                }
                bitIndex++;
                mask <<= 1;
            }
        }
        return bits;
    }

    @Override
    public int allele1(int marker) {
        return allele(bits1, marker);
    }

    @Override
    public int allele2(int marker) {
        return allele(bits2, marker);
    }

    private int allele(long[] alleleBits, int marker) {
        int start = markers.sumHaplotypeBits(marker);
        int end = markers.sumHaplotypeBits(marker+1);
        if (end==(start+1)) {
            int wordIndex =  start >> LOG2_BITS_PER_WORD;
            return (int) (alleleBits[wordIndex] >> start) & 1;
        }
        int allele = 0;
        int mask = 1;
        for (int j=start; j<end; ++j) {
            int wordIndex =  j >> LOG2_BITS_PER_WORD;
            if ((alleleBits[wordIndex] & (1L << j)) != 0) {
                allele |= mask;
            }
            mask <<= 1;
        }
        return allele;
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public int idIndex() {
        return idIndex;
    }

    /**
     * Returns a string representation of {@code this}.  The
     * exact details of the representation are unspecified and subject
     * to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("BitHapPair: idIndex=");
        sb.append(idIndex);
        return sb.toString();
    }
}
