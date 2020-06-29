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

import static org.ngseq.sparkbeagle.vcf.LowMafRefGT.checkIndicesAndReturnMajorAllele;

/**
 * <p>Class {@code LowMafRefDiallelicGT} represent represents phased,
 * non-missing genotypes for a list of reference samples at a single diallelic
 * marker.  Genotype emission probabilities are determined by the sample
 * genotypes.
 * </p>
 * <p>
 * Class {@code LowMafRefDiallelicGT} stores the minor allele indices.
 * </p>
 * <p>Instances of class {@code LowMemRefDiallelicGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LowMafRefDiallelicGT implements RefGTRec {

    private final Marker marker;
    private final Samples samples;
    private final int nHaps;
    private final int majorAllele;
    private final int minorAllele;
    private final int[] minorAlleles;

    /**
     * Constructs a new {@code LowMafRefDiallelicGT} instance with phased
     * non-missing genotypes from the specified data.
     *
     * @param rec the phased, non-missing genotype data
     * @throws IllegalArgumentException if {@code rec.marker().nAlleles() != 2}
     * @throws NullPointerException if {@code rec == null}
     */
    public LowMafRefDiallelicGT(RefGTRec rec) {
        if (rec.marker().nAlleles()!=2) {
            throw new IllegalArgumentException(
                    String.valueOf(rec.marker().nAlleles()!=2));
        }
        int[][] hapIndices = rec.hapIndices();
        int majAllele = 0;
        while (hapIndices[majAllele]!=null) {
            ++majAllele;
        }
        this.marker = rec.marker();
        this.samples = rec.samples();
        this.nHaps = rec.size();
        this.majorAllele = majAllele;
        this.minorAllele = 1 - majAllele;
        this.minorAlleles = hapIndices[minorAllele];
    }

    /**
     * Constructs a new {@code LowMafRefDiallelicGT} instance with phased
     * non-missing genotypes from the specified data.
     *
     * @param gtp a VCF record parser that extracts sample genotypes
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased genotype or missing allele
     * @throws IllegalArgumentException if {@code gtp.nAlleles() != 2}
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws NullPointerException if {@code gtp == null}
     */
    public LowMafRefDiallelicGT(VcfRecGTParser gtp) {
        if (gtp.nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(gtp.nAlleles()));
        }
        int[][] nonMajIndices = gtp.nonMajRefIndices();
        this.marker = gtp.marker();
        this.samples = gtp.samples();
        this.nHaps = 2*gtp.nSamples();
        this.majorAllele = nonMajIndices[0]==null ? 0 : 1;
        this.minorAllele = 1 - majorAllele;
        this.minorAlleles = nonMajIndices[minorAllele];
    }

    /**
     * Constructs a new {@code LowMafRefDiallelicGT} instance from the
     * specified data. The contract for this class is unspecified if a
     * haplotype index is duplicated in the specified {@code hapIndices} array.
     *
     * @param marker the marker
     * @param samples the samples
     * @param hapIndices an array whose {@code j}-th element is {@code null}
     * if {@code j} is the major allele with lowest index, and otherwise is
     * an array of indices of haplotypes that carry the {@code j}-th allele
     * sorted in increasing order
     *
     * @throws IllegalArgumentException if {@code marker.nAlleles() != 2}
     * @throws IllegalArgumentException if the {@code (hapIndices[j] == null)}
     * and {@code j} is not the major allele with lowest index, or if
     * {@code (hapIndices[j] != null)} and {@code j} is the major allele with
     * lowest index
     * @throws IllegalArgumentException if any element of {@code hapIndices}
     * is not not a sorted array of distinct haplotype indices
     * between 0 (inclusive) and {@code 2*samples.nSamples()} (exclusive)
     * @throws IllegalArgumentException if
     * {@code marker.nAlleles() != hapIndices.length}
     * @throws NullPointerException if
     * {@code marker == null || samples == null || hapIndices == null}
     */
    public LowMafRefDiallelicGT(Marker marker, Samples samples,
            int[][] hapIndices) {
        if (marker.nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(marker.nAlleles()));
        }
        this.marker = marker;
        this.samples = samples;
        this.nHaps = 2*samples.nSamples();
        this.majorAllele = checkIndicesAndReturnMajorAllele(hapIndices, nHaps);
        this.minorAllele = 1 - majorAllele;
        this.minorAlleles = hapIndices[minorAllele].clone();
    }

    @Override
    public int[][] hapIndices() {
        int[][] hapIndices = new int[2][];
        hapIndices[minorAllele] = minorAlleles.clone();
        return hapIndices;
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
    public Samples samples() {
        return samples;
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public int size() {
        return nHaps;
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
        if (allele1 != 0 &&  allele1 != 1) {
            throw new IndexOutOfBoundsException(String.valueOf(allele1));
        }
        if (allele2 != 0 &&  allele2 != 1) {
            throw new IndexOutOfBoundsException(String.valueOf(allele2));
        }
        boolean matches = (allele1==allele1(sample) && allele2==allele2(sample));
        return matches ? 1.0f : 0.0f;
    }

    @Override
    public int allele1(int sample) {
        return get(2*sample);
    }

    @Override
    public int allele2(int sample) {
        return get(2*sample + 1);
    }

    @Override
    public int get(int hap) {
        if (hap < 0 || hap >= nHaps) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        if (Arrays.binarySearch(minorAlleles, hap) >= 0) {
            return minorAllele;
        }
        else {
            return majorAllele;
        }
    }

    @Override
    public int[] alleles() {
        int[] ia = IntStream.range(0, nHaps)
                .map(h -> majorAllele)
                .toArray();
        for (int h : minorAlleles) {
            ia[h] = minorAllele;
        }
        return ia;
    }

    @Override
    public int nAlleles() {
        return this.marker().nAlleles();
    }

    @Override
    public boolean isAlleleCoded() {
        return true;
    }

    @Override
    public int majorAllele() {
        return majorAllele;
    }

    @Override
    public int alleleCount(int allele) {
        if (allele==majorAllele) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return minorAlleles.length;
        }
    }

    @Override
    public int hapIndex(int allele, int copy) {
        if (allele==majorAllele) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return minorAlleles[copy];
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
        return 1;
    }

    @Override
    public IntArray[] maps() {
        return new IntArray[] {toIntArray()};
    }

    @Override
    public IntArray map(int index) {
        if (index!=0) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return toIntArray();
    }

    private IntArray toIntArray() {
        int[] ia = IntStream.range(0, nHaps)
                .map(i -> majorAllele)
                .toArray();
        for (int i : minorAlleles) {
            ia[i] = minorAllele;
        }
        return IntArray.create(ia, 2);
    }
}
