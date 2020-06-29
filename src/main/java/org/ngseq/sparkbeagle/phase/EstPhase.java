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
package org.ngseq.sparkbeagle.phase;

import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.ints.IntArray;
import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.ints.LongArray;
import org.ngseq.sparkbeagle.vcf.GT;
import org.ngseq.sparkbeagle.vcf.Marker;
import org.ngseq.sparkbeagle.vcf.Markers;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicReferenceArray;

/**
 * <p>Class {@code EstPhase} stores haplotype pairs, missing genotypes,
 * and unphased, nonmissing heteroygote genotypes for a list of samples.
 * </p>
 * <p>Instances of class {@code EstPhase} are thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class EstPhase {

    private final Samples samples;
    private final Markers markers;
    private final AtomicReferenceArray<LongArray> phase;
    private final AtomicReferenceArray<IntArray> unphased;
    private final AtomicReferenceArray<IntArray> missing;


    /**
     * Constructs a new {@code EstPhase} instance from the specified data.
     * The haplotypes for the {@code k}-th sample must be stored in
     * entries {@code haps.get(2*k}} and {@code haps.get(2*k + 1)}
     * @param gt the observed genotype data
     * @param haps the initial haplotypes for each sample
     * @throws IllegalArgumentException if
     * {@code gt.nSamples() != hapPairs.nSamples()}
     * @throws IllegalArgumentException if
     * {@code (hapPairs.get(s).idIndex!= this.gt().samples().idIndex(s))}
     * for any {@code s} satisfying {@code (0 <= s && s < gt.nSamples())}
     * @throws NullPointerException if {@code gt == null || hapPairs == null}
     * @throws NullPointerException if {@code (hapPairs.get(j) == null)} for any
     * {@code j} satisfying {@code (0 <= j && j < hapPairs.size())}
     */
    public EstPhase(GT gt, List<LongArray> haps) {
        if (gt.nHaps()!=haps.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        int nSamples = gt.nSamples();
        this.samples = gt.samples();
        this.markers = gt.markers();
        this.phase = new AtomicReferenceArray<>(gt.nHaps());
        this.unphased = new AtomicReferenceArray<>(nSamples);
        this.missing = new AtomicReferenceArray<>(nSamples);

        IntList unphList = new IntList(gt.nMarkers()/10);
        IntList missList = new IntList(gt.nMarkers()/100);
        for (int s=0; s<nSamples; ++s) {
            unphList.clear();
            missList.clear();
            boolean firstHetFound = false;
            if (gt.isPhased(s)==false) {
                for (int m=0, n=gt.nMarkers(); m<n; ++m) {
                    int a1 = gt.allele1(m, s);
                    int a2 = gt.allele2(m, s);
                    if (a1<0 || a2<0) {
                        missList.add(m);
                    }
                    else if (a1!=a2) {
                        if (firstHetFound) {
                            unphList.add(m);
                        }
                        else {
                            firstHetFound = true;
                        }
                    }
                }
            }
            int h1 = 2*s;
            int h2 = 2*s + 1;
            phase.set(h1, haps.get(h1));
            phase.set(h2, haps.get(h2));
            missing.set(s, IntArray.create(missList, gt.nMarkers()));
            unphased.set(s, IntArray.create(unphList, gt.nMarkers()));
        }
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return samples.nSamples();
    }
    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */

    /**
     * Returns the list of target samples.
     * @return the list of target samples
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return markers.nMarkers();
    }

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    public Markers markers() {
        return markers;
    }

    /**
     * Returns {@code true} if the specified sample has an unphased
     * heterozygote genotype, and returns {@code false} otherwise.
     * @param sample the sample index
     * @return {@code true} if the specified sample has an unphased
     * heterozygote genotype
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    public boolean hasUnphasedHets(int sample) {
        return unphased.get(sample).size()>0;
    }

    /**
     * Returns the number of unphased, non-missing heterozygote
     * genotypes.
     * @param sample the sample index
     * @return the number of unphased, non-missing heterozygote
     * genotypes
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    public int nUnphasedHets(int sample) {
        return unphased.get(sample).size();
    }

    /**
     * Returns a list of marker indices in increasing order for which
     * the specified sample has an unphased, non-missing heterozygote genotype.
     * @param sample the sample index
     * @return a list of marker indices in increasing order for which
     * the specified sample has an unphased heterozygote genotype
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    public IntArray getUnphasedHets(int sample) {
        return unphased.get(sample);
    }

    /**
     * Sets the array with marker indices of unphased, nonmissing heterozygote
     * genotypes to the specified array.
     * @param sample the sample index
     * @param newUnphased the marker indices of unphased heterozygote genotypes
     * @throws IllegalArgumentException if the specified {@code newUnphased}
     * list is not sorted in increasing order, contains a duplicate elements,
     * or is not a subset of {@code this.getUnphased(sample)}.
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws NullPointerException if {@code newUnphased == null}
     */
    public void setUnphasedHets(int sample, IntArray newUnphased) {
        IntArray oldUnphased = unphased.get(sample);
        int oldSize = oldUnphased.size();
        int nextOldIndex = 0;
        int oldMkr = -1;
        for (int j=0, n=newUnphased.size(); j<n; ++j) {
            int newMkr = newUnphased.get(j);
            if (nextOldIndex<oldSize) {
                oldMkr = oldUnphased.get(nextOldIndex++);
            }
            while (oldMkr!=newMkr && nextOldIndex<oldSize) {
                oldMkr = oldUnphased.get(nextOldIndex++);
            }
            if (oldMkr!=newMkr) {
                throw new IllegalArgumentException(newUnphased.toString());
            }
        }
        unphased.set(sample, newUnphased);
    }

    /**
     * Sets the haplotype pair for the specified sample to the specified
     * haplotypes.
     * @param  sample the sample index
     * @param hap1 an array whose {@code k}-th entry is the estimated allele
     * carried by the sample's first haplotype.
     * @param hap2 an array whose {@code k}-th entry is the estimated allele
     * carried by the sample's second haplotype.
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IllegalArgumentException if
     * {@code hap1.length != this.nMarkers() || hap2.length != this.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code (hap1[k] < 0 || hap1[k] >= this.markers().marker(k).nAlleles())}
     * for any {@code k} satisfying {@code (0 <= k && k < this.nMarkers())}
     * @throws IllegalArgumentException if
     * {@code (hap2[k] < 0 || hap2[k] >= this.markers().marker(k).nAlleles())}
     * for any {@code k} satisfying {@code (0 <= k && k < this.nMarkers())}
     * @throws NullPointerException if {@code hap1 == null || hap2 == null}
     */
    public void setHapPair(int sample, int[] hap1, int[] hap2) {
        int h1 =  2*sample;
        phase.set(h1, markers.allelesToBits(hap1));
        phase.set(h1 + 1, markers.allelesToBits(hap2));
    }

    /**
     * Sets the {@code k}-th entry of the specified {@code isUnphased} array
     * to {@code true} if the specified sample's genotype at the
     * {@code k}-th marker is an unphased, nonmissing heterozygote, and
     * to {@code false} otherwise.
     * @param sample the sample index
     * @param isUnphased an array whose {@code k}-th entry
     * will be set to {@code true} if the specified sample's genotype at the
     * {@code k}-th marker is an unphased, nonmissing heterozygote,
     * and to {@code false} otherwise.
     * @throws IllegalArgumentException if
     * {@code isUnphased.length != this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws NullPointerException if {@code isUnphased == null}
     */
    public void getUnphasedHets(int sample, boolean[] isUnphased) {
        if (isUnphased.length != markers.nMarkers()) {
            throw new IllegalArgumentException(String.valueOf(isUnphased.length));
        }
        Arrays.fill(isUnphased, false);
        IntArray unph = unphased.get(sample);
        for (int j=0, n=unph.size(); j<n; ++j) {
            isUnphased[unph.get(j)] = true;
        }
    }

    /**
     * Returns {@code true} if the specified sample has a missing genotype,
     * and returns {@code false} otherwise.
     * @param sample the sample index
     * @return {@code true} if the specified sample has a missing genotype
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    public boolean hasMissingGT(int sample) {
        return missing.get(sample).size()>0;
    }

    /**
     * Returns a list of marker indices in increasing order for which
     * the specified sample has a missing genotype.
     * @param sample the sample index
     * @return a list of marker indices in increasing order for which
     * the specified sample has a missing genotype
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    public IntArray getMissing(int sample) {
        return missing.get(sample);
    }

    /**
     * Sets the {@code k}-th entry of the specified {@code isMissing} array
     * to {@code true} if the specified sample's genotype at the
     * {@code k}-th marker is missing, and to {@code false} otherwise.
     * @param sample the sample index
     * @param isMissing an array whose {@code k}-th entry will be
     * set to {@code true} if the specified sample's genotype at the
     * {@code k}-th marker is missing and to {@code false} otherwise
     * @throws IllegalArgumentException if
     * {@code isMissing.length != this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws NullPointerException if {@code isMissing == null}
     */
    public void getMissing(int sample, boolean[] isMissing) {
        if (isMissing.length != markers.nMarkers()) {
            throw new IllegalArgumentException(String.valueOf(isMissing.length));
        }
        Arrays.fill(isMissing, false);
        IntArray miss = missing.get(sample);
        for (int j=0, n=miss.size(); j<n; ++j) {
            isMissing[miss.get(j)] = true;
        }
    }

    /**
     * Sets the {@code k}-th element of the specified {@code hap1} and
     * {@code hap2} arrays to the specified sample's phased genotype
     * at the {@code k}-th marker if the genotype is nonmissing, and sets the
     * {@code k}-th element of {@code hap1} and {@code hap2} to
     * {@code -1} otherwise.
     * @param sample the sample index
     * @param hap1 an array whose {@code k}-th entry will be set to the
     * current estimated allele carried by the sample's first haplotype if the
     * observed genotype is non-missing and to -1 otherwise.
     * @param hap2 an array whose {@code k}-th entry will be set to the
     * current estimated allele carried by the sample's second haplotype if the
     * observed genotype is non-missing and to -1 otherwise.
     * @throws IllegalArgumentException if
     * {@code hap1.length != this.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code hap2.length != this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws NullPointerException if {@code hap1 == null || hap2 == null}
     */
    public void setHapsWithMaskedMissingAlleles(int sample, int[] hap1, int[] hap2) {
        if (hap1.length != markers.nMarkers()) {
            throw new IllegalArgumentException(String.valueOf(hap1.length));
        }
        if (hap2.length != markers.nMarkers()) {
            throw new IllegalArgumentException(String.valueOf(hap2.length));
        }
        LongArray bits1 = phase.get(2*sample);
        LongArray bits2 = phase.get(2*sample + 1);
        for (int m=0, nMarkers=markers.nMarkers(); m<nMarkers; ++m) {
            hap1[m] = markers.bitsToAllele(bits1, m);
            hap2[m] = markers.bitsToAllele(bits2, m);
        }
        IntArray ia = missing.get(sample);
        for (int j=0; j<ia.size(); ++j) {
            int m = ia.get(j);
            hap1[m] = hap2[m] = -1;
        }
    }

    public HapsGT hapsGT() {
        LongArray[] haps = new LongArray[phase.length()];
        for (int j=0; j<haps.length; ++j) {
            haps[j] = phase.get(j);
        }
        return new HapsGT(markers, samples, haps);
    }

    public final class HapsGT implements GT {

        private final Markers markers;
        private final Samples samples;
        private final LongArray[] haps;

        private HapsGT(Markers markers, Samples samples, LongArray[] haps) {
            this.markers = markers;
            this.samples = samples;
            this.haps = haps;
        }

        @Override
        public boolean isPhased(int marker, int sample) {
            if (marker < 0 || marker >= this.nMarkers()) {
                throw new IndexOutOfBoundsException(String.valueOf(marker));
            }
            if (sample < 0 || sample >= this.nSamples()) {
                throw new IndexOutOfBoundsException(String.valueOf(sample));
            }
            return true;
        }

        @Override
        public boolean isPhased(int sample) {
            if (sample < 0 || sample >= this.nSamples()) {
                throw new IndexOutOfBoundsException(String.valueOf(sample));
            }
            return true;
        }

        @Override
        public boolean isPhased() {
            return true;
        }

        @Override
        public int allele1(int marker, int hapPair) {
            return markers.bitsToAllele(haps[hapPair<<1], marker);
        }

        @Override
        public int allele2(int marker, int hapPair) {
            return markers.bitsToAllele(haps[(hapPair<<1) + 1], marker);
        }

        @Override
        public int allele(int marker, int haplotype) {
            return markers.bitsToAllele(haps[haplotype], marker);
        }

        @Override
        public int nMarkers() {
            return markers.nMarkers();
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
        public int nHaps() {
            return haps.length;
        }

        @Override
        public int nSamples() {
            return samples.nSamples();
        }

        @Override
        public Samples samples() {
            return samples;
        }
    }
}
