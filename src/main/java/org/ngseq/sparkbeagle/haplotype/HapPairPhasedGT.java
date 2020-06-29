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

import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.vcf.GT;
import org.ngseq.sparkbeagle.vcf.Marker;
import org.ngseq.sparkbeagle.vcf.Markers;

import java.util.Collections;
import java.util.List;

/**
 * <p>Class {@code HapPairPhasedGT} stores a list of samples and a
 * haplotype pair for each sample.
 * </p>
 * <p>Instance of class {@code HapPairPhasedGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HapPairPhasedGT implements GT {

    private final Markers markers;
    private final Samples samples;
    private final BitHapPair[] hapPairs;

    /**
     * Constructs a new {@code BasicPhasedGT} instance from the specified data.
     * @param samples a list of samples
     * @param hapPairList a list of haplotype pairs corresponding to the
     * specified list of samples
     *
     * @throws IllegalArgumentException if
     * {@code hapPairList.isEmpty() == true}
     * @throws IllegalArgumentException if
     * {@code hapPairList.get(j).markers().equals(hapPairList.get(k).markers())
     * == false}
     * for any indices {@code j, k} satisfying
     * {@code 0 <= j && j < k && k < hapPairList.size()}
     * @throws IllegalArgumentException if the list of samples does not
     * match the list of samples determined by {@code hapPairList}
     * @throws IllegalArgumentException if
     * {@code samples.nSamples() != hapPairs.size()}
     * @throws IllegalArgumentException if
     * {@code (samples.idIndex(j) != hapPairs.get(j).idIndex())} for any
     * {@code j} satisfying {@code 0 <= j && j < hapPairList.size()}
     * @throws NullPointerException if {@code samples == null}
     * @throws NullPointerException if
     * {@code (hapPairList == null || hapPairList(j) == null)}
     * for any {@code j} satisfying {@code (0 <= j && j < hapPairList.size())}
     */
    public HapPairPhasedGT(Samples samples, List<BitHapPair> hapPairList) {
        if (hapPairList.isEmpty()) {
            throw new IllegalArgumentException("haps.isEmpy()==true");
        }
        Collections.sort(hapPairList, HapPair.comparator(samples));
        checkSamples(samples, hapPairList);
        this.markers = checkAndExtractMarkers(hapPairList);
        this.samples = samples;
        this.hapPairs = hapPairList.toArray(new BitHapPair[0]);
    }

    private static void checkSamples(Samples samples, List<BitHapPair> hapPairs) {
        if (samples.nSamples()!=hapPairs.size()) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        for (int j=0, n=hapPairs.size(); j<n; ++j) {
            if (samples.idIndex(j)!=hapPairs.get(j).idIndex()) {
                throw new IllegalArgumentException("inconsistent samples");
            }
        }
    }

    /**
     * Checks that all haplotype pairs have alleles for the same list of
     * markers, and returns the list of markers.
     * @param hapPairList a list of haplotype pairs
     * @return the list of markers shared by the specified haplotype pairs
     * @throws IllegalArgumentException if
     * {@code hapPiarList.get(j).markers().equals(hapPairList.get(k).markers())
     * == false}
     * for any indices {@code j, k} satisfying
     * {@code 0 <= j && j < k && k < hapPairList.size()}
     * @throws NullPointerException if
     * {@code hapPairList == null || hapPairList(j) == null}
     * for any {@code j} satisfying {@code 0 <= j && j < hapPairList.size()}
     */
    static Markers checkAndExtractMarkers(List<BitHapPair> hapPairList) {
        if (hapPairList.isEmpty()) {
            return Markers.create(new Marker[0]);
        }
        else {
            Markers m = hapPairList.get(0).markers();
            for (int j=1, n=hapPairList.size(); j<n; ++j) {
                if (hapPairList.get(j).markers().equals(m)==false) {
                    throw new IllegalArgumentException("inconsistent markers");
                }
            }
            return m;
        }
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
        return hapPairs[hapPair].allele1(marker);
    }

    @Override
    public int allele2(int marker, int hapPair) {
        return hapPairs[hapPair].allele2(marker);
    }

    @Override
    public int allele(int marker, int haplotype) {
        int pairIndex = haplotype/2;
        if ((haplotype & 1)==0) {
            return hapPairs[pairIndex].allele1(marker);
        }
        else {
            return hapPairs[pairIndex].allele2(marker);
        }
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
        return 2*hapPairs.length;
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
