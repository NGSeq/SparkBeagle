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

/**
 * <p>Class {@code SplicedGT} represents genotype emission probabilities
 * for a set of samples. The genotype emission probabilities are determined
 * by a {@code RefGTWindow} instance for the initial markers, and are
 * determined by a {@code GT} instance for the remaining markers.
 * The {@code isPhased()} and {@code isPhased(sample)} methods of the
 * {@code SplicedGT} class return the same values as the corresponding
 * methods of the {@code GT} instance.
 * </p>
 * <p>Instances of class {@code SplicedGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class SplicedGT implements GT {

    private final int overlap;
    private final GT phasedGT;
    private final GT gt;

    /**
     * Constructs a new {@code SplicedGL} instance.
     * @param phasedOverlap sample haplotype pairs for the initial
     * markers
     * @param gt genotype emission probabilities for all markers
     * @throws IllegalArgumentException if
     * {@code phasedOverlaps.nMarkers() >= gt.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code phasedOverlaps.marker(j).equals(gt.marker(j)) == false} for any {@code j}
     * satisfying {@code 0 <= j && j < phasedOverlaps.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code phasedOverlaps.samples().equals(gt.samples()) == false}
     * @throws IllegalArgumentException if
     * {@code phasedOverlap.isPhased() == false}
     * @throws NullPointerException if {@code phasedOverlap == null || gt == null}
     */
    public SplicedGT(GT phasedOverlap, GT gt) {
        if (phasedOverlap.nMarkers()>=gt.nMarkers()) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        for (int j=0, n=phasedOverlap.nMarkers(); j<n; ++j) {
            if (phasedOverlap.marker(j).equals(gt.marker(j))==false) {
                throw new IllegalArgumentException("inconsistent markers");
            }
        }
        if (phasedOverlap.samples().equals(gt.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (phasedOverlap.isPhased()==false) {
            throw new IllegalArgumentException("unphased data");
        }
        this.overlap = phasedOverlap.nMarkers();
        this.phasedGT = phasedOverlap;
        this.gt = gt;
    }

    @Override
    public boolean isPhased() {
        return gt.isPhased();
    }

    @Override
    public boolean isGTData() {
        return gt.isGTData();
    }

    @Override
    public boolean isPhased(int sample) {
        return gt.isPhased(sample);
    }

    @Override
    public float gl(int marker, int sample, int allele1, int allele2) {
        if (marker<overlap) {
            int a1 = phasedGT.allele1(marker, sample);
            int a2 = phasedGT.allele2(marker, sample);
            return (allele1==a1 && allele2==a2) ? 1.0f : 0.0f;
        }
        else {
            return gt.gl(marker, sample, allele1, allele2);
        }
    }

    @Override
    public boolean isPhased(int marker, int sample) {
        if (marker<overlap) {
            return true;
        }
        else {
            return gt.isPhased(marker, sample);
        }
    }

    @Override
    public int allele1(int marker, int sample) {
        if (marker<overlap) {
            return phasedGT.allele1(marker, sample);
        }
        else {
            return gt.allele1(marker, sample);
        }
    }

    @Override
    public int allele2(int marker, int sample) {
        if (marker<overlap) {
            return phasedGT.allele2(marker, sample);
        }
        else {
            return gt.allele2(marker, sample);
        }
    }

    @Override
    public int allele(int marker, int hap) {
        if (marker<overlap) {
            return phasedGT.allele(marker, hap);
        }
        else {
            return gt.allele(marker, hap);
        }
    }

    @Override
    public Marker marker(int marker) {
        return gt.marker(marker);
    }

    @Override
    public Markers markers() {
        return gt.markers();
    }

    @Override
    public int nMarkers() {
        return gt.nMarkers();
    }

    @Override
    public int nHaps() {
        return gt.nHaps();
    }

    @Override
    public int nSamples() {
        return gt.nSamples();
    }

    @Override
    public Samples samples() {
        return gt.samples();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(10000);
        sb.append("SplicedGL: nSamples=");
        sb.append(this.nSamples());
        return sb.toString();
    }
}
