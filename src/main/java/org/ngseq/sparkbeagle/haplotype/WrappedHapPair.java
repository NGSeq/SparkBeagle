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

import org.ngseq.sparkbeagle.vcf.GT;
import org.ngseq.sparkbeagle.vcf.Marker;
import org.ngseq.sparkbeagle.vcf.Markers;

/**
 * Class {@code WrappedHapPair} is a {@code HapPair} instance
 * that wraps a {@code RefGTWindow} object.

* @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class WrappedHapPair implements HapPair {

    private final GT phasedGT;
    private final int hapPair;

    /**
     * Creates a {@code WrappedHapPair} instance representing
     * the specified haplotype pair.
     * @param phasedGT the {@code RefGTWindow} object that
     * will be "wrapped" by {@code this}
     * @param hapPair a haplotype pair index
     * @throws IllegalArgumentException if {@code phasedGT.isPhased() == false}
     * @throws IllegalArgumentException if
     * {@code hapPair < 0 || hapPair >= sampleHapPairs.nHapPairs()}
     * @throws NullPointerException if {@code sampleHapPairs == null}
     */
    public WrappedHapPair(GT phasedGT, int hapPair) {
        if (phasedGT.isPhased()==false) {
            throw new IllegalArgumentException("unphased data");
        }
        if (hapPair < 0 || hapPair >= phasedGT.nSamples()) {
            throw new IllegalArgumentException("hapPair: " + hapPair);
        }
        this.phasedGT = phasedGT;
        this.hapPair = hapPair;
    }

    @Override
    public int allele1(int marker) {
        return phasedGT.allele1(marker, hapPair);
    }

    @Override
    public int allele2(int marker) {
        return phasedGT.allele2(marker, hapPair);
    }

    @Override
    public Markers markers() {
        return phasedGT.markers();
    }

    @Override
    public Marker marker(int marker) {
        return phasedGT.marker(marker);
    }

    @Override
    public int nMarkers() {
        return phasedGT.nMarkers();
    }

    @Override
    public int idIndex() {
        return phasedGT.samples().idIndex(hapPair);
    }
}
