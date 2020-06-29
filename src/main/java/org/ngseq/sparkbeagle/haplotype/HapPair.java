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
import org.ngseq.sparkbeagle.vcf.Marker;
import org.ngseq.sparkbeagle.vcf.Markers;

import java.util.Comparator;

/**
 * <p>Interface {@code HapPair} represents a pair of haplotypes for a sample.
 * The pair of haplotypes are guaranteed to have non-missing alleles at each
 * marker.
 * </p>
 * All instances of {@code HapPair} are required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface HapPair {

    /**
     * Returns the first allele for the specified marker.
     * @param marker a marker index
     * @return the first allele for the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    int allele1(int marker);

    /**
     * Returns the second allele for the specified marker.
     * @param marker a marker index
     * @return the second allele for the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    int allele2(int marker);

    /**
     * Returns the markers.
     * @return the markers
     */
    Markers markers();

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
     Marker marker(int marker);

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the sample identifier index.
     * @return the sample identifier index
     */
    int idIndex();

    /**
     * Returns a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method returns -1, 0, or 1
     * depending on whether {@code samples.index(hp1.idIndex())} is
     * less than, equal, or greater than
     * {@code samples.index(hp2.idIndex())}.
     * @param samples the list of samples used to compare {@code HapsPair}
     * objects
     * @return a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method compares two
     * haplotype pairs for order
     * @throws NullPointerException if {@code samples == null}
     */
    static Comparator<HapPair> comparator(final Samples samples) {
        return (hp1, hp2) -> Integer.compare(
                samples.index(hp1.idIndex()),
                samples.index(hp2.idIndex()));
    }
}
