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

/**
 * Interface {@code MarkerContainer} represents an object that stores
 * a unique {@code vcf.Marker} instance.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface MarkerContainer {

    /**
     * Returns the marker.
     * @return the marker
     */
    Marker marker();

    /**
     * Returns the number of marker alleles.
     * @return the number of marker alleles.
     */
    int nAlleles();
}
