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
package org.ngseq.sparkbeagle.imp;

/**
 * <p>Class {@code StateProbs} stores a subset of Li and Stephens HMM states
 * and associated probabilities for a target haplotype.
 * </p>
 * <p>All instances of interface {@code StateProbs} are required to be
 * immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface StateProbs {

    /**
     * Returns the target haplotype index.
     * @return the target haplotype index
     */
    int targHap();

    /**
     * Returns the number of target markers.
     * @return the number of target markers
     */
    int nTargMarkers();

    /**
     * Returns the number of stored HMM states at the specified
     * target marker.
     * @param targMarker a target marker index
     * @return the number of stored HMM states at the specified target marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nTargMarkers()}
     */
    int nStates(int targMarker);

    /**
     * Returns the specified reference haplotype index.
     * @param targMarker a target marker index
     * @param index a stored state index at the specified target marker
     * @return the specified reference haplotype index
     * @throws IndexOutOfBoundsException if
     * {@code targMarker < 0 || targMarker >= this.nTargMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nStates(targMarker)}
     */
    int refHap(int targMarker, int index);

    /**
     * Returns the probability of the specified state at the specified target
     * marker.
     * @param targMarker a target marker index
     * @param index a stored state index
     * @return the probability of the specified state at the specified target
     * marker
     * @throws IndexOutOfBoundsException if
     * {@code targMarker < 0 || targMarker >= this.nTargMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nStates(targMarker)}
     */
    float probs(int targMarker, int index);

    /**
     * Returns the probability of the specified state at the marker following
     * the specified target marker.  If
     * {@code (targMarker + 1 == this.nTargMarkers())}, the probability of
     * the specified state at the specified target marker is returned.
     * @param targMarker a target marker index
     * @param index a stored state index
     * @return the probability of the specified state at the marker following
     * the specified target marker
     * @throws IndexOutOfBoundsException if
     * {@code targMarker < 0 || targMarker >= this.nTargMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nStates(targMarker)}
     */
    float probsP1(int targMarker, int index);
}
