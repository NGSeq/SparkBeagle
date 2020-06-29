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

import org.ngseq.sparkbeagle.Pedigree;
import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.haplotype.HapPair;

import java.io.Closeable;
import java.util.List;

/**
 * Interface {@code Data} represents a sliding window of target VCF records
 * or a sliding window of reference and target VCF records.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface Data extends Closeable {

    /**
     * Returns the pedigree.
     * @return the pedigree
     */
    Pedigree ped();

    /**
     * Returns the genetic map.
     * @return the genetic map
     */
    GeneticMap genMap();

    /**
     * Returns {@code true} if the current window of VCF records is the last
     * window for the chromosome and returns {@code false} otherwise.
     * @return {@code true} if the current window of VCF records is the last
     * window for the chromosome
     */
    boolean lastWindowOnChrom();

    /**
     * Returns {@code true} if the sliding window of VCF records can advance
     * and returns {@code false} otherwise.
     * @return {@code true} if the sliding window of VCF records can advance
     */
    boolean canAdvanceWindow();

    /**
     * Advances the sliding window of VCF records.
     *
     * @throws IllegalArgumentException if a format error in the input data
     * is detected
     * @throws IllegalStateException if
     * {@code this.canAdvanceWindow() == false}
     */
    void advanceWindowCm();

    /**
     * Returns the current window index.  The first window has index 1.
     * @return the current window index
     */
    public int windowIndex();

     /**
     * Returns the number of target data markers in the overlap between
     * the current marker window and the previous marker window.
     * Returns 0 if the current marker window is the first marker window.
     *
     * @return the number of target data markers in the overlap between
     * the current marker window and the previous marker window
     */
    int targetOverlap();

    /**
     * Returns the number of VCF records in the overlap between the current
     * window and the previous window.  Returns 0 if the current window
     * is the first window.
     *
     * @return the number of VCF records in the overlap between the current
     * window and the previous window
     */
    public int overlap();

    /**
     * Returns the first marker index in the overlap between this
     * marker window and the next marker window. Returns
     * {@code this.nMarkers()} if the next marker window is from a
     * different chromosome.
     * @return the first marker index in the overlap between this
     * marker window and the next marker window
     */
    public int nextOverlapStart();

     /**
     * Returns the number of target data markers in the current window.
     * @return the number of target data markers in the current window
     */
    int nTargetMarkers();

    /**
     * Returns the number of target VCF records in the union of the
     * current window and all previous windows.
     * @return the number of target VCF records in the union of the
     * current window and all previous windows
     */
    int nTargetMarkersSoFar();

    /**
     * Returns the list of target data markers in the current window.
     * @return the list of target data markers in the current window
     */
    Markers targetMarkers();

    /**
     * Returns the number of markers in the current window.
     * @return the number of markers in the current window
     */
     int nMarkers();

    /**
     * Returns the number of markers in the union of the current window
     * and all previous windows.
     * @return the number of markers in the union of the current window
     * and all previous windows
     */
    int nMarkersSoFar();

    /**
     * Returns the list of markers in the current window.
     * @return the list of markers in the current window
     */
     Markers markers();

     /**
      * Returns the target data marker index corresponding to the specified
      * marker, or returns -1 if no corresponding  target data marker exists.
      * Indices are with respect to the current window.
      * @param marker a marker index
      * @return the target data marker index corresponding to the specified
      * marker, or returns -1 if no corresponding  target data marker exists
      * @throws IndexOutOfBoundsException if
      * {@code marker < 0 || marker >= this.nMarkers()}
      */
     int targetMarkerIndex(int marker);

     /**
      * Returns the marker index corresponding to the
      * specified target data marker.  Indices are with
      * respect to the current window.
      * @param targetMarker a target data marker index
      * @return the marker index corresponding to the specified
      * target data marker
      * @throws IndexOutOfBoundsException if
      * {@code targetMarker < 0 || targetMarker >= this.nTargetMarkers()}
      */
     int markerIndex(int targetMarker);

    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */
    int nTargetSamples();

    /**
     * Returns the list of target samples.
     * @return the list of target samples
     */
    Samples targetSamples();

    /**
     * Returns the number of reference samples.
     * @return the number of reference samples
     */
    int nRefSamples();

    /**
     * Returns the list of reference samples, or {@code null} if
     * there are no reference samples.
     * @return the list of reference samples, or {@code null} if
     * there are no reference samples
     */
    Samples refSamples();

     /**
      * Returns the total number of reference and target samples.
      * @return the total number of reference and target samples
      */
     int nAllSamples();

     /**
      * Returns a list of all target and reference samples.
      * Target samples are listed first in the same order as the list returned
      * by {@code this.targetSamples()}. Reference samples are listed last
      * in the same order as the list returned by {@code this.refSamples()}.
      * @return a list of all target and reference samples
      */
     Samples allSamples();

    /**
     * Returns the genotype likelihoods for the target samples
     * restricted to the target data markers in the current window.
     * The returned {@code GL} instance will contain no markers if
     * {@code this.advanceWindow()} has not yet been invoked.
     * @return the genotype likelihoods for the target samples
     * restricted to the target data markers in the current window
     */
    GT targGT();

    /**
     * Returns a list of reference haplotype pairs that are restricted
     * to the target data markers in the current window.
     * The returned list will be empty if there are no reference samples
     * or if {@code this.advanceWindow()} has not yet been invoked.
     * @return a list of reference haplotype pairs that are restricted
     * to the target data markers
     */
    List<HapPair> restrictedRefHapPairs();

    /**
     * Returns a list of the reference haplotype pairs for the current
     * window.  The returned list will be empty if there are no reference
     * samples or if {@code this.advanceWindow()} has not yet been invoked.
     * @return a list of the reference haplotype pairs
     */
    List<HapPair> refHapPairs();

    /**
     * Returns the  phased, nonmissing reference genotype data
     * for the current window, or {@code null} if there are no reference data
     * @return the reference genotype data for the current window or
     * {@code null} if there are no reference data
     */
    RefGT refGT();

    /**
     * Returns the  phased, nonmissing reference genotype data
     * for the target data markers in the current window.
     * Returns {@code null} if there are no reference data
     * @return the reference genotype data for the target data markers
     * in the current window or {@code null} if there are no reference data
     */
    RefGT restrictRefGT();

    /**
     * Releases any I/O resources controlled by this object.
     */
    @Override
    void close();
}
