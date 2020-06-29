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
package org.ngseq.sparkbeagle;

import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.blbutil.FloatArray;
import org.ngseq.sparkbeagle.haplotype.BitHapPair;
import org.ngseq.sparkbeagle.haplotype.HapPair;
import org.ngseq.sparkbeagle.vcf.*;

import java.util.ArrayList;
import java.util.List;

/**
 * <p>Class {@code CurrentData} represents input data for the current marker
 * window.  All marker indices returned my methods of class {@code CurrentData}
 * are indexed with respect to the current marker window.
 * </p>
 * <p>Instances of class {@code CurrentData} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class CurrentData {

    private final Par par;
    private final int window;
    private final int prevSpliceStart;
    private final int nextOverlapStart;
    private final int nextSpliceStart;
    private final int prevTargSpliceStart;
    private final int nextTargetSpliceStart;
    private final int nextTargetOverlapStart;

    private final GT targGT;
    private final Pedigree ped;

    private final Samples refSamples;
    private final Samples targetSamples;
    private final Samples allSamples;
    private final int nHaps;

    private final Markers markers;
    private final Markers targetMarkers;
    private final int[] targMarkerIndex;
    private final int[] markerIndex;

    private final List<HapPair> restRefHapPairs;
    private final RefGT refGT;
    private final RefGT restrictRefGT;

    private final MarkerMap targMarkerMap;
    private final FloatArray genDist;
    private final double intensity;

    /**
     * Constructs a new {@code CurrentData} instance from the specified
     * data.
     *
     * @param par the analysis parameters
     * @param genMap the genetic map or {@code null} if no
     * genetic map is specified
     * @param data input data for the current marker window
     * @param overlapHaps haplotype constraints in the overlap with previous
     * window or {@code null} if no such constraints exist
     *
     * @throws IllegalArgumentException if
     * {@code (overlapHaps != null
     * && data.targetSamples().equals(overlapHaps.samples()) == false)}
     * @throws IllegalArgumentException if
     * {@code (overlapHaps != null &&
     * overlapHaps.marker(j).equals(data.targGT().marker(j) == false)}
     * for some {@code j} satisfying
     * {@code (0 <= j && j <= overlapHaps.nMarkers())}
     * @throws IllegalArgumentException if
     * {@code overlapHaps != null && overlapHaps.isPhased() == false}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public CurrentData(Par par, GeneticMap genMap, Data data, GT overlapHaps) {
        if (overlapHaps!=null) {
            if  (data.targetSamples().equals(overlapHaps.samples())==false) {
                throw new IllegalArgumentException("inconsistent samples");
            }
            if (overlapHaps.isPhased()==false) {
                throw new IllegalArgumentException("unphased data");
            }
        }
        this.par = par;
        this.window = data.windowIndex();
        this.prevSpliceStart = data.overlap()/2;
        this.nextOverlapStart = data.nextOverlapStart();
        this.nextSpliceStart = (data.nMarkers() + nextOverlapStart) >>> 1;
        this.prevTargSpliceStart = overlapHaps==null ? 0 : overlapHaps.nMarkers();
        this.nextTargetOverlapStart = targetIndex(data, nextOverlapStart);
        this.nextTargetSpliceStart = targetIndex(data, nextSpliceStart);

        this.ped = data.ped();
        this.targGT = (overlapHaps==null) ? data.targGT() :
                new SplicedGT(overlapHaps, data.targGT());

        this.refSamples = data.refSamples();
        this.targetSamples = data.targetSamples();
        this.allSamples = data.allSamples();
        this.nHaps = 2*allSamples.nSamples();
        this.markers = data.markers();
        this.targetMarkers = data.targetMarkers();
        this.targMarkerIndex = refToTargetMarker(data);
        this.markerIndex = targetToRefMarker(data);

        this.restRefHapPairs = data.restrictedRefHapPairs();
        this.refGT = data.refGT();
        this.restrictRefGT = data.restrictRefGT();

        this.targMarkerMap = MarkerMap.create(genMap, data.targetMarkers());
        this.genDist = genDist(targMarkerMap);
        this.intensity = 0.04*par.ne()/(2*data.nAllSamples());
    }

    private static FloatArray genDist(MarkerMap map) {
        float minCmDist = 1e-7f;
        float[] da = new float[map.nMarkers()];
        for (int j=1; j<da.length; ++j) {
            da[j] = (float) (map.genPos(j) - map.genPos(j-1));
            if (da[j] < minCmDist) {
                da[j] = minCmDist;
            }
        }
        return new FloatArray(da);
    }

    /* Returns the index of the first marker in the overlap */
    private static int nextOverlapStart(Data data, int targetOverlap) {
        if (targetOverlap < 0) {
            throw new IllegalArgumentException(String.valueOf(targetOverlap));
        }
        if (targetOverlap==0 || data.lastWindowOnChrom()) {
            return data.nMarkers();
        }
        Markers markers = data.markers();
        int nextOverlap = Math.max(0, data.nMarkers() - targetOverlap);
        while (nextOverlap>0
                && markers.marker(nextOverlap).pos()
                    == markers.marker(nextOverlap - 1).pos()) {
            --nextOverlap;
        }
        return nextOverlap;
    }

    /* Returns the index of the first marker after the next splice point */
    private static int nextSpliceStart(Data data, int overlap) {
        if (data.canAdvanceWindow() && data.lastWindowOnChrom()==false) {
            return data.nMarkers() - overlap + (overlap/2);
        }
        else {
            return data.nMarkers();
        }
    }

    /* first target index on or after specified ref index */
    private static int targetIndex(Data data, int refIndex) {
        int i=0;
        while (i<data.nTargetMarkers() && data.markerIndex(i)<refIndex) {
            ++i;
        }
        return i;
    }

    /**
     * Return the analysis parameters.
     * @return the analysis parameters
     */
    public Par par() {
        return par;
    }

    /**
     * Returns the marker window index.
     * @return the marker window index
     */
    public int window() {
        return window;
    }

    /**
     * Returns the first marker index in the overlap between this
     * marker window and the next marker window, or
     * returns {@code this.nMarkers()} the next marker window is from
     * a different chromosome.
     * @return the first marker index in the overlap between this
     * marker window and the next marker window
     */
    public int nextOverlapStart() {
        return nextOverlapStart;
    }

    /**
     * Returns the first target marker index in the overlap between this
     * marker window and the next marker window, or
     * returns {@code this.nMarkers()} if there is no overlap or if there are
     * no target markers in the overlap.
     * @return the first target marker index in the overlap between this
     * marker window and the next marker window
     */
    public int nextTargetOverlapStart() {
        return nextTargetOverlapStart;
    }

    /**
     * Returns the first marker index after the splice point with
     * the previous marker window. Returns 0 if the current marker window
     * is the first marker window.
     * @return the first marker index after the splice point with
     * the previous marker window
     */
    public int prevSpliceStart() {
        return prevSpliceStart;
    }

    /**
     * Returns the first marker index after the splice point between this
     * marker window and the next marker window, or returns
     * {@code this.nMarkers()} if there is no overlap or if there are
     * no markers after the splice point.
     * @return the first marker index after the next splice point
     */
    public int nextSpliceStart() {
        return nextSpliceStart;
    }

    /**
     * Returns the first target marker index after the splice point with
     * the previous marker window. Returns 0 if the current marker window
     * is the first marker window.
     * @return the first target marker index after the splice point with
     * the previous marker window
     */
    public int prevTargetSpliceStart() {
        return prevTargSpliceStart;
    }

    /**
     * Returns the first target marker index after the splice point between this
     * marker window and the next marker window, or returns
     * {@code this.nTargMarkers()} if there is no overlap or if there are
     * no target markers after the splice point
     * @return the first target marker index after the next splice point
     */
    public int nextTargetSpliceStart() {
        return nextTargetSpliceStart;
    }

    private int[] refToTargetMarker(Data data) {
        int[] ia = new int[data.nMarkers()];
        for (int j=0; j<ia.length; ++j) {
            ia[j] = data.targetMarkerIndex(j);
        }
        return ia;
    }

    private static int[] targetToRefMarker(Data data) {
        int[] ia = new int[data.nTargetMarkers()];
        for (int j=0; j<ia.length; ++j) {
            ia[j] = data.markerIndex(j);
        }
        return ia;
    }

    /**
     * Returns the parent-offspring relationships.
     * @return the parent-offspring relationships
     */
    public Pedigree ped() {
        return ped;
    }

    /**
     * Returns the number of reference haplotypes.
     * @return the number of reference haplotypes
     */
    public int nRefHaps() {
        return refGT==null ? 0 : refGT.nHaps();
    }

    /**
     * Returns the number of reference samples.
     * @return the number of reference samples
     */
    public int nRefSamples() {
        return refSamples == null ? 0 : refSamples.nSamples();
    }

    /**
     * Returns the list of reference samples, or {@code null} if
     * there are no reference samples.
     * @return the list of reference samples, or {@code null} if
     * there are no reference samples
     */
    public Samples refSamples() {
        return refSamples;
    }

    /**
     * Returns the number of target haplotypes.
     * @return the number of target haplotypes
     */
    public int nTargHaps() {
        return targGT.nHaps();
    }

    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */
    public int nTargetSamples() {
        return targetSamples.nSamples();
    }

    /**
     * Returns the list of target samples.
     * @return the list of target samples
     */
    public Samples targSamples() {
        return targetSamples;
    }

    public int nHaps() {
        return nHaps;
    }

    /**
     * Returns the number of reference and target samples.
     * @return the number of reference and target samples
     */
    public int nAllSamples() {
        return allSamples.nSamples();
    }

     /**
      * Returns a list of all target and reference samples.
      * Target samples are listed first in the same order as the list returned
      * by {@code this.targetSamples()}. Reference samples are listed last
      * in the same order as the list returned by {@code this.refSamples()}.
      * @return a list of all target and reference samples
      */
    public Samples allSamples() {
        return allSamples;
    }

    /**
     * Returns the number of target data markers.
     * @return the number of target data markers
     */
    public int nTargMarkers() {
        return targetMarkers.nMarkers();
    }

    /**
     * Returns the list of target data markers.
     * @return the list of target data markers
     */
    public Markers targMarkers() {
        return targetMarkers;
    }

    /**
     * Returns the number of reference data markers.
     * @return the number of reference data markers
     */
     public int nMarkers() {
         return markers.nMarkers();
     }

    /**
     * Returns the list of reference data markers.
     * @return the list of reference data markers
     */
     public Markers markers() {
         return markers;
     }

    /**
     * Returns the index of the specified marker in the reference data markers.
     * @param targetMarker index of a marker in the list of target data markers
     * @return the index of the specified marker in the reference data markers
     * @throws IndexOutOfBoundsException if
     * {@code targetMarker < 0 || targetMarker >= this.nTargMarkers()}
     */
    public int markerIndex(int targetMarker) {
        return markerIndex[targetMarker];
    }

    /**
     * Returns an array of length {@code this.nTargMarkers()} which maps
     * the {@code k}-th marker in the list of target data markers to the
     * index of the marker in the list of reference data markers.
     * @return an array of length {@code this.nTargMarkers()} which maps
     * the {@code k}-th marker in the list of target data markers to the
     * index of the marker in the list of reference data markers
     */
    public int[] markerIndices() {
        return markerIndex.clone();
    }

    /**
     * Returns the index of the specified marker in the target data, or
     * returns -1 if the marker is not present in the target data.
     * @param marker index of a marker in the reference data
     * @return the index of the specified marker in the target data, or
     * returns -1 if the marker is not present in the target data
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}.
     */
     public int targMarkerIndex(int marker) {
         return targMarkerIndex[marker];
     }

    /**
     * Returns an array of length {@code this.nMarkers()} whose {@code k}-th
     * element is the index of the {@code k}-th marker in the list of target
     * markers or is -1 if the marker is not present in the target data.
     * @return an array of length {@code this.nMarkers()} whose {@code k}-th
     * element is the index of the {@code k}-th marker in the list of target
     * markers or is -1 if the marker is not present in the target data
     */
     public int[] targMarkerIndices() {
         return targMarkerIndex.clone();
     }

    /**
     * Returns a list with the specified haplotypes following by the
     * reference haplotype pairs that are restricted to the target data markers.
     * @param list a list of haplotype pairs for target data markers
     * @return a list with the specified haplotypes following by the
     * reference haplotype pairs that are restricted to the target data markers
     * @throws NullPointerException if {@code list == null}
     */
    public List<HapPair> addRestrictedRefHapPairs(List<BitHapPair> list) {
        List<HapPair> newList = new ArrayList<HapPair>(list);
        newList.addAll(restRefHapPairs);
        return newList;
    }

    /**
     * Returns the phased, nonmissing reference genotype data
     * or {@code null} if there are no reference data.
     * @return the reference genotype data or {@code null} if there are no
     * reference data
     */
    public RefGT refGT() {
        return refGT;
    }

     /**
     * Returns the phased, nonmissing reference genotype data for
     * the target data markers or {@code null} if there are no reference data.
     * @return the reference genotype data for the target data markers or
     * {@code null} if there are no reference data
     */
    public RefGT restrictRefGT() {
        return restrictRefGT;
    }

    /**
     * Returns the genotype likelihoods for the
     * target samples at the target data markers.
     * @return the genotype likelihoods for the
     * target samples at the target data markers.
     */
    public GT targGT() {
        return targGT;
    }

    /**
     * Returns the intensity used to generate the pRecomb values.
     * @return the intensity used to generate the pRecomb values
     */
    public float intensity() {
        return (float) intensity;
    }

    /**
     * Returns the genetic map.
     * @return the genetic map
     */
    public MarkerMap map() {
        return targMarkerMap;
    }

    /**
     * Return a {@code FloatArray} of size {@code this.nTargMarkers()}
     * whose {@code k}-th element is the genetic distance between the
     * {@code k}-th target marker and the previous marker, or {@code 0.0}
     * if {@code (k == 0)}.
     * @return a {@code FloatArray} of size {@code this.nTargMarkers()}
     * whose {@code k}-th element is the genetic distance between the
     * {@code k}-th target marker and the previous marker,
     */
    public FloatArray genDist() {
        return genDist;
    }
}
