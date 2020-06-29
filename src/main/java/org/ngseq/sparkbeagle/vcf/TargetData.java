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

import org.ngseq.sparkbeagle.Par;
import org.ngseq.sparkbeagle.Pedigree;
import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.blbutil.SampleFileIt;
import org.ngseq.sparkbeagle.haplotype.HapPair;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

/**
 * <p>Class {@code TargetData} represents a sliding window of
 * target VCF records.
 * </p>
 * <p>Instances of class {@code TargetData} are not thread-safe.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class TargetData implements Data {

    private final Pedigree ped;
    private int window = 0;
    private Window<GTRec> currentWindow;
    private Markers markers;
    private GTRec[] emData;
    private GT gt;

    private final WindowIt<GTRec> targWindow;
    private int cumMarkerCnt = 0;

    /**
     * Constructs and returns a new {@code TargetData} instance from
     * VcfRecords returned by the specified {@code SampleFileIt} objects.
     *
     * @param supplier a supplier for the sample file iterator
     * @param par the command line parameters
     * @return a new {@code TargetData} instance
     *
     * @throws IllegalArgumentException if the data returned by
     * the specified iterator contains no samples
     * @throws IllegalArgumentException if a format error is detected
     * in a string VCF record
     * @throws IllegalArgumentException if
     * {@code overlap < 0 || Float.isFinite(overlap) == false}
     * @throws IllegalArgumentException if
     * {@code window <= overlap || Float.isFinite(window) == false}
     * @throws NullPointerException if
     * {@code it == null || ped == null || genMap == null}
     */
    public static TargetData targetData(Par par,
            Supplier<SampleFileIt<GTRec>> supplier) {
        GeneticMap genMap = GeneticMap.geneticMap(par.map(), par.chromInt());
        WindowIt<GTRec> targWindow = WindowIt.newInstance(supplier, genMap,
                par.window(), par.overlap());
        TargetData targetData = new TargetData(par, targWindow);
        assert targetData.canAdvanceWindow();
        targetData.advanceWindowCm();
        return targetData;
    }

    private TargetData(Par par, WindowIt<GTRec> targWindow) {
        this.ped = new Pedigree(targWindow.samples(), par.ped());
        this.targWindow = targWindow;
        this.currentWindow = null;
        this.markers = Markers.create(new Marker[0]);
        this.emData = new GTRec[0];
        this.gt = targGT(targWindow.samples(), emData, ped);
    }

    private static GT targGT(Samples samples, GTRec[] targData,
            Pedigree ped) {
        GT gt = new BasicGT(samples, targData);
        return gt.isGTData() ? new XBasicGT(gt, ped) : gt;
    }

    @Override
    public Pedigree ped() {
        return ped;
    }

    @Override
    public GeneticMap genMap() {
        return targWindow.genMap();
    }

    @Override
    public boolean lastWindowOnChrom() {
        return currentWindow.lastWindowOnChrom();
    }

    @Override
    public boolean canAdvanceWindow() {
        return targWindow.hasNext();
    }

    @Override
    public void advanceWindowCm() {
        currentWindow = targWindow.next();
        cumMarkerCnt += (currentWindow.nMarkers() - currentWindow.overlapEnd());
        emData = currentWindow.recList().toArray(new GTRec[0]);
        markers = extractMarkers(emData);
        gt = targGT(targWindow.samples(), emData, ped);
        ++window;
    }

    private static Markers extractMarkers(GTRec[] markerData) {
        Marker[] ma = new Marker[markerData.length];
        for (int j=0; j<ma.length; ++j) {
            ma[j] = markerData[j].marker();
        }
        return Markers.create(ma);
    }

    @Override
    public int windowIndex() {
        return window;
    }


    @Override
    public int targetOverlap() {
        return currentWindow.overlapEnd();
    }

    @Override
    public int overlap() {
        return currentWindow.overlapEnd();
    }

    @Override
    public int nextOverlapStart() {
        return currentWindow.overlapStart();
    }

    @Override
    public int nTargetMarkers() {
        return markers.nMarkers();
    }

    @Override
    public int nTargetMarkersSoFar() {
        return cumMarkerCnt;
    }

    @Override
    public Markers targetMarkers() {
        return markers;
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public int nMarkersSoFar() {
        return cumMarkerCnt;
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public int targetMarkerIndex(int refIndex) {
        if (refIndex < 0 || refIndex >= markers.nMarkers()) {
            throw new ArrayIndexOutOfBoundsException(refIndex);
        }
        return refIndex;
    }

    @Override
    public int markerIndex(int nonRefIndex) {
        if (nonRefIndex < 0 || nonRefIndex >= markers.nMarkers()) {
            throw new ArrayIndexOutOfBoundsException(nonRefIndex);
        }
        return nonRefIndex;
    }

    @Override
    public int nTargetSamples() {
        return targWindow.samples().nSamples();
    }

    @Override
    public Samples targetSamples() {
       return targWindow.samples();
    }

    @Override
    public int nRefSamples() {
        return 0;
    }

    @Override
    public Samples refSamples() {
        return null;
    }

    @Override
    public int nAllSamples() {
        return nTargetSamples();
    }

    @Override
    public Samples allSamples() {
        return targetSamples();
    }

    @Override
    public GT targGT() {
        return gt;
    }

    @Override
    public List<HapPair> restrictedRefHapPairs() {
        // no reference haplotypes to add
        return new ArrayList<>();
    }

    @Override
    public List<HapPair> refHapPairs() {
        // no reference haplotypes to return
        return new ArrayList<>();
    }


    @Override
    public RefGT refGT() {
        // no reference genotype data to return
        return null;
    }

    @Override
    public RefGT restrictRefGT() {
        // no reference genotype data to return
        return null;
    }

    @Override
    public void close() {
       targWindow.close();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("vcf.NonRefData");
        return sb.toString();
    }
}
