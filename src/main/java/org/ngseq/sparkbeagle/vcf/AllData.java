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
import org.ngseq.sparkbeagle.haplotype.WrappedHapPair;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

/**
 * <p>Class {@code AllData} represents a sliding window of
 reference and target VCF recList.
 </p>
 * <p>Instances of class {@code AllData} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AllData implements Data {

    private final Pedigree ped;
    private int window = 0;
    private Window<RefGTRec> currentRefWindow;
    private RefGT refGT;
    private RefGT restrictRefGT;
    private GTRec[] targData;  // missing markers as null entries
    private int[] refIndices;
    private int[] targIndices;
    private GT targGTWindow;
    private final Samples allSamples;

    private final List<HapPair> refHapPairs;
    private final List<HapPair> targRefHapPairs; // at target markers
    private final WindowIt<RefGTRec> refWindow;
    private final RestrictedVcfWindow targWindow;

    private int cumMarkerCnt = 0;

    /**
     * Constructs and returns a new {@code AllData} instance from VCF
     * recList returned by the specified {@code SampleFileIt} objects.
     *
     * @param supplier an object to supply the reference file iterator
     * @param targIt an iterator that returns target VCF recList
     * @param par the command line parameters
     * @return a new {@code AllData} instance
     *
     * @throws IllegalArgumentException if either the reference data or
     * target data contain no samples
     * @throws IllegalArgumentException if a format error is detected
     * in a string VCF record
     * @throws IllegalArgumentException if
     * {@code overlap < 0 || Float.isFinite(overlap) == false}
     * @throws IllegalArgumentException if
     * {@code window <= overlap || Float.isFinite(window) == false}
     * @throws NullPointerException if
     * {@code refIt == null || targetIt == null || genMap == null}
     */
    public static AllData allData(Supplier<SampleFileIt<RefGTRec>> supplier,
            SampleFileIt<GTRec> targIt, Par par) {
        GeneticMap genMap = GeneticMap.geneticMap(par.map(), par.chromInt());
        WindowIt<RefGTRec> refWindow = WindowIt.newInstance(supplier, genMap,
                par.window(), par.overlap());
        if (refWindow.samples().nSamples()==0 || targIt.samples().nSamples()==0) {
            throw new IllegalArgumentException("nSamples==0");
        }
        RestrictedVcfWindow targetWindow = new RestrictedVcfWindow(targIt);
        AllData allData = new AllData(par, refWindow, targetWindow);
        assert allData.canAdvanceWindow();
        allData.advanceWindowCm();
        return allData;
    }

    private AllData(Par par, WindowIt<RefGTRec> refWind,
            RestrictedVcfWindow targWind) {
        this.ped = new Pedigree(targWind.samples(), par.ped());
        this.refWindow = refWind;
        this.targWindow = targWind;

        this.currentRefWindow = null;
        this.refGT = null;
        this.restrictRefGT = null;
        this.targData = new GTRec[0];
        this.refIndices = new int[0];
        this.targIndices = new int[0];
        this.targGTWindow = targGTWindow(targWind.samples(), targData, ped);
        this.allSamples = allSamples(refWind.samples(), targWind.samples());

        this.refHapPairs = new ArrayList<>(0);
        this.targRefHapPairs = new ArrayList<>(0);
    }

    private static Samples allSamples(Samples ref, Samples target) {
        /*
           Target samples are listed first so that sample indices agree
           with sample indices in target data genotype likelihoods.
        */
        int nRef = ref.nSamples();
        int nTarget = target.nSamples();
        int[] idIndices = new int[nRef + nTarget];
        for (int j=0; j<nTarget; ++j) {
            idIndices[j] = target.idIndex(j);
        }
        for (int j=0; j<nRef; ++j) {
            idIndices[nTarget + j] = ref.idIndex(j);
        }
        return new Samples(idIndices);
    }

    private static GT targGTWindow(Samples samples, GTRec[] targData,
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
        return refWindow.genMap();
    }

    @Override
    public boolean lastWindowOnChrom() {
        return currentRefWindow.lastWindowOnChrom();
    }

    @Override
    public boolean canAdvanceWindow() {
        return refWindow.hasNext();
    }

    @Override
    public void advanceWindowCm() {
        currentRefWindow = refWindow.next();
        cumMarkerCnt += (currentRefWindow.nMarkers() - currentRefWindow.overlapEnd());
        RefGTRec[] refRecords = currentRefWindow.recList().toArray(new RefGTRec[0]);
        refGT = new BasicRefGT(refRecords);
        targData = targWindow.advanceWindow(refGT.markers());
        refIndices = refIndices(targData);
        targIndices = targetIndices(targData);
        targGTWindow = targGTWindow(targWindow.samples(), targData, refIndices, ped);
        ++window;
        setRefHaplotypes(refGT);
        setTargRefData(targGTWindow.markers(), refRecords, refIndices);
    }

    @Override
    public int windowIndex() {
        return window;
    }

    private static int[] refIndices(GTRec[] vma) {
        int nonNullCnt = 0;
        for (GTRec vm : vma) {
            if (vm!=null) {
                ++nonNullCnt;
            }
        }
        int[] inclusionMap = new int[nonNullCnt];
        int index = 0;
        for (int j=0; j<vma.length; ++j) {
            if (vma[j]!=null) {
                inclusionMap[index++] = j;
            }
        }
        if (index != inclusionMap.length) {
            throw new IllegalStateException("vma modification detected");
        }
        return inclusionMap;
    }

    private static int[] targetIndices(GTRec[] vma) {
        int[] inclusionMap = new int[vma.length];
        int index = 0;
        for (int j=0; j<inclusionMap.length; ++j) {
            if (vma[j]!=null) {
                inclusionMap[j] = index++;
            }
            else {
                inclusionMap[j] = -1;
            }
        }
        return inclusionMap;
    }

    private static GT targGTWindow(Samples samples, GTRec[] vma,
            int[] refMarkerIndex, Pedigree ped) {
        GTRec[] restricted = new GTRec[refMarkerIndex.length];
        for (int j=0; j<refMarkerIndex.length; ++j) {
            restricted[j] = vma[refMarkerIndex[j]];
        }
        GT gt = new BasicGT(samples, restricted);
        return gt.isGTData() ? new XBasicGT(gt, ped) : gt;
    }

    private void setRefHaplotypes(RefGT refGT) {
        refHapPairs.clear();
        for (int j=0, n=refGT.nSamples(); j<n; ++j) {
            refHapPairs.add(new WrappedHapPair(refGT, j));
        }
    }

    private void setTargRefData(Markers targMarkers,
            RefGTRec[] refData, int[] refMarkerIndices) {
        assert targMarkers.nMarkers()==refMarkerIndices.length;
        targRefHapPairs.clear();
        RefGTRec[] vma = new RefGTRec[refMarkerIndices.length];
        for (int j=0; j<refMarkerIndices.length; ++j) {
            vma[j] = refData[refMarkerIndices[j]];
        }
        restrictRefGT = new BasicRefGT(targMarkers, refWindow.samples(), vma);
        for (int j=0, n=restrictRefGT.nSamples(); j<n; ++j) {
            targRefHapPairs.add(new WrappedHapPair(restrictRefGT, j));
        }
    }

    @Override
    public int targetOverlap() {
        return targWindow.overlap();
    }

    @Override
    public int overlap() {
        return currentRefWindow.overlapEnd();
    }

    @Override
    public int nextOverlapStart() {
        return currentRefWindow.overlapStart();
    }

    @Override
    public int nTargetMarkers() {
        return targGTWindow.markers().nMarkers();
    }

    @Override
    public int nTargetMarkersSoFar() {
        return targWindow.cumMarkerCnt();
    }

    @Override
    public Markers targetMarkers() {
        return targGTWindow.markers();
    }


    @Override
    public int nMarkers() {
        return refGT.nMarkers();
    }

    @Override
    public int nMarkersSoFar() {
        return cumMarkerCnt;
    }

    @Override
    public Markers markers() {
        return refGT.markers();
    }

    @Override
    public int targetMarkerIndex(int refIndex) {
        return targIndices[refIndex];
    }

    @Override
    public int markerIndex(int nonRefIndex) {
        return refIndices[nonRefIndex];
    }

    @Override
    public int nTargetSamples() {
        return targGTWindow.nSamples();
    }

    @Override
    public Samples targetSamples() {
        return targGTWindow.samples();
    }

    @Override
    public int nRefSamples() {
        return refWindow.samples().nSamples();
    }

    @Override
    public Samples refSamples() {
        return refWindow.samples();
    }

    @Override
    public int nAllSamples() {
        return allSamples.nSamples();
    }

    @Override
    public Samples allSamples() {
        return allSamples;
    }

    @Override
    public GT targGT() {
       return targGTWindow;
    }

    @Override
    public List<HapPair> restrictedRefHapPairs() {
        return new ArrayList<>(targRefHapPairs);
    }

    @Override
    public List<HapPair> refHapPairs() {
        return new ArrayList<>(refHapPairs);
    }

    @Override
    public RefGT refGT() {
        return refGT;
    }

    @Override
    public RefGT restrictRefGT() {
        return restrictRefGT;
    }

    @Override
    public void close() {
        refWindow.close();
        targWindow.close();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(20);
        sb.append(this.getClass().toString());
        return sb.toString();
    }
}
