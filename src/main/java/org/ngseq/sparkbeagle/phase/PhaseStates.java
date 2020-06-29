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
package org.ngseq.sparkbeagle.phase;

import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.ints.IntMap;
import org.ngseq.sparkbeagle.ints.IntSet;

/**
 * <p>Class {@code PhaseStates} identifies a rolling window of reference
 * haplotypes for a target sample.
 * </p>
 * <p>Instances of {@code PhaseStates} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class PhaseStates {

    private final PhaseIbs ibsHaps;
    private final int nStates;

    private final IntMap<CopyHap> hap2Slots;
    private final CopyHap[] slotHeap;
    private final IntList[] copyHaps;
    private final IntList[] copyEnds;

    /**
     * Constructs a new {@code PhaseIbs} object from the specified data.
     * @param ibsHaps IBS haplotype segments
     * @throws IllegalArgumentException if {@code nHapsPerStep < 1}
     * @throws NullPointerException if {@code ibsStates == null}
     */
    public PhaseStates(PhaseIbs ibsHaps) {
        this.ibsHaps = ibsHaps;
        this.nStates = ibsHaps.nStates();
        this.hap2Slots = new IntMap<>(nStates);
        this.slotHeap = new CopyHap[nStates];
        this.copyHaps = new IntList[nStates];
        this.copyEnds = new IntList[nStates];
        for (int j=0; j<nStates; ++j) {
            slotHeap[j] = new CopyHap(j);
            copyHaps[j] = new IntList(8);
            copyEnds[j] = new IntList(8);
        }
    }

    /**
     * Returns the number of HMM states per marker.
     * @return the number of HMM states per marker
     */
    public int nStates() {
        return ibsHaps.nStates();
    }

    /**
     * Identifies the HMM state alleles for the specified sample.  The
     * {@code j}-th state allele for the {@code m}-th marker will be
     * stored in {@code stateAlleles[m][j]}.
     * @param sample the sample index
     * @param stateAlleles the two-dimensional array in which
     * state alleles will be stored
     * @return the number of state alleles at each marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.hapPairs().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code stateAlleles.length < this.hapPairs().nMarkers()} or if
     * {@code stateAlleles[j].length < this.nStates()} for any
     * {@code j} satisfying {@code (0 <= j && j < this.hapPairs().nMarkers())}
     * @throws NullPointerException if
     * {@code stateAlleles == null} or if {@code stateAlleles[j] == null} for
     * any {@code j} satisfying
     * {@code (0 <= j && j < this.hapPairs().nMarkers())}
     */
    public int ibsStates(int sample,  int[][] stateAlleles) {
        int h1 = 2*sample;
        int h2 = 2*sample + 1;

        initializeFields();
        for (int w=0, n=ibsHaps.nSteps(); w<n; ++w) {
            IntSet intSet = ibsHaps.ibsHaps(h1, h2, w);
            for (int j=0, m=intSet.size(); j<m; ++j) {
                updateFields(intSet.elementWithIndex(j), w);
            }
            intSet = ibsHaps.ibsHaps(h2, h1, w);
            for (int j=0, m=intSet.size(); j<m; ++j) {
                updateFields(intSet.elementWithIndex(j), w);
            }
        }
        int numStates = copyData(stateAlleles);
        if (numStates<2) {
            numStates = naiveStates(sample, stateAlleles);
        }
        return numStates;
    }

    private void initializeFields() {
        hap2Slots.clear();
        for (int j=0; j<nStates; ++j) {
            slotHeap[j].hap = -1;
            slotHeap[j].end = -1;
            slotHeap[j].copyIndex = slotHeap[j].heapIndex;
            copyHaps[j].clear();
            copyEnds[j].clear();
        }
    }

    private void updateFields(int hap, int end) {
        CopyHap data = hap2Slots.get(hap);
        if (data!=null) {
            data.end = end;
            heapifyDown(slotHeap, data.heapIndex);
        }
        else {
            int prevHap = slotHeap[0].hap;
            int prevEnd = slotHeap[0].end;
            int copy = slotHeap[0].copyIndex;
            if (prevEnd >= 0) {
                hap2Slots.remove(prevHap);
                copyEnds[copy].add( ibsHaps.stepStart((prevEnd + end) >> 1) );
            }
            copyHaps[copy].add(hap);
            hap2Slots.put(hap, slotHeap[0]);
            slotHeap[0].hap = hap;
            slotHeap[0].end = end;
            heapifyDown(slotHeap, 0);
        }
    }

    private static void heapifyDown(CopyHap[] heap, int index) {
        CopyHap tmp = heap[index];
        int child = 2*index + 1;
        while (child < heap.length) {
            if ((child+1<heap.length) && (heap[child+1].end<heap[child].end)) {
                child++;
            }
            if (heap[child].end < tmp.end) {
                heap[index] = heap[child];
                heap[index].heapIndex = index;
                index = child;
                child = 2*index + 1;
            }
            else {
                break;
            }
        }
        heap[index] = tmp;
        heap[index].heapIndex = index;
    }

    private int copyData(int[][] stateAlleles) {
        PhaseData phaseData = ibsHaps.phaseData();
        int nMarkers = phaseData.nMarkers();
        IntList nonEmptyList = new IntList(nStates);
        for (int j=0; j<nStates; ++j) {
            if (copyHaps[j].size()>0) {
                copyEnds[j].add(nMarkers);
                nonEmptyList.add(j);
            }
        }
        int[] usedIndex = nonEmptyList.toArray();
        int[] cpIndex = new int[usedIndex.length];
        int[] hap = new int[usedIndex.length];
        int[] end = new int[usedIndex.length];
        for (int j=0; j<usedIndex.length; ++j) {
            hap[j] = copyHaps[usedIndex[j]].get(0);
            end[j] = copyEnds[usedIndex[j]].get(0);
        }
        for (int m=0; m<nMarkers; ++m) {
            for (int j=0; j<usedIndex.length; ++j) {
                if (m==end[j]) {
                    ++cpIndex[j];
                    hap[j] = copyHaps[usedIndex[j]].get(cpIndex[j]);
                    end[j] = copyEnds[usedIndex[j]].get(cpIndex[j]);
                }
                stateAlleles[m][j] = phaseData.allele(m, hap[j]);
            }
        }
        return usedIndex.length;
    }

    private int naiveStates(int sample, int[][] stateAlleles) {
        PhaseData phaseData = ibsHaps.phaseData();
        int nMarkers = phaseData.nMarkers();
        int nHaps = phaseData.nHaps();
        int maxNStates = Math.min(nStates, nHaps-2);
        int h1 = 2*sample;
        int h2 = h1 + 1;
        int hap = h2;     // start coping immediately after h2
        int stateCnt = 0;
        for (int j=0; j<maxNStates; ++j) {
            if (++hap >= nHaps) {
                hap -= nHaps;
            }
            if (hap!=h1 && hap!=h2) {
                for (int m=0; m<nMarkers; ++m) {
                    stateAlleles[m][stateCnt] = phaseData.allele(m, hap);
                }
                ++stateCnt;
            }
        }
        return stateCnt;
    }

    private static class CopyHap {
        private int hap;
        private int end;
        private int heapIndex;
        private int copyIndex;

        public CopyHap(int heapIndex) {
            this.hap = -1;
            this.end = -1;
            this.heapIndex = heapIndex;
            this.copyIndex = heapIndex;
        }
    }
}
