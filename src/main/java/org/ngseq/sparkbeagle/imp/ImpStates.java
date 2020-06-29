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

import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.ints.IntMap;

/**
 * <p>Class {@code ImpStates} identifies a list of pseudo-reference haplotypes
 * for a target haplotype. Each pseudo-reference haplotype is a
 * one-dimensional mosaic of reference haplotype segments.
 * </p>
 * <p>Instances of {@code ImpStates} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ImpStates {

    private final ImpIbs ibsHaps;
    private final int nStates;

    private final IntMap<CopyHap> hap2Slots;
    private final CopyHap[] slotHeap;
    private final IntList[] copyHaps;

    /**
     * Constructs a new {@code ImpStates} object from the specified data.
     * @param ibsHaps the IBS haplotype segments
     * @throws NullPointerException if
     * {@code impData == null || ibsStates == null}
     */
    public ImpStates(ImpIbs ibsHaps) {
        this.ibsHaps = ibsHaps;
        this.nStates = ibsHaps.nStates();
        this.hap2Slots = new IntMap<>(nStates);
        this.slotHeap = new CopyHap[nStates];
        this.copyHaps = new IntList[nStates];
        for (int j=0; j<nStates; ++j) {
            slotHeap[j] = new CopyHap(j);
            copyHaps[j] = new IntList(10);
        }
    }

    /**
     * Returns the number of HMM states per marker.
     * @return the number of HMM states per marker
     */
    public int nStates() {
        return nStates;
    }

    /**
     * Stores the reference haplotype for the {@code j}-th state
     * at the {@code m}-th marker in {@code hapIndices[m][j]}, and stores
     * the equality of the allele carried by the reference haplotype for
     * the {@code j}-th state and the allele carried by the target haplotype
     * at the {@code m}-th marker in {@code alMatch[m][j]}.  The number of
     * HMM states states at each marker is returned.
     * @param targHap the haplotype index
     * @param hapIndices the two-dimensional array in which
     * reference haplotype indices for each HMM state will be stored
     * @param alMatch the two-dimensional array in which allele match status
     * between the target haplotype and HMM state will be stored
     * @return the number of HMM states at each marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code targHap < 0 || targHap >= this.nTargHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code (hapIndices.length < this.impData().nClusters())} or if any
     * {@code hapIndices[j].length} is less than the number of HMM states
     * at a marker for any {@code j} satisfying
     * {@code (0 <= j && j < this.nClusters())}
     * @throws IndexOutOfBoundsException if
     * {@code (alMatch.length < this.impdata().nClusters())} or if any
     * {@code alMatch[j].length} is less than the number of HMM states
     * at a marker for any {@code j} satisfying
     * {@code (0 <= j && j < this.nClusters())}
     * @throws NullPointerException if
     * {@code (hapIndices == null)} or if {@code (hapIndices[j] == null)} for
     * any {@code j} satisfying {@code (0 <= j && j < this.nClusters())}
     * @throws NullPointerException if
     * {@code (alMatch == null)} or if {@code (alMatch[j] == null)} for
     * any {@code j} satisfying {@code (0 <= j && j < this.nClusters())}
     */
    public int ibsStates(int targHap, int[][] hapIndices, boolean[][] alMatch) {
        initializeFields();
        for (int j=0, n=ibsHaps.nSteps(); j<n; ++j) {
            int[] ibs = ibsHaps.ibsHaps(targHap, j);
            for (int hap : ibs) {
                updateFields(hap, j);
            }
        }
        int numStates = copyData(targHap, hapIndices, alMatch);
        return numStates;
    }

    private void initializeFields() {
        hap2Slots.clear();
        for (int j=0; j<nStates; ++j) {
            slotHeap[j].hap = -1;
            slotHeap[j].end = -1;
            slotHeap[j].slotIndex = slotHeap[j].heapIndex;
            copyHaps[j].clear();
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
            int slot = slotHeap[0].slotIndex;
            if (prevEnd >= 0) {
                hap2Slots.remove(prevHap);
                copyHaps[slot].add( ibsHaps.stepStart((prevEnd + end) >>> 1) );
            }
            copyHaps[slot].add(hap);
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

    private int copyData(int targHap, int[][] hapIndices, boolean[][] alMatch) {
        ImpData impData = ibsHaps.impData();
        int stateCnt = 0;
        for (int j=0; j<nStates; ++j) {
            int start = 0;
            IntList il = copyHaps[j];
            il.add(ibsHaps.impData().nClusters());
            if (il.size()>1) {
                for (int k=0, n=il.size(); k<n; k+=2) {
                    int hap = il.get(k);
                    int end = il.get(k+1);
                    for (int m=start; m<end; ++m) {
                        hapIndices[m][stateCnt] = hap;
                    }
                    start = end;
                }
                stateCnt++;
            }
        }
        int shiftedTargHap = impData.nRefHaps() + targHap;
        for (int m=0; m<hapIndices.length; ++m) {
            int targAllele = impData.allele(m, shiftedTargHap);
            for (int k=0; k<stateCnt; ++k) {
                int refHap = hapIndices[m][k];
                alMatch[m][k] = impData.allele(m, refHap)==targAllele;
            }
        }
        return stateCnt;
    }

    private static class CopyHap {
        private int hap;
        private int end;
        private int heapIndex;
        private int slotIndex;

        public CopyHap(int heapIndex) {
            this.hap = -1;
            this.end = -1;
            this.heapIndex = heapIndex;
            this.slotIndex = heapIndex;
        }
    }
}
