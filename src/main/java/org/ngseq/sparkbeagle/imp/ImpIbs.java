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

import org.ngseq.sparkbeagle.Par;
import org.ngseq.sparkbeagle.blbutil.Utilities;
import org.ngseq.sparkbeagle.ints.IndexArray;
import org.ngseq.sparkbeagle.ints.IntArray;
import org.ngseq.sparkbeagle.ints.IntList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * <p>Class {@code ImpIbs} identifies haplotypes that share a long
 * IBS segment with a specified haplotype.
 * </p>
 * <p>Instances of {@code ImpIbs} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ImpIbs {

    private final ImpData impData;
    private final int nRefHaps;
    private final long seed;
    private final int nStates;
    private final int nSteps;
    private final int nHapsPerStep;

    private final int[] stepStarts;
    private final int[][][] ibsHaps; //[window][targ hap][ibs_set]

    /**
     * Constructs a new {@code ImpIbs} object from the specified data.
     * @param impData the input data for genotype imputation
     *
     * @throws NullPointerException if {@code impData == null}
     */
    public ImpIbs(ImpData impData) {
        Par par = impData.par();
        this.impData = impData;
        this.seed = par.seed();
        this.nRefHaps = impData.nRefHaps();
        this.nStates = par.imp_states();
        this.nSteps = par.nsteps();
        int nStepsPerSegment = Math.round(par.imp_segment()/par.step());
        this.nHapsPerStep = (par.imp_states() / nStepsPerSegment);

        this.stepStarts = stepStarts(impData);
        IndexArray[] codedSteps = IntStream.range(0, stepStarts.length)
                .parallel()
                .mapToObj(j -> codeStep(impData, stepStarts, j))
                .toArray(IndexArray[]::new);
        this.ibsHaps = IntStream.range(0, codedSteps.length)
                .parallel()
                .mapToObj(j -> getIbsHaps(codedSteps, j))
                .toArray(int[][][]::new);
    }

    private static int[] stepStarts(ImpData impData) {
        double[] pos = impData.pos();
        double step = impData.par().step();
        IntList indices = new IntList(pos.length/10);
        indices.add(0);
        double nextPos =  pos[0] + step/2;  // make first step be half-length
        int index = nextIndex(pos, 0, nextPos);
        while (index < pos.length) {
            indices.add(index);
            nextPos = pos[index] + step;
            index = nextIndex(pos, index, nextPos);
        }
        return indices.toArray();
    }

    private static int nextIndex(double[] pos, int start, double targetPos) {
        int nextIndex = Arrays.binarySearch(pos, start, pos.length, targetPos);
        return (nextIndex<0) ? -nextIndex-1 : nextIndex;
    }

    private static IndexArray codeStep(ImpData impData, int[] starts, int startIndex) {
        int nRefHaps = impData.nRefHaps();
        int nHaps = impData.nHaps();
        int[] hapToSeq = IntStream.range(0, nHaps).map(i -> 1).toArray();
        int start = starts[startIndex];
        int end = (startIndex+1) < starts.length ?  starts[startIndex+1] : impData.nClusters();

        int seqCnt = 2; // seq 0 is reserved for sequences not found in target
        for (int m=start; m<end; ++m) {
            IndexArray h2s = impData.hapToSeq(m);
            IntArray codedHaps = h2s.intArray();
            int nAlleles = h2s.valueSize();
            int[] seqMap = new int[seqCnt*nAlleles];

            seqCnt = 1;
            for (int h=nRefHaps; h<nHaps; ++h) {
                int index = nAlleles*hapToSeq[h] + codedHaps.get(h);
                if (seqMap[index] == 0) {
                    seqMap[index] = seqCnt++;
                }
                hapToSeq[h] = seqMap[index];
            }
            for (int h=0; h<nRefHaps; ++h) {
                if (hapToSeq[h] != 0) {
                    int index = hapToSeq[h]*nAlleles + codedHaps.get(h);
                    hapToSeq[h] = seqMap[index];
                }
            }
        }
        IntArray intArray = IntArray.create(hapToSeq, seqCnt);
        return new IndexArray(intArray, seqCnt);
    }

    private int[][] getIbsHaps(IndexArray[] codedSteps, int index) {
        int nTargHaps = impData.nTargHaps();
        int[][] results = new int[nTargHaps][];
        int nStepsToMerge = Math.min(nSteps, codedSteps.length - index);
        List<IntList> children = initPartition(codedSteps[index]);
        List<IntList> lastIbs = new ArrayList<>(children.size());
        initUpdateResults(children, results, lastIbs);
        for (int i=1; i<nStepsToMerge; ++i) {
            int initCapacity = Math.min(nTargHaps, 2*lastIbs.size());
            List<IntList> nextIbs = new ArrayList<>(initCapacity);
            IndexArray codedStep = codedSteps[index+i];
            for (int j=0, nj=lastIbs.size(); j<nj; ++j) {
                IntList parent = lastIbs.get(j);
                children = partition(parent, codedStep);
                updateResults(parent, children, results, nextIbs);
            }
            lastIbs = nextIbs;
        }
        finalUpdateResults(lastIbs, results);
        return results;
    }

    private List<IntList> initPartition(IndexArray codedStep) {
        IntList[] list = new IntList[codedStep.valueSize()];
        IntArray hap2Seq = codedStep.intArray();
        int nHaps = hap2Seq.size();
        List<IntList> children = new ArrayList<>();
        for (int h=nRefHaps; h<nHaps; ++h) {
            int seq = hap2Seq.get(h);
            if (list[seq]==null) {
                list[seq] = new IntList();
                children.add(list[seq]);
            }
        }
        for (int h=0; h<nHaps; ++h) {
            int seq = hap2Seq.get(h);
            if (list[seq]!=null) {
                list[seq].add(h);
            }
        }
        return children;
    }

    private List<IntList> partition(IntList parent, IndexArray codedStep) {
        IntList[] list = new IntList[codedStep.valueSize()];
        IntArray hap2Seq = codedStep.intArray();
        int nHaps = parent.size();
        List<IntList> children = new ArrayList<>();
        int targStart = insertionPoint(parent, nRefHaps);
        for (int k=targStart; k<nHaps; ++k) {
            int hap = parent.get(k);
            int seq = hap2Seq.get(hap);
            if (list[seq]==null) {
                list[seq] = new IntList();
                children.add(list[seq]);
            }
        }
        for (int k=0; k<nHaps; ++k) {
            int hap = parent.get(k);
            int seq = hap2Seq.get(hap);
            if (list[seq]!=null) {
                list[seq].add(hap);
            }
        }
        return children;
    }

    private void initUpdateResults(List<IntList> lists, int[][] result,
            List<IntList> nextIbs) {
        for (int j=0, n=lists.size(); j<n; ++j) {
            IntList hapList = lists.get(j);
            int nRef = insertionPoint(hapList, nRefHaps);
            if (nRef <= nHapsPerStep) {
                setResult(hapList, nRef, result, hapList.copyOf(nRef));
            }
            else {
                nextIbs.add(hapList);
            }
        }
    }

    private void updateResults(IntList parent, List<IntList> children,
            int[][] results, List<IntList> nextIbs) {
        for (int k=0, nk=children.size(); k<nk; ++k) {
            IntList child = children.get(k);
            int nRef = insertionPoint(child, nRefHaps);
            if (nRef <= nHapsPerStep) {
                int[] ibsList = combineLists(parent, child, nRef);
                setResult(child, nRef, results, ibsList);
                child.clear();
            }
            else {
                nextIbs.add(child);
            }
        }
    }

    private void finalUpdateResults(List<IntList> children, int[][] results) {
        for (int j=0, n=children.size(); j<n; ++j) {
            IntList child = children.get(j);
            int nRef = insertionPoint(child, nRefHaps);
            int[] ibsList = child.copyOf(nRef);
            if (nHapsPerStep < ibsList.length) {
                Random rand = new Random(seed + child.get(0));
                Utilities.shuffle(ibsList, rand);
                ibsList = Arrays.copyOf(ibsList, nHapsPerStep);
                Arrays.sort(ibsList);
            }
            setResult(child, nRef, results, ibsList);
        }
    }

    private int[] combineLists(IntList parent, IntList child, int nChildRef) {
        int[] ibsList = child.copyOf(nHapsPerStep);
        if (nChildRef<nHapsPerStep) {
            Random rand = new Random(seed + parent.get(0) + child.get(0));
            int nParentRef = insertionPoint(parent, nRefHaps);
            int[] parentRefHaps = parent.copyOf(nParentRef);
            Utilities.shuffle(parentRefHaps, rand);
            for (int i=0, ibsIndex=nChildRef; ibsIndex<nHapsPerStep; ++i) {
                if (Arrays.binarySearch(ibsList, 0, nChildRef, parentRefHaps[i]) < 0) {
                    ibsList[ibsIndex++] = parentRefHaps[i];
                }
            }
        }
        Arrays.sort(ibsList);
        return ibsList;
    }

    private void setResult(IntList hapIndices, int firstTargIndex,
            int[][] result, int[] ibsHaps) {
        for (int j=firstTargIndex, nj=hapIndices.size(); j<nj; ++j) {
            result[hapIndices.get(j) - nRefHaps] = ibsHaps;
        }
    }

    private static int insertionPoint(IntList il, int nRefHaps) {
        int index = il.binarySearch(nRefHaps);
        return index >= 0 ? index : -index - 1;
    }

    /**
     * Returns an array containing reference haplotype indices that
     * that are IBS with the specified target haplotype in an interval
     * beginning with the specified step. The returned array will contain fewer
     * than {@code this.nHapsPerStep()} haplotypes if the number of reference
     * haplotypes that are IBS with specified target haplotype in the specified
     * step is less than {@code this.nHapsPerStep()}.
     * @param hap a haplotype index
     * @param step a step index
     * @return an array containing reference haplotype indices that
     * that are IBS with the specified target haplotype
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.hapPairs().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public int[] ibsHaps(int hap, int step) {
        return ibsHaps[step][hap].clone();
    }

    /**
     * Return the data for genotype imputation in the marker window.
     * @return the data for genotype imputation in the marker window
     */
    public ImpData impData() {
        return impData;
    }

    /**
     * Returns the number of HMM states per marker.
     * @return the number of HMM states per marker
     */
    public int nStates() {
        return nStates;
    }

    /**
     * Returns the number of IBS steps in the marker window.
     * @return the number of IBS steps in the marker window
     */
    public int nSteps() {
        return ibsHaps.length;
    }

    /**
     * Returns the index of the first marker in the specified step.
     * @param step a step index
     * @return the index of the first marker in the specified step
     * @throws IllegalArgumentException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public int stepStart(int step) {
        return stepStarts[step];
    }
}
