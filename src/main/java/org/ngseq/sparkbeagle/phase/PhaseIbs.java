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

import org.ngseq.sparkbeagle.Par;
import org.ngseq.sparkbeagle.blbutil.Utilities;
import org.ngseq.sparkbeagle.ints.IndexArray;
import org.ngseq.sparkbeagle.ints.IntArray;
import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.ints.IntSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * <p>Class {@code PhaseIbs} identifies haplotypes that share a long
 * IBS segment with a specified haplotype.
 * </p>
 * <p>Instances of {@code PhaseIbs} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class PhaseIbs {

    private final PhaseData phaseData;
    private final int nTargHaps;
    private final long seed;
    private final int nStates;
    private final int nSteps;
    private final int nHapsPerStep;
    private final int ibsThreshold;

    private final int[] stepStarts;
    private final int[][][] ibsHaps; //[step][targ hap][ibs_set]

    /**
     * Constructs a new {@code PhaseIbs} object from the specified data.
     * @param phaseData the input data for an iteration of genotype phasing
     *
     * @throws IllegalArgumentException if
     * {@code hapPairs.markers().equals(map.markers()) == false}
     * @throws IllegalArgumentException if
     * {@code stepLength <= 0.0 || nStepsToMerge < 1 || nStates < 1}
     * @throws NullPointerException if {@code hapPairs == null || map == null}
     */
    public PhaseIbs(PhaseData phaseData) {
        Par par = phaseData.par();
        this.phaseData = phaseData;
        this.seed = phaseData.seed();
        this.nTargHaps = phaseData.nTargHaps();
        this.nStates = par.phase_states();
        this.nSteps = par.nsteps();
        int nStepsPerSegment = Math.round(par.phase_segment()/par.step());
        this.nHapsPerStep = (par.phase_states()/nStepsPerSegment)/2;  // for each of the 2 haps
        this.ibsThreshold = phaseData.burnin() ? 20*nHapsPerStep : nHapsPerStep + 2;

        this.stepStarts = stepStarts(phaseData);
        IndexArray[] codedSteps = IntStream.range(0, stepStarts.length)
                .parallel()
                .mapToObj(j -> codeStep(phaseData, stepStarts, j))
                .toArray(IndexArray[]::new);
        this.ibsHaps = IntStream.range(0, codedSteps.length)
                .parallel()
                .mapToObj(j -> getIbsHaps(codedSteps, j))
                .toArray(int[][][]::new);
    }

    private static int[] stepStarts(PhaseData phaseData) {
        double[] pos = phaseData.pos();
        double step = phaseData.par().step();
        IntList indices = new IntList(pos.length/10);
        indices.add(0);
        Random random = new Random(phaseData.seed());
        double nextPos =  pos[0] + random.nextDouble()*step;
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

    private static IndexArray codeStep(PhaseData phaseData, int[] starts,
            int startIndex) {
        int nTargHaps = phaseData.nTargHaps();
        int nHaps = phaseData.nHaps();
        int[] hapToSeq = IntStream.range(0, nHaps).map(i -> 1).toArray();
        int start = starts[startIndex];
        int end = (startIndex+1) < starts.length ?  starts[startIndex+1] : phaseData.nMarkers();

        int seqCnt = 2; // seq 0 is reserved for sequences not found in target
        for (int m=start; m<end; ++m) {
            int nAlleles = phaseData.marker(m).nAlleles();
            int[] seqMap = new int[seqCnt*nAlleles];

            seqCnt = 1;
            for (int h=0; h<nTargHaps; ++h) {
                int allele = phaseData.allele(m, h);
                int index = nAlleles*hapToSeq[h] + allele;
                if (seqMap[index] == 0) {
                    seqMap[index] = seqCnt++;
                }
                hapToSeq[h] = seqMap[index];
            }
            for (int h=nTargHaps; h<nHaps; ++h) {
                if (hapToSeq[h] != 0) {
                    int allele = phaseData.allele(m, h);
                    int index = hapToSeq[h]*nAlleles + allele;
                    hapToSeq[h] = seqMap[index];
                }
            }
        }
        IntArray intArray = IntArray.create(hapToSeq, seqCnt);
        return new IndexArray(intArray, seqCnt);
    }

    private int[][] getIbsHaps(IndexArray[] codedSteps, int step) {
        Random rand = new Random(seed*(step+1)); // thread-confined
        int[][] results = new int[nTargHaps][];
        int nStepsToMerge = Math.min(nSteps, codedSteps.length - step);
        List<IntList> children = initPartition(codedSteps[step]);
        List<IntList> nextParents = new ArrayList<>(children.size());
        initUpdateResults(children, nextParents, results);
        for (int i=1; i<nStepsToMerge; ++i) {
            int initCapacity = Math.min(nTargHaps, 2*nextParents.size());
            List<IntList> parents = nextParents;
            nextParents = new ArrayList<>(initCapacity);
            IndexArray codedStep = codedSteps[step+i];
            for (int j=0, nj=parents.size(); j<nj; ++j) {
                IntList parent = parents.get(j);
                children = partition(parent, codedStep);
                updateResults(parent, children, nextParents, results, rand);
            }
        }
        finalUpdateResults(nextParents, results, rand);
        return results;
    }

    private List<IntList> initPartition(IndexArray codedStep) {
        IntList[] list = new IntList[codedStep.valueSize()];
        IntArray hap2Seq = codedStep.intArray();
        int nHaps = hap2Seq.size();
        List<IntList> children = new ArrayList<>();
        for (int h=0; h<nTargHaps; ++h) {
            int seq = hap2Seq.get(h);
            if (list[seq]==null) {
                list[seq] = new IntList();
                children.add(list[seq]);
            }
            list[seq].add(h);
        }
        for (int h=nTargHaps; h<nHaps; ++h) {
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
        int nTarg = insertionPoint(parent, nTargHaps);
        for (int k=0; k<nTarg; ++k) {
            int hap = parent.get(k);
            int seq = hap2Seq.get(hap);
            if (list[seq]==null) {
                list[seq] = new IntList();
                children.add(list[seq]);
            }
            list[seq].add(hap);
        }
        for (int k=nTarg; k<nHaps; ++k) {
            int hap = parent.get(k);
            int seq = hap2Seq.get(hap);
            if (list[seq]!=null) {
                list[seq].add(hap);
            }
        }
        return children;
    }

    private void initUpdateResults(List<IntList> children,
            List<IntList> nextParents, int[][] result) {
        for (int j=0, n=children.size(); j<n; ++j) {
            IntList child = children.get(j);
            if (child.size() < ibsThreshold) {
                setResult(child.toArray(), child, result);
            }
            else {
                nextParents.add(child);
            }
        }
    }

    private void updateResults(IntList parent, List<IntList> children,
            List<IntList> nextParents, int[][] results, Random rand) {
        int[] ia = null;
        for (int j=0, n=children.size(); j<n; ++j) {
            IntList child = children.get(j);
            if (child.size() < ibsThreshold) {
                if (ia==null) {
                    ia = randomizedArray(parent, rand);
                }
                setResult(ia, child, results);
            }
            else {
                nextParents.add(child);
            }
        }
    }

    private void finalUpdateResults(List<IntList> children, int[][] results,
            Random rand) {
        for (int j=0, n=children.size(); j<n; ++j) {
            IntList child = children.get(j);
            int[] ia = randomizedArray(child, rand);
            setResult(ia, child, results);
        }
    }

    private static int[] randomizedArray(IntList list, Random rand) {
        int[] ia = list.toArray();
        Utilities.shuffle(ia, rand);
        return ia;
    }

    private void setResult(int[] parent, IntList child, int[][] result) {
        int nTarg = insertionPoint(child, nTargHaps);
        for (int j=0; j<nTarg; ++j) {
            result[child.get(j)] = parent;
        }
    }

    private static int insertionPoint(IntList il, int nRefHaps) {
        int index = il.binarySearch(nRefHaps);
        return index >= 0 ? index : -index - 1;
    }

    /**
     * Returns a set containing the specified number of haplotype
     * indices that are IBS with haplotype {@code h1} in the specified
     * step. The returned set is guaranteed to not contain indices
     * {@code h1} or {@code h2}.  The returned set will contain fewer than
     * {@code this.nHapsPerStep()} haplotypes if the number of haplotypes
     * that are IBS with haplotype {@code h1} is less than
     * {@code this.nHapsPerStep()}.
     * @param h1 a haplotype index
     * @param h2 a haplotype index
     * @param step a step index
     * @return a set containing IBS haplotype indices
     * @throws IndexOutOfBoundsException if
     * {@code h1 < 0 || h1 >= this.hapPairs().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code h2 < 0 || h2 >= this.hapPairs().nHaps()}
      * @throws IndexOutOfBoundsException if
     * {@code step < 0 || step >= this.nSteps()}
     * @throws IllegalArgumentException if
     * {@code nHaps < 0 || nHaps > (1 << 30)}
     */
    public IntSet ibsHaps(int h1, int h2, int step) {
        IntSet intSet = new IntSet(nHapsPerStep);
        Random rand = new Random(seed*(step+1) + h1);
        int[] ia = ibsHaps[step][h1];
        int index = rand.nextInt(ia.length);
        for (int j=0; j<ia.length; ++j) {
            if (index==ia.length) {
                index -= ia.length;
            }
            int h = ia[index];
            if (h!=h1 && h!=h2 && intSet.add(h)) {
                if (intSet.size()==nHapsPerStep) {
                    return intSet;
                }
            }
            ++index;
        }
        return intSet;
    }

    /**
     * Return the data used to phase genotypes in a marker window.
     * @return the data used to phase genotypes in a marker window
     */
    public PhaseData phaseData() {
        return phaseData;
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
     * Returns the first marker index in the specified step.
     * @param step a step index
     * @return the first marker index in the specified step
     * @throws IllegalArgumentException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public int stepStart(int step) {
        return stepStarts[step];
    }
}
