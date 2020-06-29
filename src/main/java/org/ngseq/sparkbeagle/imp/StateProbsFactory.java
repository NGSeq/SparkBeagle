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

import org.ngseq.sparkbeagle.blbutil.FloatList;
import org.ngseq.sparkbeagle.ints.IntList;

/**
 * <p>Class {@code StateProbsFactory} stores HMM state probabilities
 * that that can be used to imputed impute missing HMM state probabilities
 * using linear interpolation.
 * </p>
 * <p>Instances of class {@code StateProbsFactory} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class StateProbsFactory {

    private final int nTargMarkers;
    private final int nTargMarkersM1;

    private final IntList hapList;
    private final FloatList probList;
    private final FloatList probP1List;

    private final int[][] haps;
    private final float[][] probs;
    private final float[][] probsP1;

    /**
     * Creates an {@code StateProbsFactory} instance from the specified data.
     *
     * @param nTargMarkers the number of target markers
     * @throws IllegalArgumentException if {@code nTargMarkers <= 0}
     */
    public StateProbsFactory(int nTargMarkers) {
        if (nTargMarkers <= 0) {
            throw new IllegalArgumentException(String.valueOf(nTargMarkers));
        }
        this.nTargMarkers = nTargMarkers;
        this.nTargMarkersM1 = nTargMarkers - 1;

        this.hapList = new IntList(50);
        this.probList = new FloatList(50);
        this.probP1List = new FloatList(50);

        this.haps = new int[nTargMarkers][];
        this.probs = new float[nTargMarkers][];
        this.probsP1 = new float[nTargMarkers][];
    }

    /**
     * Returns the number of target marker.
     * @return the number of target marker
     */
    public int nTargMarkers() {
        return nTargMarkers;
    }

    /**
     * Stores the reference haplotypes and state probabilities that
     * will be used to imputed missing state probabilities using linear
     * interpolation.  If a HMM state probability at a marker or the
     * following marker are greater than or equal to
     * {@code Math.min(0.005f, 1.0f/nStates)}, then the reference haplotype,
     * the state probability at the marker, and the state probability at
     * the succeeding marker are stored.
     * @param targHap a target data haplotype index
     * @param nStates the number of HMM states at each marker
     * @param hapIndices the reference haplotype indices corresponding to
     * each state at each marker
     * @param stateProbs the HMM state probabilities at each marker
     * @throws IllegalArgumentException if
     * {@code hapIndices.length != this.nTargMarkers()}
     * @return the state probabilities
     * @throws IllegalArgumentException if
     * {@code stateProbs.length != this.nTargMarkers()}
     * @throws IndexOutOfBoundsException if there exists {@code j} satisfying
     * {@code (0 <= j && j <= this.nTargMarkers())} such that
     * {@code (hapIndices[j].length < nStates || stateProbs[j].length < nStates)}
     * @throws NullPointerException if
     * {@code hapIndices == null || stateProbs == null}
     * @throws NullPointerException if there exists {@code j} satisfying
     * {@code (0 <= j && j <= this.nTargMarkers())} such that
     * {@code (hapIndices[j] == null || stateProbs[j] == null)}
     */
    public StateProbs stateProbs(int targHap, int nStates, int[][] hapIndices,
            float[][] stateProbs) {
        if (hapIndices.length != nTargMarkers) {
            throw new IllegalArgumentException(String.valueOf(hapIndices.length));
        }
        if (stateProbs.length != nTargMarkers) {
            throw new IllegalArgumentException(String.valueOf(stateProbs.length));
        }
        float threshold = threshold(nStates);
        for (int m=0; m<nTargMarkers; ++m) {
            hapList.clear();
            probList.clear();
            probP1List.clear();
            int mP1 = m<nTargMarkersM1 ? m+1 : m;
            for (int j=0; j<nStates; ++j) {
                if (stateProbs[m][j] > threshold || stateProbs[mP1][j] > threshold) {
                    hapList.add(hapIndices[m][j]);
                    probList.add(stateProbs[m][j]);
                    probP1List.add(stateProbs[mP1][j]);
                }
            }
            haps[m] = hapList.toArray();
            probs[m] = probList.toArray();
            probsP1[m] = probP1List.toArray();
        }
        return new BasicStateProbs(targHap, haps, probs, probsP1);
    }

    private float threshold(int nStates) {
        return Math.min(0.005f, 0.9999f/nStates);
    }

    private static class BasicStateProbs implements StateProbs {
        private final int targHap;
        private final int[][] haps;
        private final float[][] probs;
        private final float[][] probsP1;

        public BasicStateProbs(int targHap, int[][] haps, float[][] probs,
                float[][] probsP1) {
            this.targHap = targHap;
            this.haps = haps.clone();
            this.probs = probs.clone();
            this.probsP1 = probsP1.clone();
        }

        @Override
        public int nTargMarkers() {
            return haps.length;
        }

        @Override
        public int nStates(int marker) {
            return haps[marker].length;
        }

        @Override
        public int targHap() {
            return targHap;
        }

        @Override
        public int refHap(int marker, int index) {
            return haps[marker][index];
        }

        @Override
        public float probs(int marker, int index) {
            return probs[marker][index];
        }

        @Override
        public float probsP1(int marker, int index) {
            return probsP1[marker][index];
        }
    }
}
