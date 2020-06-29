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

import java.util.Arrays;

/**
 * <p>Class {@code ImpLSBaum} implements a Baum hidden Markov model
 * forward and backward algorithms for computing HMM state probabilities
 * at genotyped markers using IBS-matched reference haplotypes.
 * </p>
 * <p>Instances of class {@code ImpLSBaum} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ImpLSBaum {

    private final ImpData impData;
    private final ImpStates states;
    private final int nMarkers;
    private final int[][] hapIndices;
    private final boolean[][] allelesMatch;
    private final float[][] fwdVal;
    private final float[] bwdVal;
    private final StateProbsFactory alProbsFactory;
    private final int[] targAllele;

    /**
     * Creates a {@code LSHapBaum} instance from the specified data.
     *
     * @param impData the input data for genotype imputation
     * @param ibsHaps the IBS haplotype segments
     *
     * @throws NullPointerException if
     * {@code impData == null || ibsStates == null}
     */
    public ImpLSBaum(ImpData impData, ImpIbs ibsHaps) {
        this.impData = impData;
        this.states = new ImpStates(ibsHaps);
        this.nMarkers = impData.nClusters();
        int maxStates = states.nStates();
        this.hapIndices = new int[nMarkers][maxStates];
        this.allelesMatch = new boolean[nMarkers][maxStates];
        this.fwdVal = new float[nMarkers][maxStates];
        this.bwdVal = new float[maxStates];
        this.alProbsFactory = new StateProbsFactory(nMarkers);
        this.targAllele = new int[impData.nClusters()];
    }

    /**
     * <p>Returns HMM state probabilities at genotyped markers for the
     * specified  target haplotype. States with probabilities that are small
     * and inconsequential are excluded from the returned state probabilities.
     * </p>
     *
     * @param targHap a target haplotype index
     * @return HMM state probabilities at genotyped markers for the specified
     * target haplotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code targHap < 0 || hap >= this.impData().nTargHaps()}
     */
    public StateProbs impute(int targHap) {
        int lastMarker = impData.nClusters() - 1;
        int nStates = states.ibsStates(targHap, hapIndices, allelesMatch);
        setFwdValues(targHap, nStates);
        Arrays.fill(bwdVal, 0, nStates, 1.0f/nStates);
        float lastSum = 1.0f;
        for (int m=lastMarker; m>=0; --m) {
            lastSum = setBwdValue(m, nStates, lastSum);
        }
        return alProbsFactory.stateProbs(targHap, nStates, hapIndices, fwdVal);
    }

    /**
     * Return the input data for genotype imputation
     * @return the input data for genotype imputation
     */
    public ImpData impData() {
        return impData;
    }

    private void setFwdValues(int targHap, int nHaps) {
        int nRefHaps = impData.nRefHaps();
        float lastSum = 1.0f;
        for (int m=0; m<fwdVal.length; ++m) {
            float pRecomb = impData.pRecomb(m);
            float pErr = impData.errProb(m);
            float pNoErr = 1.0f - pErr;
            float shift = pRecomb/nHaps;
            float scale = (1.0f - pRecomb)/lastSum;
            float sum = 0.0f;
            targAllele[m] = impData.allele(m, nRefHaps + targHap);
            for (int j=0; j<nHaps; ++j) {
                float em = allelesMatch[m][j] ? pNoErr : pErr;
                fwdVal[m][j] = m==0 ? em : em*(scale*fwdVal[m-1][j] + shift);
                sum += fwdVal[m][j];
            }
            lastSum = sum;
        }
    }

    private float setBwdValue(int m, int nStates, float lastSum) {
        int mP1 = m + 1;
        float pRecomb = (mP1 < nMarkers) ? impData.pRecomb(mP1) : 0.0f;
        float pErr = impData.errProb(m);
        float pNoErr = 1.0f - pErr;
        float scale = (1.0f - pRecomb)/lastSum;
        float shift = pRecomb/nStates;
        float bwdValSum = 0f;
        float stateSum = 0f;
        for (int j=0; j<nStates; ++j) {
            bwdVal[j] = scale*bwdVal[j] + shift; // finish calculating bwd value
            fwdVal[m][j] *= bwdVal[j];  // store state probabilties in fwdVal[m]
            stateSum += fwdVal[m][j];

            float em = allelesMatch[m][j] ? pNoErr : pErr;
            bwdVal[j] *= em;
            bwdValSum += bwdVal[j];
        }
        for (int j=0; j<nStates; ++j) {
            fwdVal[m][j] /= stateSum;   // normalize state probabilities
        }
        return bwdValSum;
    }
}