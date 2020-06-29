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

import org.ngseq.sparkbeagle.blbutil.FloatArray;
import org.ngseq.sparkbeagle.blbutil.FloatList;
import org.ngseq.sparkbeagle.ints.IntArray;
import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.math.Regress;

import java.util.Arrays;

/**
 * <p>Class {@code PhaseBaum1} implements the forward and backward algorithms
 * for a haploid Li and Stephens hidden Markov model.  It evaluates
 * probabilities of diplotypes in a window defined by two heterozygote
 * genotypes.
 * </p>
 * <p>Instances of class {@code PhaseBaum1} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PhaseBaum1 {

    private final PhaseData phaseData;
    private final FloatArray genDist;
    private final FloatArray pRecomb;
    private final int[][] stateAlleles ;
    private final PhaseStates states;
    private final FloatList lrList;
    private final int nMarkers;
    private final float pErr;
    private final float pNoErr;

    private final int[] hap1;
    private final int[] hap2;

    private final float[] fwd;
    private final float[] fwd1;
    private final float[] fwd2;

    private final float[] bwd;
    private final float[] bwd1;
    private final float[] bwd2;
    private final float[][] savedBwd1;
    private final float[][] savedBwd2;

    private int nStates;
    private float shift;
    private float scale;
    private float scale1;
    private float scale2;

    private float sum;
    private float sum1;
    private float sum2;

    private final boolean[] switchHaps;

    /**
     * Creates a {@code PhaseLSBaum} instance from the specified data.
     * The contract for this class is unspecified if any element of the
     * {@code unphased} {@code AtomicReferenceArray} is {@code null}.
     *
     * @param phaseData the input data for an iteration of genotype phasing
     * @param phaseIbs the IBS haplotype segments

     * @throws IllegalArgumentException if {@code nItsRemaining < 1}
     * @throws NullPointerException if any input parameter is {@code null}
     */
    public PhaseBaum1(PhaseData phaseData, PhaseIbs phaseIbs) {
        this.phaseData = phaseData;
        this.genDist = phaseData.genDist();
        this.pRecomb = phaseData.pRecomb();
        this.lrList = new FloatList(200);
        this.nMarkers = phaseData.nMarkers();
        this.states = new PhaseStates(phaseIbs);
        this.pErr = phaseData.par().err();
        this.pNoErr = 1 - pErr;

        int maxStates = states.nStates();
        this.stateAlleles = new int[phaseData.nMarkers()][maxStates];

        this.hap1 = new int[nMarkers];
        this.hap2 = new int[nMarkers];

        this.bwd = new float[maxStates];
        this.bwd1 = new float[maxStates];
        this.bwd2 = new float[maxStates];
        this.savedBwd1 = new float[nMarkers][maxStates];
        this.savedBwd2 = new float[nMarkers][maxStates];

        this.fwd = new float[maxStates];
        this.fwd1 = new float[maxStates];
        this.fwd2 = new float[maxStates];

        this.switchHaps = new boolean[nMarkers];
    }

    /**
     * Estimates and returns phased haploytpes for the specified sample
     * @param estPhase the estimated sample phase
     * @param sample a sample index
     * @param recombRegress object for storing data points for regression of
     * state-switch probability on inter-marker genetic distance or {@code null}
     * if no regression is to be performed
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= samplePhase.gl().nSamples()}
     */
    public void phase(EstPhase estPhase, int sample, Regress recombRegress) {
        boolean missingGTs = estPhase.hasMissingGT(sample);
        boolean hasUnphasedHet = estPhase.hasUnphasedHets(sample);
        if (missingGTs || hasUnphasedHet) {
            this.nStates = states.ibsStates(sample, stateAlleles);
            estPhase.setHapsWithMaskedMissingAlleles(sample, hap1, hap2);

            phaseHets(estPhase, sample);
            if (missingGTs || recombRegress!=null) {
                if (recombRegress==null) {
                    imputeMissingAlleles(hap1);
                    imputeMissingAlleles(hap2);
                }
                else {
                    imputeMissingAlleles(hap1, recombRegress);
                    imputeMissingAlleles(hap2, recombRegress);
                }
            }
            estPhase.setHapPair(sample, hap1, hap2);
        }
    }

    private void phaseHets(EstPhase estPhase, int sample) {
        lrList.clear();
        Arrays.fill(switchHaps, false); // can optimize this later
        IntArray unph = estPhase.getUnphasedHets(sample);
        if (unph.size()>0) {
            setBwd(unph);
            Arrays.fill(fwd, 0, nStates, 1.0f/nStates);
            sum = 1.0f;
            for (int j=0, n=unph.size(); j<n; ++j) {
                int het1 = j==0 ? 0 : unph.get(j-1);
                int het2 = unph.get(j);
                setFwd(het1, het2);
                phaseHet(het2);
            }
            updatePhase();
            if (phaseData.burnin()==false) {
                updateUnphased(estPhase, sample);
            }
        }
    }

    private void updatePhase() {
        boolean switchAlleles = false;
        for (int m=0; m<nMarkers; ++m) {
            if (switchHaps[m]) {
                switchAlleles = !switchAlleles;
            }
            if (switchAlleles) {
                int tmp = hap1[m];
                hap1[m] = hap2[m];
                hap2[m] = tmp;
            }
        }
    }

    private void setBwd(IntArray unph) {
        Arrays.fill(bwd, 0, nStates, 1.0f/nStates);
        for (int j=unph.size()-1; j>=0; --j) {
            final int het2 = unph.get(j);
            final int het3 = (j+1)==unph.size() ? nMarkers : unph.get(j+1);
            System.arraycopy(bwd, 0, bwd1, 0, nStates);
            System.arraycopy(bwd, 0, bwd2, 0, nStates);
            for (int m=(het3-2); m>(het2-2); --m) {
                bwdUpdate(m);
            }
            System.arraycopy(bwd1, 0, savedBwd1[het2-1], 0, nStates);
            System.arraycopy(bwd2, 0, savedBwd2[het2-1], 0, nStates);
        }
    }

    private void bwdUpdate(int m) {
        int mP1 = m + 1;
        sum1 = sum2 = sum = 0.0f;
        for (int k=0; k<nStates; ++k) {
            float em1 = hap1[mP1]==-1 ? 1.0f : (stateAlleles[mP1][k]==hap1[mP1] ? pNoErr : pErr);
            float em2 = hap2[mP1]==-1 ? 1.0f : (stateAlleles[mP1][k]==hap2[mP1] ? pNoErr : pErr);
            float em = em1==em2 ? em1 : 1.0f;
            bwd[k]  *= em;
            bwd1[k] *= em1;
            bwd2[k] *= em2;
            sum += bwd[k];
            sum1 += bwd1[k];
            sum2 += bwd2[k];
        }
        setScaleAndShift1(mP1);
        for (int k=0; k<nStates; ++k) {
            bwd[k] = scale*bwd[k] + shift;
            bwd1[k] = scale1*bwd1[k] + shift;
            bwd2[k] = scale2*bwd2[k] + shift;
        }
    }

    private void setFwd(int het1, int het2) {
        System.arraycopy(fwd, 0, fwd1, 0, nStates);
        System.arraycopy(fwd, 0, fwd2, 0, nStates);
        sum1 = sum2 = sum;
        for (int m=het1; m<het2; ++m)  {
            fwdUpdate(m);
        }
    }

    private void fwdUpdate(int m) {
        setScaleAndShift1(m);
        sum = sum1 = sum2 = 0f;
        for (int k=0; k<nStates; ++k) {
            float em1 = hap1[m]==-1 ? 1.0f : (stateAlleles[m][k]==hap1[m] ? pNoErr : pErr);
            float em2 = hap2[m]==-1 ? 1.0f : (stateAlleles[m][k]==hap2[m] ? pNoErr : pErr);
            float em = em1==em2 ? em1 : 1.0f;
            fwd[k] = em*(scale*fwd[k] + shift);
            fwd1[k] = em1*(scale1*fwd1[k] + shift);
            fwd2[k] = em2*(scale2*fwd2[k] + shift);
            sum += fwd[k];
            sum1 += fwd1[k];
            sum2 += fwd2[k];
        }
    }

    private void imputeMissingAlleles(int[] obs, Regress recombRegress) {
        Arrays.fill(bwd, 1.0f);
        System.arraycopy(bwd, 0, savedBwd1[nMarkers-1], 0, nStates);
        for (int m=nMarkers-2; m>=0; --m) {
            int mP1 = m + 1;
            sum = 0.0f;
            for (int k=0; k<nStates; ++k) {
                if (obs[mP1]>=0) {
                    bwd[k] *= (stateAlleles[mP1][k]==obs[mP1] ? pNoErr : pErr);
                }
                sum += bwd[k];
            }
            setScaleAndShift(mP1);
            for (int k=0; k<nStates; ++k) {
                bwd[k] = scale*bwd[k] + shift;
            }
            System.arraycopy(bwd, 0, savedBwd1[m], 0, nStates);
        }

        Arrays.fill(fwd, 0, nStates, 1.0f/nStates);
        sum = 1.0f;
        updateFwd(obs, 0);
        if (obs[0]<0) {
            obs[0] = maxIndex(alProbs(0));
        }
        for (int m=1; m<nMarkers; ++m) {
            updateFwdAndRecomb(obs, m, recombRegress);
            if (obs[m]<0) {
                obs[m] = maxIndex(alProbs(m));
            }
        }
    }

    private void imputeMissingAlleles(int[] obs) {
        Arrays.fill(bwd, 0, nStates, 1.0f);
        if (obs[nMarkers-1]<0) {
            Arrays.fill(savedBwd1[nMarkers-1], 0, nStates, 1.0f);
        }
        for (int m=nMarkers-2; m>=0; --m) {
            int mP1 = m + 1;
            sum = 0.0f;
            for (int k=0; k<nStates; ++k) {
                if (obs[mP1]>=0) {
                    bwd[k] *= (stateAlleles[mP1][k]==obs[mP1] ? pNoErr : pErr);
                }
                sum += bwd[k];
            }
            setScaleAndShift(mP1);
            for (int k=0; k<nStates; ++k) {
                bwd[k] = scale*bwd[k] + shift;
            }
            if (obs[m]<0) {
                System.arraycopy(bwd, 0, savedBwd1[m], 0, nStates);
            }
        }

        Arrays.fill(fwd, 0, nStates, 1.0f/nStates);
        sum = 1.0f;
        for (int m=0; m<nMarkers; ++m) {
            updateFwd(obs, m);
            if (obs[m]<0) {
                obs[m] = maxIndex(alProbs(m));
            }
        }
    }

    private void updateFwdAndRecomb(int[] obs, int m, Regress recombRegress) {
        setScaleAndShift(m);
        float factor = nStates/(nStates-1.0f);
        float numerator = 0.0f;
        float denominator = 0.0f;
        float lastSum = sum;
        float termDenominator = lastSum/shift;
        sum = 0.0f;
        for (int k=0; k<nStates; ++k) {
            float term = (lastSum - fwd[k])/termDenominator;  // must calculate before fwd[k] update
            fwd[k] = scale*fwd[k] + shift;
            if (obs[m]>=0) {
                float em = stateAlleles[m][k]==obs[m] ? pNoErr : pErr;
                term *=  em;
                fwd[k] *= em;
            }
            numerator += term*savedBwd1[m][k];
            denominator += fwd[k]*savedBwd1[m][k];
            sum += fwd[k];
        }
        recombRegress.add(genDist.get(m), factor*(numerator/denominator));
    }

    private void updateFwd(int[] obs, int m) {
        setScaleAndShift(m);
        sum = 0.0f;
        for (int k=0; k<nStates; ++k) {
            fwd[k] = scale*fwd[k] + shift;
            if (obs[m]>=0) {
                fwd[k] *= (stateAlleles[m][k]==obs[m] ? pNoErr : pErr);
            }
            sum += fwd[k];
        }
    }

    private float[] alProbs(int m) {
        float[] alleleFreq = new float[phaseData.marker(m).nAlleles()];
        for (int k=0; k<nStates; ++k) {
            alleleFreq[stateAlleles[m][k]] += fwd[k]*savedBwd1[m][k];
        }
        return alleleFreq;
    }

    private void updateUnphased(EstPhase estPhase, int sample) {
        IntArray prevUnph = estPhase.getUnphasedHets(sample);
        IntList nextUnph = new IntList();
        float threshold = threshold(lrList, phaseData.nItsRemaining());
        for (int j=0, n=prevUnph.size(); j<n; ++j) {
            if (lrList.get(j) < threshold) {
                nextUnph.add(prevUnph.get(j));
            }
        }
        estPhase.setUnphasedHets(sample, IntArray.create(nextUnph, nMarkers));
    }

    private static float threshold(FloatList lrList, int itsRemaining) {
        float[] lra = lrList.toArray();
        Arrays.sort(lra);
        if (itsRemaining==1) {
            return lra[0];
        }
        else {
            int nUnphasedHets = lra.length + 1;
            double prop = Math.pow(1.0/nUnphasedHets, 1.0/itsRemaining);
            int rank = (int) Math.floor(prop*lra.length + 0.5);
            return lra[rank<lra.length ? rank : (lra.length-1)];
        }
    }

    private void phaseHet(int m) {
        int mM1 = m - 1;
        float[] b1 = savedBwd1[mM1];
        float[] b2 = savedBwd2[mM1];
        float p11 = 0.0f;
        float p12 = 0.0f;
        float p21 = 0.0f;
        float p22 = 0.0f;
        for (int k=0; k<nStates; ++k) {
            p11 += fwd1[k]*b1[k];
            p12 += fwd1[k]*b2[k];
            p21 += fwd2[k]*b1[k];
            p22 += fwd2[k]*b2[k];
        }
        float lr = (p11*p22)/(p12*p21);
        if (lr>=1.0f) {
            lrList.add(lr);
        }
        else {
            lrList.add(1.0f/lr);
            switchHaps[m] = true;
        }
    }

    private void setScaleAndShift(int m) {
        float switchProb = pRecomb.get(m);
        shift = (switchProb/nStates);
        scale = (1.0f - switchProb)/sum;
    }

    private void setScaleAndShift1(int m) {
        float switchProb = pRecomb.get(m);
        shift = (switchProb/nStates);
        scale = (1.0f - switchProb)/sum;
        scale1 = (1.0f - switchProb)/sum1;
        scale2 = (1.0f - switchProb)/sum2;
    }

    private static int maxIndex(float[] da) {
        int maxIndex = 0;
        for (int j=1; j<da.length; ++j) {
            if (da[j] > da[maxIndex]) {
                maxIndex = j;
            }
        }
        return maxIndex;
    }

    private void switchHaps(int start, int end) {
        for (int m=start; m<end; ++m) {
            int tmp = hap1[m];
            hap1[m] = hap2[m];
            hap2[m] = tmp;
        }
    }
}
