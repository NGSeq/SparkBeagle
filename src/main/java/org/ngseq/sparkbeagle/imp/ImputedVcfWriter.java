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

import org.ngseq.sparkbeagle.CurrentData;
import org.ngseq.sparkbeagle.blbutil.FloatList;
import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.vcf.Markers;
import org.ngseq.sparkbeagle.vcf.RefGT;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.stream.IntStream;

/**
 * <p>Class {@code ImputedVcfWriter} writes observed and imputed genotypes
 * to a VCF output file.
 * </p>
 * <p>Instances of class {@code ImputedVcfWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ImputedVcfWriter {

    private final ImpData impData;
    private final int targCluster;
    private final int refStart;
    private final int clustEnd;
    private final int refEnd;

    private final IntList indices = new IntList(4);
    private final IntList hashes = new IntList(4);
    private final FloatList seqProbs = new FloatList(4);
    private final FloatList seqProbsP1 = new FloatList(4);

    /**
     * Constructs a new {@code ImputedVcfWriter} instance from
     * the specified data.
     * @param impData the input data for genotype imputation
     * @param targCluster the index of the target marker cluster in the
     * interval of reference markers that will be printed
     * @param refStart an lower bound (inclusive) on the reference markers
     * that will be printed
     * @param refEnd an upper bound (exclusive) on the reference markers
     * that will be printed
     * @throws IndexOutOfBoundsException if
     * {@code (refStart < 0 || refEnd > impData.refGT().nMarkers())}
     * @throws IndexOutOfBoundsException if
     * {@code targCluster < 0 || targCluster >= impData.nClusters()}
     * @throws NullPointerException if {@code impData == null}
     */
    public ImputedVcfWriter(ImpData impData, int targCluster, int refStart,
            int refEnd) {
        if (refStart < 0) {
            throw new IndexOutOfBoundsException(String.valueOf(refStart));
        }
        if (refEnd > impData.refGT().nMarkers()) {
            throw new IndexOutOfBoundsException(String.valueOf(refEnd));
        }
        this.impData = impData;
        this.targCluster = targCluster;
        if (targCluster==0) {
            this.refStart = refStart;
        }
        else {
            int rf = impData.refClusterStart(targCluster);
            this.refStart = Math.max(refStart, rf);
        }
        if (targCluster < impData.nClusters()-1) {
            int tmpClustEnd = Math.max(refStart, impData.refClusterEnd(targCluster));
            this.clustEnd = Math.min(tmpClustEnd, refEnd);
            this.refEnd = Math.min(impData.refClusterStart(targCluster + 1), refEnd);
        }
        else {
            this.clustEnd = refEnd;
            this.refEnd = refEnd;
        }
    }

    /**
     * Writes the VCF records to the specified {@code PrintWriter}.
     * @param stateProbs the imputed state probabilities at genotyped
     * markers in the target samples
     * @param out the {@code PrintWriter} to which VCF records will be written
     * @throws NullPointerException if {@code stateProbs == null || out == null}
     * @throws NullPointerException if there exists a {@code j} satisfying
     * {@code (0 <= j && j < stateProbs.size()) && (stateProbs.get(j) == null)}
     */
    public void appendRecords(AtomicReferenceArray<StateProbs> stateProbs,
            PrintWriter out) {
        if (refStart >= refEnd) {
            return;
        }
        RefHapHash refHapHash = new RefHapHash(stateProbs, targCluster,
                impData.refGT(), refStart, refEnd);

        ImputedRecBuilder[] recBuilders = recBuilders();
        float[][] a1Probs = alProbs();
        float[][] a2Probs = alProbs();
        boolean[] isImputed = isImputed();
        for (int h=0, n=stateProbs.length(); h<n; h+=2) {
            setAlProbs(stateProbs.get(h), refHapHash, a1Probs);
            setAlProbs(stateProbs.get(h+1), refHapHash, a2Probs);

            for (int m=0; m<a1Probs.length; ++m) {
                if (isImputed[m]==false) {
                    setToObsAlleles(a1Probs, a2Probs, m, h);
                }
                recBuilders[m].addSampleData(a1Probs[m], a2Probs[m]);
                Arrays.fill(a1Probs[m], 0f);
                Arrays.fill(a2Probs[m], 0f);
            }
        }
        for (int m=0; m<a1Probs.length; ++m) {
            recBuilders[m].printRec(out, isImputed[m]);
        }
    }

    private void setAlProbs(StateProbs stateProbs, RefHapHash rhh,
            float[][] alProbs) {
        clearLists(indices, hashes, seqProbs, seqProbsP1);
        int[] alleles = new int[rhh.end() - rhh.start()];
        for (int j=0, n=stateProbs.nStates(targCluster); j<n; ++j) {
            int hap = stateProbs.refHap(targCluster, j);
            float val = stateProbs.probs(targCluster, j);
            float valP1 = stateProbs.probsP1(targCluster, j);
            int index = rhh.hap2Index(hap);
            int hash = rhh.hash(index);
            int i = 0;
            while (i<hashes.size() && hashes.get(i)!=hash) {
                ++i;
            }
            if (i == hashes.size()) {
                indices.add(index);
                hashes.add(hash);
                seqProbs.add(val);
                seqProbsP1.add(valP1);
            }
            else {
                seqProbs.addToElement(i, val);
                seqProbsP1.addToElement(i, valP1);
            }
        }
        setAlProbs(alProbs, rhh, alleles);
    }

    private void setAlProbs(float[][] alProbs, RefHapHash rhh, int[] alleles) {
        int nSeq = seqProbs.size();
        if (nSeq==1) {
            int index = indices.get(0);
            rhh.setAlleles(index, alleles);
            for (int m=refStart; m<refEnd; ++m) {
                int mm = m - refStart;
                alProbs[mm][alleles[mm]] = 1.0f;
            }
        }
        else {
            for (int j=0; j<nSeq; ++j) {
                int index = indices.get(j);
                rhh.setAlleles(index, alleles);
                float prob = seqProbs.get(j);
                float probP1 = seqProbsP1.get(j);
                for (int m=refStart; m<clustEnd; ++m) {
                    int mm = m - refStart;
                    alProbs[mm][alleles[mm]] += prob;
                }
                for (int m=clustEnd; m<refEnd; ++m) {
                    double wt = impData.weight(m);
                    int mm = m - refStart;
                    alProbs[mm][alleles[mm]] += (wt*prob + (1-wt)*probP1);
                }
            }
        }
    }

    private void setToObsAlleles(float[][] a1Probs, float[][] a2Probs, int m,
            int targHap) {
        Arrays.fill(a1Probs[m], 0f);
        Arrays.fill(a2Probs[m], 0f);
        int preClustIndex = impData.cd().targMarkerIndex(refStart + m);
        int a1 = impData.targGT().allele(preClustIndex, targHap);
        int a2 = impData.targGT().allele(preClustIndex, targHap+1);
        a1Probs[m][a1] = 1f;
        a2Probs[m][a2] = 1f;
    }

    private ImputedRecBuilder[] recBuilders() {
        RefGT refGT = impData.refGT();
        boolean gp = impData.par().gp();
        boolean ap = impData.par().ap();
        int nTargSamples = impData.nTargSamples();
        return IntStream.range(refStart, refEnd)
                .mapToObj(m -> new ImputedRecBuilder(refGT.marker(m),
                        nTargSamples, ap, gp))
                .toArray(ImputedRecBuilder[]::new);
    }

    private float[][] alProbs() {
        Markers refMarkers = impData.refGT().markers();
        return IntStream.range(refStart, refEnd)
                .mapToObj(m -> new float[refMarkers.marker(m).nAlleles()])
                .toArray(float[][]::new);
    }

    private boolean[] isImputed() {
        CurrentData cd = impData.cd();
        boolean[] isImputed = new boolean[refEnd - refStart];
        for (int j=0; j<isImputed.length; ++j) {
            if (cd.targMarkerIndex(refStart + j) == -1) {
                isImputed[j] = true;
            }
        }
        return isImputed;
    }

    private static void clearLists(IntList hapIndices, IntList seqIndices,
            FloatList seqProbs, FloatList seqProbsP1) {
        hapIndices.clear();
        seqIndices.clear();
        seqProbs.clear();
        seqProbsP1.clear();
    }
}
