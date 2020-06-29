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

import org.ngseq.sparkbeagle.blbutil.Utilities;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReferenceArray;

/**
 * <p>Class {@code ImpLS} computes HMM state probabilities
 * at genotyped markers in the target haplotypes.
 * </p>
 * <p>Instances of class {@code ImpLS} are thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ImpLS {

    private ImpLS() {
        // private constructor to prevent instantiation
    }

    /**
     * Returns estimated HMM state probabilities at genotyped markers
     * for each target haplotype.  The returned array maps each target
     * haplotype index to the target haplotype's HMM state probabilities.
     * @param impData the input data for genotype imputation
     * @return a list of HMM state probabilities at genotyped markers
     * for each target haplotype.
     * @throws NullPointerException if {@code impData == null}
     */
    public static AtomicReferenceArray<StateProbs> stateProbs(ImpData impData)  {
        int nThreads = impData.par().nthreads();
        int nTargHaps = impData.targGT().nHaps();
        final AtomicReferenceArray<StateProbs> stProbs
                = new AtomicReferenceArray<>(nTargHaps);
        final AtomicInteger hapIndices = new AtomicInteger(0);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        ImpIbs ibsHaps = new ImpIbs(impData);
        try {
            for (int j=0; j<nThreads; ++j) {
                ImpLSBaum baum = new ImpLSBaum(impData, ibsHaps);
                es.submit(() -> {
                    try {
                        int hap = hapIndices.getAndIncrement();
                        while (hap>=0 && hap<nTargHaps) {
                            StateProbs stateProbs = baum.impute(hap);
                            stProbs.set(hap, stateProbs);
                            hap = hapIndices.getAndIncrement();
                        }
                    }
                    catch(Throwable t) {
                        Utilities.exit(t);
                    }
                } );
            }
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable t) {
            Utilities.exit(t);
        }
        return stProbs;
    }
}
