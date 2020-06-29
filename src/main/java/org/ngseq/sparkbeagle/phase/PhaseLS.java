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

import org.ngseq.sparkbeagle.blbutil.Utilities;
import org.ngseq.sparkbeagle.math.Regress;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * <p>Class {@code PhaseLS} estimated genotypes phase using
 * a haploid Li and Stephens hidden Markov model.  It
 * uses a rolling window of reference haplotypes for phasing each sample.
 * </p>
 * <p>Instances of class {@code PhaseLS} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PhaseLS {

    private PhaseLS() {
        // private constructor to prevent instantiation
    }

    /**
     * Estimates and stores phased haplotypes for the target samples.
     * The returned haplotypes are ordered by increasing sample index.
     * @param phaseData the input data for an iteration of genotype phasing
     * @param estPhase the estimated sample phase
     * @param recombRegress used for regression of state-switch probability
     * on inter-marker genetic distance
     * @throws IllegalArgumentException if {@code nItsRemaining < 1}
     * @throws NullPointerException if any parameter is {@code null}
     * or if any element of {@code phase} is {@code null}
     */
    public static void run(PhaseData phaseData, EstPhase estPhase,
            Regress recombRegress) {
        int nThreads = phaseData.par().nthreads();
        int nSamples = phaseData.nTargSamples();

        PhaseIbs phaseIbs = new PhaseIbs(phaseData);

        AtomicInteger sampleIndices = new AtomicInteger(0);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        try {
            for (int j=0; j<nThreads; ++j) {
                PhaseBaum1 baum = new PhaseBaum1(phaseData, phaseIbs);
                es.submit(() -> {
                    try {
                        int sample = sampleIndices.getAndIncrement();
                        while (sample>=0 && sample<nSamples) {
                            baum.phase(estPhase, sample, recombRegress);
                            sample = sampleIndices.getAndIncrement();
                        }
                    }
                    catch (Throwable t) {
                        Utilities.exit(t);
                    }
                } );
            }
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable e) {
            Utilities.exit("ERROR", e);
        }
    }
}
