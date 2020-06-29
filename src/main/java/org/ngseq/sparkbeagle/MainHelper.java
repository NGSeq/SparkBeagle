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
package org.ngseq.sparkbeagle;

import org.ngseq.sparkbeagle.ints.LongArray;
import org.ngseq.sparkbeagle.math.Regress;
import org.ngseq.sparkbeagle.phase.EstPhase;
import org.ngseq.sparkbeagle.phase.PhaseData;
import org.ngseq.sparkbeagle.phase.PhaseLS;
import org.ngseq.sparkbeagle.sample.InitTargHapPairs;
import org.ngseq.sparkbeagle.vcf.GT;
import org.ngseq.sparkbeagle.vcf.GeneticMap;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.atomic.AtomicReferenceArray;

/**
 * Class {@code MainHelper} is an auxiliary class with methods called by
 * the {@code main.Main} class.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MainHelper {

    private final Par par;
    private final RunStats runStats;

    /**
     * Constructs a new {@code MainHelper} instance.
     * @param par the command line parameters
     * @param genMap the genetic map
     * @param runStats the class for collecting and printing run-time statistics
     * @throws NullPointerException
     * if {@code (par == null || genMap == null || runStarts == null)}
     */
    MainHelper(Par par, GeneticMap genMap, RunStats runStats) {
        this.par = par;
        this.runStats = runStats;
    }

    /**
     * Phases the current window of genotype data.
     * @param cd the current window of data
     * @return the phased genotype data
     * @throws IllegalArgumentException if
     * {@code gv != null && gv.markers().equals(cd.targMarkers() == false)}
     * @throws IllegalArgumentException if
     * {@code gv != null && gv.samples().equals(cd.targetSamples() == false)}
     * @throws NullPointerException if {@code cd == null}
     */
    GT phase(CurrentData cd) {
        return cd.targGT();
        //TODO: check this out
        /*if (cd.targGT().isPhased()) {
            return cd.targGT();
        }
        EstPhase.HapsGT hapsGT = lsPhaseSingles(cd);
        //phaseParentOffspring(cd, hapPairs);
        return hapsGT;*/
    }

    private EstPhase.HapsGT lsPhaseSingles(CurrentData cd) {
        float minAlleleFreq = 0.0001f;
        List<LongArray> initHaps = InitTargHapPairs.run(cd.targGT(),
                cd.restrictRefGT(), minAlleleFreq, par.seed());
        EstPhase estPhase = new EstPhase(cd.targGT(), initHaps);
        Random rand = new Random(par.seed());
        int nBurninIts = par.burnin();
        int nIts = nBurninIts + par.iterations();

        GT gt = cd.targGT();
        double recombFactor = (float) 0.04f*par.ne()/(2*cd.nAllSamples());
        for (int it=0; it<nIts; ++it) {
            long t0 = System.nanoTime();
            boolean burnin = it<nBurninIts;
            boolean updateRecomb = (it>(nBurninIts-3) && it<nBurninIts);
            Regress recombRegress = updateRecomb ? new Regress() : null;
            PhaseData phaseData = new PhaseData(cd, estPhase, recombFactor, it,
                    rand.nextLong());
            PhaseLS.run(phaseData, estPhase, recombRegress);
            if (recombRegress!=null && recombRegress.cnt()>=100) {
                recombFactor = recombFactor(cd, recombRegress);
            }
            long elapsedNanos = System.nanoTime() - t0;
            runStats.printPhasingIterationUpdate(it+1, burnin, elapsedNanos);
        }
        return estPhase.hapsGT();
    }

    private double recombFactor(CurrentData cd, Regress recombRegress) {
        double recombFactor = recombRegress.beta();
        float initNe = cd.par().ne();
        float maxNe = Math.max(initNe, 5e7f);
        final double maxRecombFactor = (float) 0.04f*maxNe/(2*cd.nAllSamples());
        if (recombFactor <= 0.0 || recombFactor > maxRecombFactor) {
            runStats.println("");
            runStats.println("WARNING: recombination factor estimate is out-of-bounds: "
                    + recombFactor);
            runStats.println("");
            recombFactor = (float) 0.04f*initNe/(2*cd.nAllSamples());
        }
        return recombFactor;
    }

    public static <E> List<E> toList(AtomicReferenceArray<E> ara) {
        List<E> list = new ArrayList<>(ara.length());
        for (int j = 0, n = ara.length(); j < n; ++j) {
            list.add(ara.get(j));
        }
        return list;
    }
}
