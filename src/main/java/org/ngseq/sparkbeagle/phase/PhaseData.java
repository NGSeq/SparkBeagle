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

import org.ngseq.sparkbeagle.CurrentData;
import org.ngseq.sparkbeagle.Par;
import org.ngseq.sparkbeagle.blbutil.FloatArray;
import org.ngseq.sparkbeagle.vcf.Marker;
import org.ngseq.sparkbeagle.vcf.Markers;
import org.ngseq.sparkbeagle.vcf.RefGT;

/**
 * <p>Class {@code PhaseData} contains the input data for phasing
 * genotypes.
 * </p>
 * <p>Instances of class {@code PhaseData} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PhaseData {

    private final CurrentData cd;
    private final int it;
    private final int nIts;
    private final boolean burnin;
    private final RefGT ref;
    private final EstPhase.HapsGT targHaps;
    private final int nTargHaps;
    private final FloatArray pRecomb;
    private final double[] pos;
    private final long seed;

    /**
     * Constructs a new {@code ImpData} instance from the specified data.
     * @param cd the input data for the current marker window
     * @param estPhase the current estimate of phased target genotypes
     * @param recombFactor the factor multiplied by genetic distance to
     * obtain the probability of transitioning to a random HMM state.
     * @param it the current iteration (first iteration has index 0)
     * @param seed seed for random numbers
     *
     * @throws IllegalArgumentException if
     * {@code cd.targMarkers().equals(estPhase.markers() == false}
     * @throws IllegalArgumentException if
     * {@code recombFactor < 0 || Double.isFinite(recombFactor)==false}
     * @throws IllegalArgumentException if
     * {@code cd.targSamples().equals(estPhase.samples()) == false}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public PhaseData(CurrentData cd, EstPhase estPhase, double recombFactor,
            int it, long seed) {
        if (cd.targSamples().equals(estPhase.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (cd.targMarkers().equals(estPhase.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (recombFactor < 0 || Double.isFinite(recombFactor)==false) {
            throw new IllegalArgumentException(String.valueOf(recombFactor));
        }
        Par par = cd.par();
        this.cd = cd;
        this.it = it;
        this.nIts = par.burnin() + par.iterations();
        this.burnin = it < par.burnin();
        this.ref = cd.restrictRefGT();
        this.targHaps = estPhase.hapsGT();
        this.nTargHaps = cd.nTargHaps();
        this.pos = cd.map().genPos();
        this.pRecomb = pRecomb(cd.genDist(), recombFactor);
        this.seed = seed;
    }

    private static FloatArray pRecomb(FloatArray genDist, double recombFactor) {
        float[] pRecomb = new float[genDist.size()];
        for (int m=0; m<pRecomb.length; ++m) {
            pRecomb[m] = (float) -Math.expm1(-recombFactor * genDist.get(m));
        }
        return new FloatArray(pRecomb);
    }

    /**
     * Returns the command line parameters
     * @return the command line parameters
     */
    public Par par() {
        return cd.par();
    }

    /**
     * Returns {@code true} if the current phasing iteration is a burnin
     * iteration and returns {@code false} otherwise.
     * @return {@code true} if the current phasing iteration is a burnin
     * iteration
     */
    public boolean burnin() {
        return burnin;
    }

    public int allele(int marker, int hap) {
        if (hap < nTargHaps) {
            return targHaps.allele(marker, hap);
        }
        else {
            return ref.allele(marker, hap - nTargHaps);
        }
    }

    /**
     * Returns the iteration index.
     * @return the iteration index
     */
    public int iter() {
        return it;
    }

    /**
     * Returns the number of iterations remaining.
     * @return the number of iterations remaining
     */
    public int nItsRemaining() {
        return nIts - it;
    }

    /**
     * Returns the number of target markers
     * @return the number of target markers
     */
    public int nMarkers() {
        return cd.nTargMarkers();
    }

    /**
     * Returns the list of target markers.
     * @return the list of target markers
     */
    public Markers markers() {
        return cd.targMarkers();
    }

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public Marker marker(int marker) {
        return cd.targMarkers().marker(marker);
    }

    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */
    public int nTargSamples() {
        return cd.nTargetSamples();
    }

   /**
     * Return the number of reference haplotypes.
     * @return the number of reference haplotypes
     */
    public int nRefHaps() {
        return cd.nRefHaps();
    }

    /**
     * Return the number of target haplotypes.
     * @return the number of target haplotypes
     */
    public int nTargHaps() {
        return cd.nTargHaps();
    }

    /**
     * Returns the number of reference and target haplotypes.
     * @return the number of reference and target haplotypes
     */
    public int nHaps() {
        return cd.nHaps();
    }

    /**
     * Returns the probability that the allele carried by the specified
     * target marker cluster matches the allele labeling the latent HMM state.
     * @param marker index of a target marker cluster
     * @return the probability that the allele carried by the specified
     * target marker cluster matches the allele labeling the latent HMM state
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public float pErr(int marker) {
        return cd.par().err();
    }

    /**
     * Returns the array of genetic map positions whose {@code k}-th element
     * equals {@code this.pos(k)}.
     * @return the array of genetic map positions
     */
    public double[] pos() {
        return pos.clone();
    }

    /**
     * Return a {@code FloatArray} of size {@code this.nTargMarkers()}
     * whose {@code k}-th element is the genetic distance between the
     * {@code k}-th target marker and the previous marker, or {@code 0.0}
     * if {@code (k == 0)}.
     * @return a {@code FloatArray} of size {@code this.nTargMarkers()}
     * whose {@code k}-th element is the genetic distance between the
     * {@code k}-th target marker and the previous marker,
     */
    public FloatArray genDist() {
        return cd.genDist();
    }

    /**
     * Return a {@code FloatArray} of size {@code this.nTargMarkers()}
     * whose {@code k}-th element is the factor multiplied by genetic distance
     * to obtain the probability of transitioning to a random HMM state
     * between the {@code k}-th target marker and the previous marker.
     * @return a {@code FloatArray} of size {@code this.nTargMarkers()}
     * whose {@code k}-th element is the factor multiplied by genetic distance
     * to obtain the probability of transitioning to a random HMM state
     * between the {@code k}-th target marker and the previous marker
     */
    public FloatArray pRecomb() {
        return pRecomb;
    }

    /**
     * Returns the seed for generating random numbers.
     * @return the seed for generating random numbers
     */
    public long seed() {
        return seed;
    }
}
