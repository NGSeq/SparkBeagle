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
package org.ngseq.sparkbeagle.vcf;

import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.blbutil.Const;

import java.util.Arrays;

/**
 * <p>Class {@code BasicGT} represents genotype emission probabilities
 * for a set of samples.
 * </p>
 * Instances of class {@code BasicGT} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BasicGT implements GT {

    private final Samples samples;
    private final Markers markers;
    private final GTRec[] vma;
    private final boolean isRefData;
    private final boolean isGTData;
    private final boolean[] isRefSample;

    /**
     * Returns the genotype index corresponding to the
     * specified unordered alleles.
     * @param a1 the first allele index of an unordered genotype
     * @param a2 the second allele index of an unordered genotype
     * @return the genotype index corresponding to the
     * specified unordered alleles
     * @throws IllegalArgumentException if {@code a1 < 0 || a2 < 0}
     */
    public static int genotype(int a1, int a2) {
        if (a1<=a2) {
            if (a1 < 0) {
                String s = "allele < 0: " + a1 + " " + a2;
                throw new IllegalArgumentException(s);
            }
            return (a2*(a2+1))/2 + a1;
        }
        else {
            if (a2<0) {
                String s = "allele < 0: " + a1 + " " + a2;
                throw new IllegalArgumentException(s);
            }
            return (a1*(a1+1))/2 + a2;
        }
    }

    /**
     * Constructs a {@code BasicGT} instance.
     *
     * @param samples the list of samples with genotype data
     * @param vma genotype emission probabilities
     *
     * @throws IllegalArgumentException
     * if elements of {@code vma} corresponding to the same chromosome
     * are not contiguous and sorted in chromosome position order
     * @throws IllegalArgumentException if any
     * two {@code vma} elements correspond to the same genetic marker
     * @throws IllegalArgumentException if
     * {@code vma[j].samples().equals(samples) == false} for any {@code j}
     * satisfying {@code 0 <= j && j < vma.length}
     *
     * @throws NullPointerException if {@code samples == null}
     * @throws NullPointerException if {@code vma == null}
     * @throws NullPointerException if {@code vma[j] == null} any {@code j}
     * satisfying {@code 0 <= j && j < vma.length}
     */
    public BasicGT(Samples samples, GTRec[] vma) {
        checkSamples(samples, vma);
        this.markers = markers(vma);
        this.samples = samples;
        this.vma = vma.clone();
        this.isRefSample = isRefSample(samples, vma);
        this.isRefData = isRefData(vma);
        this.isGTData = isGTData(vma);
    }

    private static void checkSamples(Samples samples, GTRec[] mla) {
        for (int j=0; j<mla.length; ++j) {
            if (mla[j].samples().equals(samples)==false) {
                throw new IllegalArgumentException("inconsistent samples");
            }
        }
    }

    private static Markers markers(GTRec[] vma) {
        Marker[] markers = new Marker[vma.length];
        for (int j=0; j<markers.length; ++j) {
            markers[j] = vma[j].marker();
        }
        return Markers.create(markers);
    }

    private static boolean isRefData(GTRec[] vma) {
        boolean isRefData = true;
        for (int j=0; j<vma.length && isRefData==true; ++j) {
            if (vma[j].isPhased()==false) {
                isRefData = false;
            }
        }
        return isRefData;
    }

    private static boolean isGTData(GTRec[] vma) {
        boolean isGTData = true;
        for (int j=0; j<vma.length && isGTData==true; ++j) {
            if (vma[j].isGTData()==false) {
                isGTData = false;
            }
        }
        return isGTData;
    }

    private static boolean[] isRefSample(Samples samples, GTRec[] vea) {
        boolean[] isRefSample = new boolean[samples.nSamples()];
        Arrays.fill(isRefSample, true);
        for (int s=0; s<isRefSample.length; ++s) {
            for (GTRec ve : vea) {
                if (ve.isPhased(s)==false) {
                    isRefSample[s] = false;
                    break;
                }
            }
        }
        return isRefSample;
    }

    @Override
    public boolean isPhased() {
        return isRefData;
    }

    @Override
    public boolean isGTData() {
        return isGTData;
    }

    @Override
    public boolean isPhased(int sample) {
        return isRefSample[sample];
    }

    @Override
    public float gl(int marker, int sample, int allele1, int allele2) {
        return vma[marker].gl(sample, allele1, allele2);
    }

    @Override
    public boolean isPhased(int marker, int sample) {
        return vma[marker].isPhased(sample);
    }

    @Override
    public int allele1(int marker, int sample) {
        return vma[marker].allele1(sample);
    }

    @Override
    public int allele2(int marker, int sample) {
        return vma[marker].allele2(sample);
    }

    @Override
    public int allele(int marker, int hap) {
        return vma[marker].get(hap);
    }

    @Override
    public int nMarkers() {
        return vma.length;
    }

    @Override
    public Marker marker(int markerIndex) {
        return markers.marker(markerIndex);
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public int nHaps() {
        return 2*samples.nSamples();
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public String toString() {
        StringBuilder sb  = new StringBuilder();
        sb.append("[BasicGL: nMarkers=");
        sb.append(nMarkers());
        sb.append(" nSamples=");
        sb.append(nSamples());
        for (GTRec vm : vma) {
            sb.append(Const.nl);
            sb.append(vm);
        }
        sb.append(']');
        return sb.toString();
    }
}
