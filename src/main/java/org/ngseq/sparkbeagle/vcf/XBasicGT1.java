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

import org.ngseq.sparkbeagle.blbutil.Const;

import java.util.BitSet;

/**
 * <p>Class {@code XBasicGT1} represents genotype likelihoods for one sample.
 * </p>
 * <p>Instances of class {@code XBasicGT1} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class XBasicGT1 implements XGT1 {

    private final int idIndex;
    private final Markers markers;
    private final BitSet missing1;
    private final BitSet missing2;
    private final BitSet allele1;
    private final BitSet allele2;
    private final BitSet isPhased;
    private final boolean isRefSample;

    /**
     * Constructs a {@code XBasicGL1} instance from the specified data.
     *
     * @param gl the genotype likelihoods
     * @param sample the sample index
     *
     * @throws IndexOutOfBoundsException if
     * {@code samples < 0 || sample >= gl.nSamples()}
     * @throws NullPointerException if {@code gl == null}
     */
    public XBasicGT1(GT gl, int sample) {
        this(gl, sample, -1, -1);
    }

    /**
     * Constructs a {@code XBasicGL1} instance from the specified data.
     *
     * @param gt the genotype data
     * @param sample the sample index
     * @param father the sample index of the sample's father, or -1 if the
     * father is not genotyped
     * @param mother the sample index of the sample's mother, or -1 if the
     * mother is not genotyped
     *
     * @throws IndexOutOfBoundsException if
     * {@code samples < 0 || sample >= gl.nSamples()}
     * @throws NullPointerException if {@code gl == null}
     */
    public XBasicGT1(GT gt, int sample, int father, int mother) {
        // NB: phase information in GL is ignored if parental data is used
        int nMarkers = gt.markers().nMarkers();
        int sumHapBits = gt.markers().sumHaplotypeBits();
        this.idIndex = gt.samples().idIndex(sample);
        this.markers = gt.markers();
        this.missing1 = new BitSet(nMarkers);
        this.missing2 = new BitSet(nMarkers);
        this.allele1 = new BitSet(sumHapBits);
        this.allele2 = new BitSet(sumHapBits);
        this.isPhased = new BitSet(nMarkers);
//      // Removing use of pedigree constraints to preserve regression
//      //   after changing main.Par.ped() to return pedigree
//        if (father==-1 && mother==-1) {
//            this.isRefSample = copyGL(gl, sample);
//        }
//        else {
//            this.isRefSample = useConstraintAndGL(gl, sample, father, mother);
//        }
        this.isRefSample = copyGL(gt, sample);
    }

    private boolean copyGL(GT gl, int sample) {
        int nMarkers = markers.nMarkers();
        for (int m=0; m<nMarkers; ++m) {
            int a1 = gl.allele1(m, sample);
            int a2 = gl.allele2(m, sample);
            setBits(markers, m, a1, allele1, missing1);
            setBits(markers, m, a2, allele2, missing2);
            if (gl.isPhased(m, sample)) {
                isPhased.set(m);
            }
        }
        return gl.isPhased(sample);
    }

    private boolean useConstraintAndGL(GT gl, int sample, int father,
            int mother) {
        int nMarkers = markers.nMarkers();
        int phasedCnt = 0;
        for (int m=0; m<nMarkers; ++m) {
            int a1 = gl.allele1(m, sample);
            int a2 = gl.allele2(m, sample);
            int tr1 = transmittedAllele(gl, m, father);
            int tr2 = transmittedAllele(gl, m, mother);
            if ((tr1>=0 || tr2>=0) && isUnphasedConsistent(a1, a2, tr1, tr2)) {
                if (isPhasedConsistent(a1, a2, tr1, tr2)==false) {
                    int tmp = a1;
                    a1 = a2;
                    a2 = tmp;
                }
                if (tr1 != -1) {
                    a1 = tr1;
                }
                if (tr2 != -1) {
                    a2 = tr2;
                }
                if (a1>=0 && a2>=0) {
                    ++phasedCnt;
                    isPhased.set(m);
                }
            }
            setBits(markers, m, a1, allele1, missing1);
            setBits(markers, m, a2, allele2, missing2);
        }
        return phasedCnt==nMarkers;
}

    private static boolean isUnphasedConsistent(int a1, int a2, int b1, int b2) {
        return isPhasedConsistent(a1, a2, b1, b2)
                || isPhasedConsistent(a1, a2, b2, b1);
    }

    private static boolean isPhasedConsistent(int a1, int a2, int b1, int b2) {
        return (isConsistent(a1, b1) && isConsistent(a2, b2));
    }

    private static boolean isConsistent(int a1, int b1) {
        return (a1==-1 || b1==-1 || a1==b1);
    }

    private static int transmittedAllele(GT gl, int marker, int sample) {
        if (sample==-1) {
            return -1;
        }
        int a1 = gl.allele1(marker, sample);
        return a1>=0 && a1==gl.allele2(marker, sample) ? a1 : -1;
    }

    private static void setBits(Markers markers, int marker, int allele,
            BitSet alleles, BitSet missing) {
        if (allele == -1) {
            missing.set(marker);
        }
        else {
            int mask = 1;
            int start = markers.sumHaplotypeBits(marker);
            int end = markers.sumHaplotypeBits(marker+1);
            for (int i=start; i<end; ++i) {
                boolean b = (allele & mask)==mask;
                alleles.set(i, b);
                mask <<= 1;
            }
        }
    }

    private int allele(BitSet bitset, int marker) {
        int start = markers.sumHaplotypeBits(marker);
        int end = markers.sumHaplotypeBits(marker+1);
        if (end==(start+1)) {
            return bitset.get(start) ? 1 : 0;
        }
        int allele = 0;
        int mask = 1;
        for (int j=start; j<end; ++j) {
            if (bitset.get(j)) {
                allele += mask;
            }
            mask <<= 1;
        }
        return allele;
    }

    @Override
    public boolean isRefSample() {
        return isRefSample;
    }

    @Override
    public float gl(int marker, int a1, int a2) {
        if (a1 < 0 || a1 >= markers.marker(marker).nAlleles()) {
            throw new IndexOutOfBoundsException(String.valueOf(a1));
        }
        if (a2 < 0 || a2 >= markers.marker(marker).nAlleles()) {
            throw new IndexOutOfBoundsException(String.valueOf(a2));
        }
        int obsA1 = allele1(marker);
        int obsA2 = allele2(marker);
        boolean consistent = (obsA1==-1 || obsA1==a1) && (obsA2==-1 || obsA2==a2);
        if (consistent==false && isPhased.get(marker)==false) {
            consistent = (obsA1==-1 || obsA1==a2) && (obsA2==-1 || obsA2==a1);
        }
        return consistent ? 1.0f : 0.0f;
    }

    @Override
    public boolean isPhased(int marker) {
        return isPhased.get(marker);
    }

    @Override
    public int allele1(int marker) {
        return missing1.get(marker) ? -1 : allele(allele1, marker);
    }

    @Override
    public int allele2(int marker) {
        return missing2.get(marker) ? -1 : allele(allele2, marker);
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
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
    public int idIndex() {
        return idIndex;
    }

    @Override
    public String toString() {
        int nMarkers = markers.nMarkers();
        StringBuilder sb  = new StringBuilder();
        sb.append("[XBasicGL1: nMarkers=");
        sb.append(nMarkers);
        sb.append(Const.nl);
        for (int m=0; m<nMarkers; ++m) {
            sb.append(markers.marker(m));
            sb.append(Const.tab);
            sb.append(allele1(m));
            sb.append(isPhased(m) ? Const.phasedSep : Const.unphasedSep);
            sb.append(allele2(m));
            sb.append(Const.nl);
        }
        sb.append(']');
        return sb.toString();
    }
}
