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
package org.ngseq.sparkbeagle.sample;

import org.ngseq.sparkbeagle.ints.LongArray;
import org.ngseq.sparkbeagle.vcf.GT;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.IntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * <p>Class {@code InitTargHapPairs} has a static method for returning
 * initial target haplotype pairs.
 * </p>
 * <p>Instances of class {@code InitTargHapPairs} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class InitTargHapPairs {

    private InitTargHapPairs() {
        // private constructor to prevent instantiation
    }

    /**
     * Returns initial target haplotype pairs.  Missing alleles are sampled
     * according from the empirical allele frequency distribution, and input
     * genotypes with unknown phase will have their alleles randomly ordered
     * in the returned haplotype pairs.
     * @param targ the target genotype data
     * @param ref the reference haplotype pairs
     * @param minFreq a minimum allele frequency to be enforced
     * @param seed the seed for random number generation
     * @return initial target haplotype pairs
     * @throws IllegalArgumentException if
     * {@code minFreq <= 0.0 || minFreq >= 0.5 || Double.isNaN(minFreq) == true}
     * @throws NullPointerException if {@code gl == null || refHapPairs == null}
     */
    public static List<LongArray> run(GT targ, GT ref, double minFreq,
            long seed) {
        if (minFreq <= 0.0 || minFreq >= 0.5 || Double.isNaN(minFreq)) {
            throw new IllegalArgumentException(String.valueOf(minFreq));
        }
        double[][] alFreq = alleleFreq(targ, ref, minFreq);
        IntFunction<List<LongArray>> mapper = mapper(targ, alFreq, seed);
        return IntStream.range(0, targ.nSamples())
                .parallel()
                .mapToObj(mapper)
                .flatMap(list -> list.stream())
                .collect(Collectors.toCollection(ArrayList::new));
    }

    private static double[][] alleleFreq(GT targ, GT ref,
            double minFreq) {
        IntFunction<double[]> freqMapper = freqMapper(targ, ref, minFreq);
        return IntStream.range(0, targ.nMarkers())
                .parallel()
                .mapToObj(freqMapper)
                .toArray(double[][]::new);
    }

    private static IntFunction<double[]> freqMapper(GT targ, GT ref,
            double minFreq) {
        return (marker) -> {
            if (targ.isGTData()) {
                return alleleFreqGT(targ, ref, marker, minFreq);
            }
            else {
                return alleleFreqGL(targ, ref, marker, minFreq);
            }
        };
    }

    private static double[] alleleFreqGL(GT targ, GT ref, int marker,
            double minFreq) {
        int nAlleles = targ.marker(marker).nAlleles();
        double[] alleleFreq = new double[nAlleles];
        double[] scaledFreq = new double[nAlleles];
        for (int s=0, n=targ.nSamples(); s<n; ++s) {
            for (int a2=0; a2<nAlleles; ++a2) {
                for (int a1=0; a1<=a2; ++a1) {
                    double like = targ.gl(marker, s, a1, a2);
                    if (a1!=a2) {
                        like = Math.max(like, targ.gl(marker, s, a2, a1));
                    }
                    scaledFreq[a2] += like;
                    scaledFreq[a1] += like;
                }
            }
            divideElementsBySum(scaledFreq);
            for (int j=0; j<scaledFreq.length; ++j) {
                alleleFreq[j] += scaledFreq[j];
                scaledFreq[j] = 0.0f;
            }
        }
        if (ref != null) {
            for (int h=0, n=ref.nHaps(); h<n; ++h) {
                ++alleleFreq[ref.allele(marker, h)];
            }
        }
        divideElementsBySum(alleleFreq);
        enforceMinFreq(alleleFreq, minFreq);
        return alleleFreq;
    }

    private static double[] alleleFreqGT(GT targ, GT ref, int marker,
            double minFreq) {
        int nAlleles = targ.marker(marker).nAlleles();
        int[] cnts = new int[nAlleles];
        double[] freq = new double[nAlleles];
        for (int h=0, n=targ.nHaps(); h<n; ++h) {
            int allele = targ.allele(marker, h);
            if (allele>=0) {
                ++cnts[allele];
            }
        }
        if (ref != null) {
            for (int h=0, n=ref.nHaps(); h<n; ++h) {
                int allele = ref.allele(marker, h);
                ++cnts[allele];
            }
        }
        int sum = Arrays.stream(cnts).sum();
        for (int al=0; al<cnts.length; ++al) {
            freq[al] = (double) cnts[al] / sum;
        }
        enforceMinFreq(freq, minFreq);
        return freq;
    }

    private static void divideElementsBySum(double[] da) {
        double sum = Arrays.stream(da).sum();
        for (int j=0; j<da.length; ++j) {
            da[j] /= sum;
        }
    }

    private static void enforceMinFreq(double[] alleleFreq, double minFreq) {
        boolean changedFreq = false;
        for (int j=0; j<alleleFreq.length; ++j) {
            if (alleleFreq[j] < minFreq) {
                alleleFreq[j] = minFreq;
                changedFreq = true;
            }
        }
        if (changedFreq) {
            divideElementsBySum(alleleFreq);
        }
    }

    private static IntFunction<List<LongArray>> mapper(GT gt,
            double[][] alFreq, long seed) {
        return (sample) -> {
            int nMarkers = gt.nMarkers();
            Random rand = new Random(seed + sample);
            List<LongArray> list = new ArrayList<>(2);
            int[] hap1 = new int[gt.nMarkers()];
            int[] hap2 = new int[gt.nMarkers()];
            boolean isGT = gt.isGTData();
            for (int m=0; m<nMarkers; ++m) {
                if (isGT) {
                    sampleGT(gt, m, sample, alFreq, hap1, hap2, rand);
                }
                else {
                    sampleGL(gt, m, sample, alFreq, hap1, hap2, rand);
                }
            }
            list.add(gt.markers().allelesToBits(hap1));
            list.add(gt.markers().allelesToBits(hap2));
            return list;
        };
    }

    private static void sampleGT(GT gl, int m, int sample,
            double[][] alleleFreq, int[] hap1, int[]  hap2, Random rand) {
        boolean isPhased = gl.isPhased(m, sample);
        boolean swap = isPhased ? false : rand.nextBoolean();
        int a1 = swap ? gl.allele2(m, sample) : gl.allele1(m, sample);
        int a2 = swap ? gl.allele1(m, sample) : gl.allele2(m, sample);
        hap1[m] = a1<0 ? randAllele(alleleFreq[m], rand) : a1;
        hap2[m] = a2<0 ? randAllele(alleleFreq[m], rand) : a2;
    }

    private static void sampleGL(GT gl, int m, int sample,
            double[][] alleleFreq, int[] hap1, int[]  hap2, Random rand) {
        int a1 = randAllele(alleleFreq[m], rand);
        int a2 = randAllele(alleleFreq[m], rand);
        while (gl.gl(m, sample, a1, a2)==0.0) {
            a1 = randAllele(alleleFreq[m], rand);
            a2 = randAllele(alleleFreq[m], rand);
        }
        hap1[m] = a1;
        hap2[m] = a2;
    }

    private static int randAllele(double[] freq, Random rand) {
        double d = rand.nextDouble();
        double sum = 0.0;
        for (int j=0, n=freq.length; j<n; ++j) {
            sum += freq[j];
            if (sum >= d) {
                return j;
            }
        }
        return freq.length - 1;
    }
}
