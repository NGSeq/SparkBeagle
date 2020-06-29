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

import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.vcf.RefGT;
import org.ngseq.sparkbeagle.vcf.RefGTRec;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.stream.IntStream;

/**
 * <p>Class {@code RefHapHash} stores a hash code for each haplotype
 * in a sublist of reference haplotypes. The hash code is computed
 * from the allele sequence carried by the reference haplotype.
 * </p>
 * <p>Instances of class {@code RefHapHash} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefHapHash {

    private final int targMarker;
    private final RefGT refGT;
    private final int[] i2hap;
    private final IntList[] altAlleles; // marker offsets and ALT alleles
    private final int[] i2hash;
    private final int start;
    private final int end;

    /**
     * Constructs a new {@code RefHapHash} instance for the specified data.
     * The sublist of reference haplotypes is the ordered list of distinct
     * reference haplotypes with stored state probability data at the specified
     * target marker in the {@code stateProbs} parameter.
     * @param stateProbs HMM state probabilities at the genotyped
     * markers in the target samples
     * @param targMarker a target marker index
     * @param refHapPairs the reference haplotypes
     * @param start the starting reference marker index (inclusive) for
     * the reference haplotype allele sequences
     * @param end the ending reference marker index (exclusive) for the
     * reference haplotype allele sequences
     * @throws IndexOutOfBoundsException if {@code targMarker < 0}
     * @throws IndexOutOfBoundsException if there exists a {@code j} satisfying
     * {@code (0 <= j && j < stateProbs.size()) &&
     * (targMarker >= stateProbs.get(j).nMarkers())}
     * @throws IllegalArgumentException if
     * {@code start < 0 || start >= end || end > refHapPairs.nMarkers()}
     * @throws NullPointerException if
     * {@code stateProbs == null || refHapPairs == null}
     * @throws NullPointerException if there exists a {@code j} satisfying
     * {@code (0 <= j && j < stateProbs.size()) && (stateProbs.get(j) == null)}
     */
    public RefHapHash(AtomicReferenceArray<StateProbs> stateProbs,
            int targMarker, RefGT refHapPairs, int start, int end) {
        if (start < 0 || start >= end) {
            throw new IllegalArgumentException(String.valueOf(start));
        }
        if (end > refHapPairs.nMarkers()) {
            throw new IllegalArgumentException(String.valueOf(end));
        }
        this.targMarker = targMarker;
        this.refGT = refHapPairs;
        this.i2hap = i2hap(stateProbs, targMarker);
        this.i2hash = new int[i2hap.length];
        this.start = start;
        this.end = end;

        this.altAlleles = IntStream.range(0, i2hap.length)
                .mapToObj(i -> new IntList(6))
                .toArray(IntList[]::new);
        setHashAndAltAlleles();
    }

    private void setHashAndAltAlleles() {
        Random rand = new Random(start);
        for (int m=start; m<end; ++m) {
            RefGTRec rec = refGT.get(m);
            int markerOffset = m - start;
            if (rec.isAlleleCoded() && rec.majorAllele()==0) {
                lowAltFreqUpdate(rec, markerOffset, rand);
            }
            else {
                standardUpdate(rec, markerOffset, rand);
            }
        }
    }

    private void lowAltFreqUpdate(RefGTRec rec, int markerOffset, Random rand) {
        int nAlleles = rec.nAlleles();
        assert rec.majorAllele()==0;
        for (int al=1; al<nAlleles; ++al) {
            int hash = rand.nextInt();
            int nCopies = rec.alleleCount(al);
            if (i2hap.length < nCopies) {
                for (int i=0; i<i2hap.length; ++i) {
                    if (rec.isCarrier(al, i2hap[i])) {
                        i2hash[i] += hash;
                        altAlleles[i].add(markerOffset);
                        altAlleles[i].add(al);
                    }
                }
            }
            else {
                for (int c=0; c<nCopies; ++c) {
                    int hap = rec.hapIndex(al, c);
                    int i = Arrays.binarySearch(i2hap, hap);
                    if (i>=0) {
                        i2hash[i] += hash;
                        altAlleles[i].add(markerOffset);
                        altAlleles[i].add(al);
                    }
                }
            }
        }
    }

    private void standardUpdate(RefGTRec rec, int markerOffset, Random rand) {
        int[] alleleHash = new int[rec.nAlleles()];
        for (int a=1; a<alleleHash.length; ++a) {
            alleleHash[a] = rand.nextInt();
        }
        for (int i=0; i<i2hap.length; ++i) {
            int allele = rec.get(i2hap[i]);
            if (allele!=0) {
                i2hash[i] += alleleHash[allele];
                altAlleles[i].add(markerOffset);
                altAlleles[i].add(allele);
            }
        }
    }

    private static int[] i2hap(AtomicReferenceArray<StateProbs> stateProbs,
            int targMarker) {
        IntList list = new IntList(10*stateProbs.length());
        for (int j=0, n=stateProbs.length(); j<n; ++j) {
            StateProbs sp = stateProbs.get(j);
            for (int k=0, m=sp.nStates(targMarker); k<m; ++k) {
                list.add(sp.refHap(targMarker, k));
            }
        }
        return list.stream().sorted().distinct().toArray();
    }

    /**
     * Returns the target marker.
     * @return the target marker
     */
    public int targMarker() {
        return targMarker;
    }

    /**
     * Returns the phased reference genotypes.
     * @return the phased reference genotypes
     */
    public RefGT refGT() {
        return refGT;
    }

    /**
     * Returns the starting reference marker index (inclusive)
     * for the reference haplotype allele sequences.
     * @return the starting reference marker index (inclusive)
     * for the reference haplotype allele sequences
     */
    public int start() {
        return start;
    }

   /**
     * Returns the ending reference marker index (exclusive)
     * for the reference haplotype allele sequences.
     * @return the ending reference marker index (exclusive)
     * for the reference haplotype allele sequences
     */
    public int end() {
        return end;
    }

    /**
     * Copies the alleles between {@code this.start()} (inclusive) and
     * {@code this.end()} (exclusive) of the specified haplotype in the sublist
     * of reference haplotypes to the specified array.
     * @param index an index in the sublist of reference haplotypes
     * @param alleles the array to which haplotype alleles will be copied
     * @throws IllegalArgumentException if
     * {@code alleles.length < (this.end() - this.start())}
     * @throws NullPointerException if {@code alleles == null}
     */
    public void setAlleles(int index, int[] alleles) {
        Arrays.fill(alleles, 0, (end-start), 0);
        IntList il = altAlleles[index];
        for (int j=0, n=il.size(); j<n; j+=2) {
            alleles[il.get(j)] = il.get(j+1);
        }
    }

    /**
     * Returns the size of the sublist of reference haplotypes.
     * @return the size of the sublist of reference haplotypes
     */
    public int nHaps() {
        return i2hap.length;
    }

    /**
     * Returns the index of the specified haplotype in the sublist of
     * reference haplotypes. Returns (-insertionPoint - 1) if the
     * specified haplotype is not found in the sublist of reference
     * haplotypes.
     * @param hap a reference haplotype index
     * @return the index of the specified haplotype in the sublist of
     * reference haplotypes
     */
    public int hap2Index(int hap) {
        return Arrays.binarySearch(i2hap, hap);
    }

    /**
     * Return the specified haplotype in the sublist of reference haplotypes.
     * @param index an index in the sublist of reference haplotypes
     * @return the specified haplotype in the sublist of reference haplotypes
     */
    public int hap(int index) {
        return i2hap[index];
    }

    /**
     * Return a hash code computed from the allele sequence
     * of the specified haplotype in the sublist of reference haplotypes.
     * The hash code is computed from the allele sequence between
     * markers with index {@code this.start()} (inclusive) and
     * {@code this.end()} (exclusive).
     * @param index an index in the sublist of reference haplotypes
     * @return a hash code computed from the allele sequence
     * of the specified haplotype in the sublist of reference haplotypes
     */
    public int hash(int index) {
        return i2hash[index];
    }
}
