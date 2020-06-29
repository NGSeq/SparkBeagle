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

import org.ngseq.sparkbeagle.ints.IndexArray;
import org.ngseq.sparkbeagle.ints.IntArray;
import org.ngseq.sparkbeagle.vcf.GT;
import org.ngseq.sparkbeagle.vcf.RefGT;
import org.ngseq.sparkbeagle.vcf.RefGTRec;

import java.util.Arrays;

/**
 * <p>Class {@code HaplotypeCoder} indexes the observed allele sequences
 * in phased reference and target genotype data in a chromosome interval.
 * </p>
 * <p>Instances of class {@code HaplotypeCoder} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HaplotypeCoder {

    private final int nRefHaps;
    private final int nHaps;
    private final RefGT ref;
    private final GT targ;

    /**
     * Constructs a new {@code HaplotypeCoder} instance from the specified
     * data.
     * @param restrictRefGT the phased reference genotypes at the target
     * markers
     * @param phasedTarg the phased target genotypes
     * @throws IllegalArgumentException if
     * {@code refHapPairs.markers().equals(targetHapPairs.markers()) == false}
     * @throws NullPointerException if
     * {@code refHapPairs == null || targetHapPairs == null}
     */
    public HaplotypeCoder(RefGT restrictRefGT, GT phasedTarg) {
        if (restrictRefGT.markers().equals(phasedTarg.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        this.nRefHaps = restrictRefGT.nHaps();
        this.nHaps = nRefHaps + phasedTarg.nHaps();
        this.ref = restrictRefGT;
        this.targ = phasedTarg;
    }

    /**
     * Returns the phased reference genotypes at the target markers.
     * @return the phased reference genotypes at the target markers
     */
    public RefGT refHapPairs() {
        return ref;
    }

     /**
     * Returns the phased target genotypes.
     * @return the phased target genotypes
     */
    public GT targHapPairs() {
        return targ;
    }

    /**
     * Returns an array mapping haplotype indices to allele sequence indices
     * in the specified marker interval.
     * The target haplotype indices are shifted by the number of reference
     * haplotypes in the returned array. Thus the first target haplotype
     * index is {@code this.refHapPairs.nHaps()} and the last target
     * haplotype index is
     * {@code (this.refHapPairs.nHaps() + this.targHapPairs().nHaps() - 1)}.
     *
     * @param start the first marker index (inclusive)
     * @param end the last marker index (exclusive)
     * @return an array mapping haplotype indices to allele sequence indices
     * @throws IllegalArgumentException if {@code start >= end}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end >= this.refHapPairs.nMarkers()}
     */
    public IndexArray run(int start, int end) {
        if (start >= end) {
            throw new IllegalArgumentException("start >= end");
        }
        if (isHapCoded(ref, start, end)) {
            return codeSeqCodedRef(start, end);
        }
        else {
            return codeSeq(start, end);
        }
    }

    private static boolean isHapCoded(RefGT ref, int start, int end) {
        RefGTRec startRec = ref.get(start);
        if (startRec.isAlleleCoded()) {
            return false;
        }
        else {
            IntArray hapToSeq = ref.get(start).map(0);
            for (int m=start+1; m<end; ++m) {
                RefGTRec rec = ref.get(m);
                if (rec.isAlleleCoded() || rec.map(0)!=hapToSeq) {
                    return false;
                }
            }
            return true;
        }
    }

    private IndexArray codeTarg(int start, int[][] seqMap) {
        int nTargHaps = targ.nHaps();
        int[] hapToSeq = new int[nTargHaps];
        Arrays.fill(hapToSeq, 1);
        int seqCnt = 2;
        for (int j=0; j<seqMap.length; ++j) {
            int m = start + j;
            int nAlleles = ref.marker(m).nAlleles();
            seqMap[j] = new int[seqCnt*nAlleles];
            seqCnt = 1;
            for (int h=0; h<nTargHaps; ++h) {
                int index = nAlleles*hapToSeq[h] + targ.allele(m, h);
                if (seqMap[j][index]==0) {
                    seqMap[j][index] = seqCnt++;
                }
                hapToSeq[h] = seqMap[j][index];
            }
        }
        IntArray intArray = IntArray.create(hapToSeq, seqCnt);
        return new IndexArray(intArray, seqCnt);
    }

    private IndexArray codeSeqCodedRef(int start, int end) {
        int[][] seqMap = new int[end - start][];
        IndexArray codedTarg = codeTarg(start, seqMap);
        RefGTRec rec = ref.get(start);
        assert rec.nMaps()==2;
        IntArray seqToAllele = rec.map(1);
        int[] seq1ToSeq2 = new int[seqToAllele.size()];
        Arrays.fill(seq1ToSeq2, 1);
        for (int j=0; j<seqMap.length; ++j) {
            int m = start + j;
            int nAlleles = ref.marker(m).nAlleles();
            rec = ref.get(m);
            assert rec.nMaps()==2;
            seqToAllele = rec.map(1);
            for (int s=0; s<seq1ToSeq2.length; ++s) {
                if (seq1ToSeq2[s]>0) {
                    int index = seq1ToSeq2[s]*nAlleles + seqToAllele.get(s);
                    seq1ToSeq2[s] = seqMap[j][index];
                }
            }
        }
        return combine(ref.get(start).map(0), seq1ToSeq2, codedTarg);
    }

    private IndexArray combine(IntArray refBasicIndexArray, int[] seq1ToSeq2,
            IndexArray codedTarg) {
        IntArray targSeq = codedTarg.intArray();
        IntArray combined = new IntArray() {

            @Override
            public int size() {
                return nHaps;
            }

            @Override
            public int get(int index) {
                if (index < nRefHaps) {
                    return seq1ToSeq2[refBasicIndexArray.get(index)];
                }
                else {
                    return targSeq.get(index - nRefHaps);
                }
            }
        } ;
        return new IndexArray(combined, codedTarg.valueSize());
    }

    private IndexArray codeSeq(int start, int end) {
        int[][] seqMap = new int[end - start][];
        IndexArray codedTarg = codeTarg(start, seqMap);
        int[] codedRefHap = new int[nRefHaps];
        Arrays.fill(codedRefHap, 1);
        for (int j=0; j<seqMap.length; ++j) {
            int m = start + j;
            int nAlleles = ref.marker(m).nAlleles();
            for (int h=0; h<codedRefHap.length; ++h) {
                if (codedRefHap[h]>0) {
                    int index = codedRefHap[h]*nAlleles + ref.allele(m, h);
                    codedRefHap[h] = seqMap[j][index];
                }
            }
        }
        return combine(codedRefHap, codedTarg);
    }

    private IndexArray combine(int[] codedRef, IndexArray codedTarg) {
        IntArray targSeq = codedTarg.intArray();
        IntArray combined = new IntArray() {

            @Override
            public int size() {
                return nHaps;
            }

            @Override
            public int get(int index) {
                if (index < nRefHaps) {
                    return codedRef[index];
                }
                else {
                    return targSeq.get(index - nRefHaps);
                }
            }
        } ;
        return new IndexArray(combined, codedTarg.valueSize());
    }
}
