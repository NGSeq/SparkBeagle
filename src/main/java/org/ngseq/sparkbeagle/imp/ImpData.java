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
import org.ngseq.sparkbeagle.Par;
import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.ints.IndexArray;
import org.ngseq.sparkbeagle.ints.IntArray;
import org.ngseq.sparkbeagle.ints.IntList;
import org.ngseq.sparkbeagle.vcf.*;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code ImpData} contains the input data for imputation of
 * ungenotyped markers.
 * </p>
 * <p>Instances of class {@code ImpData} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ImpData {

    private static final double MIN_CM_DIST = 1e-7;

    private final Par par;
    private final CurrentData cd;
    private final RefGT refGT;
    private final GT phasedTarg;
    private final int[] targClustStartEnd;
    private final int[] refClusterStart;
    private final int[] refClusterEnd;
    private final IndexArray[] hapToSeq;
    private final float[] errProb;
    private final double[] pos;
    private final float[] pRecomb;
    private final float[] weight;
    private final int nClusters;
    private final int nRefHaps;
    private final int nTargHaps;
    private final int nHaps;

    /**
     * Constructs a new {@code ImpData} instance from the specified data.
     * @param par the analysis parameters
     * @param cd the input data for the current marker window
     * @param phasedTarg the phased target genotypes
     * @param map the genetic map
     *
     * @throws IllegalArgumentException if
     * {@code cd.targMarkers().equals(phasedTarg.markers() == false}
     * @throws IllegalArgumentException if
     * {@code cd.targSamples().equals(phasedTarg.samples()) == false}
     * @throws IllegalArgumentException if
     * {@code phasedTarg.isPhased() == false}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public ImpData(Par par, CurrentData cd, GT phasedTarg, GeneticMap map) {
        if (cd.targMarkers().equals(phasedTarg.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (cd.targSamples().equals(phasedTarg.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (phasedTarg.isPhased() == false) {
            throw new IllegalArgumentException("unphased data");
        }
        int[] targToRef = cd.markerIndices();
        this.par = par;
        this.cd = cd;
        this.refGT = cd.refGT();
        this.phasedTarg = phasedTarg;
        double[] targPos = cumPos(phasedTarg.markers(), map);
        int[] blockEnd = targBlockEnd(refGT, targToRef);
        this.targClustStartEnd = targClustStartEnd(targPos, blockEnd, par.cluster());
        this.pos = midPos(targPos, targClustStartEnd);
        this.nClusters = targClustStartEnd.length - 1;
        this.nRefHaps = refGT.nHaps();
        this.nTargHaps = phasedTarg.nHaps();
        this.nHaps = refGT.nHaps() + phasedTarg.nHaps();
        this.hapToSeq = hapToSeq(cd.restrictRefGT(), phasedTarg, targClustStartEnd);
        this.refClusterStart = refClustStart(targClustStartEnd, targToRef);
        this.refClusterEnd = refClustEnd(targClustStartEnd, targToRef);
        this.errProb = err(par.err(), targClustStartEnd);
        this.pRecomb = pRecomb(par.ne(), refGT.nHaps(), pos);
        this.weight = wts(refGT.markers(), refClusterStart, refClusterEnd, map);
    }

    private static double[] cumPos(Markers markers, GeneticMap map) {
        double[] cumPos = new double[markers.nMarkers()];
        double lastGenPos = map.genPos(markers.marker(0));
        cumPos[0] = 0.0;
        for (int j=1; j<cumPos.length; ++j) {
            double genPos = map.genPos(markers.marker(j));
            double genDist = Math.max(Math.abs(genPos - lastGenPos), MIN_CM_DIST);
            cumPos[j] = cumPos[j-1] + genDist;
            lastGenPos = genPos;
        }
        return cumPos;
    }

    private static int[] targBlockEnd(RefGT refGT, int[] targToRef) {
        IntList intList = new IntList(targToRef.length/4);
        IntArray lastHap2Seq = null;
        for (int j=0; j<targToRef.length; ++j) {
            int refIndex = targToRef[j];
            RefGTRec rec = refGT.get(refIndex);
            if (rec.isAlleleCoded()==false) {
                IntArray hap2Seq = rec.map(0);
                if (hap2Seq!=lastHap2Seq) {
                    if (lastHap2Seq!=null) {
                        intList.add(j);
                    }
                    lastHap2Seq = hap2Seq;
                }
            }
        }
        intList.add(targToRef.length);
        return intList.toArray();
    }

    /*
     * indices in int[] targBlockEnd are adjusted when method returns
     */
    private static int[] targClustStartEnd(double[] rawPos, int[] targBlockEnd,
            float clusterDist) {
        int[] clustStartEnd = new int[rawPos.length+1];
        int size = 1;   // clustStartEnd[0] = 0
        for (int j=0; j<targBlockEnd.length; ++j) {
            int clustStart = clustStartEnd[size - 1];
            int blockEnd = targBlockEnd[j];
            double startPos = rawPos[clustStart];
            for (int m=clustStart+1; m<blockEnd; ++m) {
                double pos = rawPos[m];
                if ((pos - startPos) > clusterDist)  {
                    clustStartEnd[size++] = m;
                    startPos = pos;
                }
            }
            clustStartEnd[size++] = blockEnd;
            targBlockEnd[j] = size-2;   // size = nClusters + 1
        }
        return Arrays.copyOf(clustStartEnd, size);
    }

    private static double[] midPos(double[] pos, int[] startEnd) {
        return IntStream.range(1, startEnd.length)
                .mapToDouble(j -> (pos[startEnd[j-1]] + pos[startEnd[j]-1])/2)
                .toArray();
    }

    private static IndexArray[] hapToSeq(RefGT restrictRef, GT phasedTarg,
            int[] targStartEnd) {
        HaplotypeCoder coder = new HaplotypeCoder(restrictRef, phasedTarg);
        return IntStream.range(1, targStartEnd.length)
                .mapToObj(j -> coder.run(targStartEnd[j-1], targStartEnd[j]))
                .toArray(IndexArray[]::new);
    }

    private static float[] err(float errRate, int[] startEnd) {
        float maxErrProb = 0.5f;
        float[] err = new float[startEnd.length - 1];
        for (int j=0; j<err.length; ++j) {
            err[j] = errRate * (startEnd[j+1] - startEnd[j]);
            if (err[j] > maxErrProb) {
                err[j] = maxErrProb;
            }
        }
        return err;
    }

    private static int[] refClustStart(int[] clustStartEnd, int[] targToRef) {
        return IntStream.range(0, clustStartEnd.length-1)
                .map(j -> targToRef[clustStartEnd[j]])
                .toArray();
    }

    private static int[] refClustEnd(int[] clustStartEnd, int[] targToRef) {
        return IntStream.range(1, clustStartEnd.length)
                .map(j -> targToRef[clustStartEnd[j] - 1] + 1)
                .toArray();
    }

    private static float[] pRecomb(float ne, int nHaps, double[] pos) {
        float[] pRecomb = new float[pos.length];
        double c = -(0.04*ne/nHaps);    // 0.04 = 4/(100 cM/M)
        for (int j=1; j<pRecomb.length; ++j) {
            pRecomb[j] = (float) -Math.expm1(c*(pos[j] - pos[j-1]));
        }
        return pRecomb;
    }

    private static float[] wts(Markers refMarkers, int[] refClusterStart,
            int[] refClusterEnd, GeneticMap map) {
        double[] cumPos = cumPos(refMarkers, map);
        int nTargMarkersM1 = refClusterStart.length - 1;
        float[] wts = new float[cumPos.length];
        Arrays.fill(wts, 0, refClusterStart[0], Float.NaN);
        for (int j=0; j<nTargMarkersM1; ++j) {
            int start = refClusterStart[j];
            int end = refClusterEnd[j];
            int nextStart = refClusterStart[j+1];
            double nextStartPos = cumPos[nextStart];
            double totalLength = nextStartPos - cumPos[end - 1];
            Arrays.fill(wts, start, end, Float.NaN);
            for (int m=end; m<nextStart; ++m) {
                wts[m] = (float) ((cumPos[nextStart] - cumPos[m]) / totalLength);
            }
        }
        Arrays.fill(wts, refClusterStart[nTargMarkersM1], refMarkers.nMarkers(),
                Float.NaN);
        return wts;
    }

    /**
     * Returns the command line parameters
     * @return the command line parameters
     */
    public Par par() {
        return par;
    }

    /**
     * Returns the current data.
     * @return the current data
     */
    public CurrentData cd() {
        return cd;
    }

    /**
     * Return the reference genotype data
     * @return the reference genotype data
     */
    public RefGT refGT() {
        return refGT;
    }

    /**
     * Return the phased target genotype data.  The {@code isPhased()} method
     * of the returned object returns {@code true}.
     * @return the phased target genotype data
     */
    public GT targGT() {
        return phasedTarg;
    }

    /**
     * Returns the target marker index corresponding to the start (inclusive)
     * of the specified marker cluster.
     * @param cluster index of a target marker cluster
     * @return the target marker index corresponding to the start (inclusive)
     * of the specified marker cluster
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public int targClusterStart(int cluster) {
        if (cluster >= nClusters) {
            throw new IndexOutOfBoundsException(String.valueOf(cluster));
        }
        return targClustStartEnd[cluster];
    }

    /**
     * Returns the target marker index corresponding to the end (exclusive) of
     * the specified marker cluster.
     * @param cluster index of a target marker cluster
     * @return the target marker index corresponding to the end (exclusive)
     * of the specified marker cluster
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public int targClusterEnd(int cluster) {
        if (cluster < 0) {
            throw new IndexOutOfBoundsException(String.valueOf(cluster));
        }
        return targClustStartEnd[cluster + 1];
    }

    /**
     * Returns the index of the reference marker corresponding to the start
     * (inclusive) of the specified target marker cluster.
     * @param cluster index of a target marker cluster
     * @return the index of the reference marker corresponding to the start
     * (inclusive) of the specified target marker cluster
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public int refClusterStart(int cluster) {
        return refClusterStart[cluster];
    }

    /**
     * Returns the index of the reference marker corresponding to the end
     * (exclusive) of the specified target marker cluster.
     * @param cluster index of a target marker cluster
     * @return the index of the reference marker corresponding to the end
     * (exclusive) of the specified target marker cluster
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public int refClusterEnd(int cluster) {
        return refClusterEnd[cluster];
    }

    /**
     * Return the number of target marker clusters.
     * @return the number of target marker clusters
     */
    public int nClusters() {
        return nClusters;
    }

    /**
     * Returns the list of target samples.
     * @return the list of target samples
     */
    public Samples targSamples() {
        return phasedTarg.samples();
    }

    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */
    public int nTargSamples() {
        return phasedTarg.nSamples();
    }

    /**
     * Return the total number of reference and target haplotypes.
     * @return the total number of reference and target haplotypes
     */
    public int nHaps() {
        return nHaps;
    }

   /**
     * Return the number of reference haplotypes.
     * @return the number of reference haplotypes
     */
    public int nRefHaps() {
        return nRefHaps;
    }

    /**
     * Return the number of target haplotypes.
     * @return the number of target haplotypes
     */
    public int nTargHaps() {
        return nTargHaps;
    }

    /**
     * Returns the specified target marker cluter alleles for the
     * reference and target haplotypes.  Alleles for the
     * reference haplotypes precede alleles for the target haplotypes.  If
     * {@code (this.nRefHaps() <= hap && hap < this.nHaps())} then
     * {@code (this.allele(marker, hap) ==
     * this.targAllele(marker, hap - this.nRefHaps())}
     * @param cluster index of a target marker cluster
     * @param hap a haplotype index
     * @return the specified target marker cluster allele for the specified
     * haplotype
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     * @throws IndexOutOfBoundsException if
     * {@code haplotype < 0 || haplotype >= this.nHaps()}
     */
    public int allele(int cluster, int hap) {
        return hapToSeq[cluster].get(hap);
    }

    /**
     * Returns the specified target marker cluster alleles for the
     * reference and target haplotypes.  Alleles for the
     * reference haplotypes precede alleles for the target haplotypes. The
     * returned value will satisfy
     * {@code (this.hapToSeq(cluster).get(hap)==this.allele(cluster, hap))}
     * for any {@code cluster} and {@code hap} satisfying
     * {@code (0 <= cluster && cluster < this.nClusters())} and
     * {@code (0 <= hap && hap < this.nHaps())}
     * @param cluster index of a target marker cluster
     * @return the specified target marker cluster alleles for the
     * reference and target haplotypes
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public IndexArray hapToSeq(int cluster) {
        return hapToSeq[cluster];
    }

    /**
     * Returns the probability that the allele carried by the specified
     * target marker cluster matches the allele labeling the latent HMM state.
     * @param cluster index of a target marker cluster
     * @return the probability that the allele carried by the specified
     * target marker cluster matches the allele labeling the latent HMM state.
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public float errProb(int cluster) {
        return errProb[cluster];
    }

    /**
     * Return the genetic map position of the specified target marker cluster.
     * @param cluster index of a target marker cluster
     * @return the genetic map position of the specified target marker cluster
     * @throws IllegalArgumentException if
     * {@code cluster < 0 || marker >= this.nClusters()}
     */
    public double pos(int cluster) {
        return pos[cluster];
    }

    /**
     * Return an array of size {@code this.nClusters()} containing the
     * the genetic map positions of the target marker clusters.
     * @return the genetic map positions of the target marker clusters
     */
    public double[] pos() {
        return pos.clone();
    }

    /**
     * Return the probability of recombination between the specified
     * target marker cluster and the previous target marker cluster.
     * Returns {@code 0} if {@code (cluster == 0)}.
     * @param cluster index of a target marker cluster
     * @return the probability of recombination between the specified
     * target marker cluster and the previous target marker cluster
     * @throws IllegalArgumentException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public float pRecomb(int cluster) {
        return pRecomb[cluster];
    }

    /**
     * Return the weight for the HMM state probability at the
     * preceding target marker cluster when estimating the HMM state
     * probability at the specified reference marker via linear interpolation
     * of HMM state probabilities at the preceding and succeeding target
     * marker clusters.
     *
     * @param refMarker a reference marker index
     * @return the weight for the HMM state probability at the preceding
     * target marker cluster when estimating the HMM state
     * probability at the specified reference marker via linear interpolation
     * @throws IllegalArgumentException if
     * {@code refMarker < 0 || refMarker >= this.refGT().nMarkers()}
     */
    public double weight(int refMarker) {
        return weight[refMarker];
    }
}
