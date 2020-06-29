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

import java.util.Arrays;

/**
 * <p>Class {@code BasicRefGT} stores a list of samples and a
 * haplotype pair for each sample.
 * </p>
 * <p>Instances of class {@code BasicRefGT} are immutable.<p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class BasicRefGT implements RefGT {

    private final Markers markers;
    private final Samples samples;
    private final RefGTRec[] recs;

    /**
     * Constructs a new {@code RefHapPairs} instance.
     * @param markers the sequence of markers
     * @param samples the sequence of samples
     * @param refVcfRecs the sequence of per-marker genotype data
     *
     * @throws IllegalArgumentException if
     * {@code markers.nMarkers() != refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].samples().equals(samples) == false} for any
     * {@code k} satisfying {@code 0 <= k  && k < refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].marker().equals(markers.marker(k)) == false}
     * for any {@code k} satisfying {@code 0 <= k && k < refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].isPhased() == false} for any {@code k}
     * satisfying {@code 0 <= k && k < refVcfRecs.length}
     * @throws NullPointerException if
     * {@code markers == null || samples == null || refVcfRecs == null
     * || refVcfRecs[k] == null} for any {@code k} satisfying
     * {@code 0 <= k && k <= refVcfRecs.length}
     */
    public BasicRefGT(Markers markers, Samples samples,
            RefGTRec[] refVcfRecs) {
        checkData(markers, samples, refVcfRecs);
        this.markers = markers;
        this.samples = samples;
        this.recs = refVcfRecs.clone();
    }

    /**
     * Constructs a new {@code RefHapPairs} instance.
     * @param refVcfRecs the sequence of per-marker genotype data
     *
     * @throws IllegalArgumentException if {@code refVcfRecs.length == 0}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].samples().equals(samples) == false} for any
     * {@code k} satisfying {@code 0 <= k  && k < refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].isPhased() == false} for any {@code k}
     * satisfying {@code 0 <= k && k < refVcfRecs.length}
     * @throws NullPointerException if
     * {@code samples == null || refVcfRecs == null}
     * @throws NullPointerException if
     * {@code (refVcfRecs[k] == null)} for any {@code k} satisfying
     * {@code (0 <= k && k <= refVcfRecs.length)}
     */
    public BasicRefGT(RefGTRec[] refVcfRecs) {
        this.samples = checkData(refVcfRecs);
        Marker[] ma = Arrays.stream(refVcfRecs)
                .parallel()
                .map(rec -> rec.marker())
                .toArray(Marker[]::new);
        this.markers = Markers.create(ma);
        this.recs = refVcfRecs.clone();
    }

    private static Samples checkData(GTRec[] refVcfRecs) {
        if (refVcfRecs.length==0) {
            String s = "Missing data in VCF file";
            throw new IllegalArgumentException(s);
        }
        Samples samples = refVcfRecs[0].samples();
        for (int j=0; j<refVcfRecs.length; ++j) {
            if (refVcfRecs[j].samples().equals(samples)==false) {
                String s = "sample inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (refVcfRecs[j].isPhased()==false) {
                String s = "non-reference data at marker index " + j;
                throw new IllegalArgumentException(s);
            }
        }
        return samples;
    }

    private static void checkData(Markers markers, Samples samples,
            GTRec[] refVcfRecs) {
        if (markers.nMarkers()!=refVcfRecs.length) {
            String s = "markers.nMarkers()=" + markers.nMarkers()
                    + " refVcfRecs.length=" + refVcfRecs.length;
            throw new IllegalArgumentException(s);
        }
        for (int j=0; j<refVcfRecs.length; ++j) {
            if (refVcfRecs[j].samples().equals(samples)==false) {
                String s = "sample inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (refVcfRecs[j].marker().equals(markers.marker(j))==false) {
                String s = "marker inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (refVcfRecs[j].isPhased()==false) {
                String s = "non-reference data at marker index " + j;
                throw new IllegalArgumentException(s);
            }
        }
    }

    @Override
    public boolean isPhased(int marker, int sample) {
        if (marker < 0 || marker >= this.nMarkers()) {
            throw new IndexOutOfBoundsException(String.valueOf(marker));
        }
        if (sample < 0 || sample >= this.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        return true;
    }

    @Override
    public boolean isPhased(int sample) {
        if (sample < 0 || sample >= this.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        return true;
    }

    @Override
    public boolean isPhased() {
        return true;
    }

    @Override
    public int allele1(int marker, int hapPair) {
        return recs[marker].allele1(hapPair);
    }

    @Override
    public int allele2(int marker, int hapPair) {
        return recs[marker].allele2(hapPair);
    }

    @Override
    public int allele(int marker, int haplotype) {
        return recs[marker].get(haplotype);
    }

    @Override
    public int nMarkers() {
       return markers.nMarkers();
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public Marker marker(int marker) {
        return markers.marker(marker);
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
    public RefGTRec get(int marker) {
        return recs[marker];
    }
}
