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

import org.ngseq.sparkbeagle.blbutil.Const;
import org.ngseq.sparkbeagle.vcf.Marker;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.stream.IntStream;

/**
 * <p>Class {@code ImputeRecBuilder} contains methods for constructing
 * and printing a VCF record in VCF 4.3 format.  The sample data in
 * the output VCF record are in the same order that the data were added
 * with the {@code addSampleData()} method.
 * </p>
 * <p>Instances of class {@code ImputeRecBuilder} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ImputedRecBuilder {

    private static final DecimalFormat DF = new DecimalFormat("#.##");
    private static final DecimalFormat DF2 = new DecimalFormat("0.00");
    private static final DecimalFormat DF4 = new DecimalFormat("0.0000");

    private static final String[] DS_VALS = IntStream.range(0, 201)
            .mapToObj(j -> DF.format(j/100.0))
            .toArray(String[]::new);
    private static final String[] R2_VALS = IntStream.range(0, 101)
            .limit(101)
            .mapToObj(i -> (DS_VALS[i].length()!=4 ? DF2.format(i/100.0) : DS_VALS[i]))
            .toArray(String[]::new);

    private static final String[] DEFAULT_HOM_REF_FIELDS = defaultHomRefFields();

    private final Marker marker;
    private final int nAlleles;
    private final int nHaps;
    private final boolean ap;
    private final boolean gp;
    private final float[] sumAlProbs;
    private final float[] sumAlProbs2;
    private final String[] homRefField;
    private StringBuilder sampleData = null;

    private int sampleCnt;

    /**
     * Constructs a new {@code ImputedRecBuilder} instance for the specified
     * number of samples.
     *
     * @param marker the marker corresponding to the VCF record
     * @param nSamples the number of samples
     * @param ap {@code true} if posterior allele probabilities are to be printed
     * @param gp {@code true} if posterior genotype probabilities are to be printed
     * @throws IllegalArgumentException if {@code nSamples < 1}
     * @throws NullPointerException if {@code marker == null}
     */
    public ImputedRecBuilder(Marker marker, int nSamples, boolean ap,
            boolean gp) {
        if (nSamples < 1) {
            throw new IllegalArgumentException(String.valueOf(nSamples));
        }
        this.marker = marker;
        this.nAlleles = marker.nAlleles();
        this.nHaps = 2*nSamples;
        this.ap = ap;
        this.gp = gp;
        this.sumAlProbs = new float[nAlleles];
        this.sumAlProbs2 = new float[nAlleles];
        try{
          this.sampleData = new StringBuilder(nSamples*4);
        }catch(OutOfMemoryError e){
            System.out.println("OUTOFMEM AT ALLELE!!!!!!!!!!!!!! "+marker.id());
            e.printStackTrace();
        }
        this.homRefField = (ap || gp) ? homRefFields(ap, gp)
                : DEFAULT_HOM_REF_FIELDS;
        this.sampleCnt = 0;
    }

    /**
     * Returns the marker in the VCF record.
     * @return the marker in the VCF record
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns the required number of samples in the VCF record.
     * @return the required number of samples in the VCF record
     */
    public int nSamples() {
        return nHaps/2;
    }

    /**
     * Returns the number of times the {@code this.addSampleData()}
     * method has been invoked.
     * @return the number of times the {@code this.addSampleData()}
     * method has been invoked
     */
    public int sampleCnt() {
        return sampleCnt;
    }

    /**
     * Scales the specified probabilities for each allele to each sum to 1.0,
     * and adds the sample data to the VCF record.  The contract
     * for this method is undefined if any element of the specified arrays is
     * not a finite non-negative number.
     * @param a1 the allele probabilities for the first allele
     * @param a2 the allele probabilities for the second allele
     * @throws IndexOutOfBoundsException if
     * {@code a1.length < this.marker().nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code a2.length < this.marker().nAlleles()}
     * @throws NullPointerException if {@code a1 == null || a2 == null}
     */
    public void addSampleData(float[] a1, float[] a2) {
        ++sampleCnt;
        if (a1[0]==1.0f && a2[0]==1.0f && a1.length < DEFAULT_HOM_REF_FIELDS.length) {
            sampleData.append(homRefField[a1.length]);
        }
        else {
            scale(a1);
            scale(a2);
            sampleData.append(Const.tab);
            sampleData.append(maxIndex(a1));
            sampleData.append(Const.phasedSep);
            sampleData.append(maxIndex(a2));
            for (int a=1; a<nAlleles; ++a) {
                float dose = a1[a] + a2[a];
                float dose2 = a1[a]*a1[a] + a2[a]*a2[a];
                sumAlProbs[a] += dose;
                sumAlProbs2[a] += dose2;
                sampleData.append( (a==1) ? Const.colon : Const.comma );
                sampleData.append(DS_VALS[(int) Math.rint(100*dose)]);
            }
            if (ap) {
                for (int a=1; a<nAlleles; ++a) {
                    sampleData.append( (a==1) ? Const.colon : Const.comma );
                    sampleData.append(DS_VALS[(int) Math.rint(100*a1[a])]);
                }
                for (int a=1; a<nAlleles; ++a) {
                    sampleData.append( (a==1) ? Const.colon : Const.comma );
                    sampleData.append(DS_VALS[(int) Math.rint(100*a2[a])]);
                }
            }
            if (gp) {
                for (int i2=0; i2<nAlleles; ++i2) {
                    for (int i1=0; i1<=i2; ++i1) {
                        float prob = a1[i1]*a2[i2];
                        if (i1!=i2) {
                            prob += a1[i2]*a2[i1];
                        }
                        sampleData.append( (i2==0) ? Const.colon : Const.comma );
                        sampleData.append(DS_VALS[(int) Math.rint(100*prob)]);
                    }
                }
            }
        }
    }

    private static void scale(float[] fa) {
        float sum = 0f;
        for (float f : fa) {
            sum += f;
        }
        for (int j=0; j<fa.length; ++j) {
            fa[j] /= sum;
        }
    }

    private static int maxIndex(float[] fa) {
        int maxIndex=0;
        for (int j=1; j<fa.length; ++j) {
            if (fa[j] > fa[maxIndex]) {
                maxIndex = j;
            }
        }
        return maxIndex;
    }

    /**
     * Prints the VCF record to the specified {@code PrintWriter}.
     * The INFO field of the VCF record will include the DR2 (dose r2) and
     * AF (ALT allele frequency) subfields.
     * @param isImputed {@code true} if the INFO field of the VCF record will
     * have an IMP flag and {@code false} otherwise
      *@param out the {@code PrintWriter} to which the VCF record will be
     * printed
     * @throws IllegalStateException if
     * {@code this.nSamples() != this.sampleCnt()}
     * @throws NullPointerException if {@code out == null}
     */
    public void printRec(PrintWriter out, boolean isImputed) {
        if (2*sampleCnt != nHaps) {
            throw new IllegalStateException("inconsistent data");
        }
        printMarkerFields(marker, out);
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR);             // QUAL
        out.print(Const.tab);
        out.print("PASS");                              // FILTER
        out.print(Const.tab);
        printInfoField(out, isImputed);                 // INFO
        out.print(Const.tab);
        out.print("GT:DS");                             // FORMAT
        if (ap) {
            out.print(":AP1:AP2");
        }
        if (gp) {
            out.print(":GP");
        }
        out.println(sampleData);
    }

    private static void printMarkerFields(Marker marker, PrintWriter out) {
        out.print(marker.chrom());
        out.print(Const.tab);
        out.print(marker.pos());
        int nIds = marker.nIds();
        if (nIds==0) {
            out.print(Const.tab);
            out.print(Const.MISSING_DATA_CHAR);
        }
        else {
            for (int j=0; j<nIds; ++j) {
                out.print(j==0 ? Const.tab : Const.semicolon);
                out.print(marker.id(j));
            }
        }
        int nAlleles = marker.nAlleles();
        if (nAlleles==1) {
            out.print(Const.tab);
            out.print(marker.allele(0));
            out.print(Const.tab);
            out.print(Const.MISSING_DATA_CHAR);
        }
        else {
            for (int j=0; j<nAlleles; ++j) {
                out.print(j<2 ? Const.tab : Const.comma);
                out.print(marker.allele(j));
            }
        }
    }

    private void printInfoField(PrintWriter out, boolean isImputed) {
        if (nAlleles==1) {
            if (isImputed) {
                out.print("IMP");
            }
        }
        else {
            for (int a=1; a<nAlleles; ++a) {
                out.print( (a==1) ? "DR2=" : Const.comma);
                out.print(R2_VALS[(int) Math.rint(100*r2(a))]);
            }
            for (int a=1; a<nAlleles; ++a) {
                out.print( (a==1) ? ";AF=" : Const.comma);
                out.print(DF4.format(sumAlProbs[a]/nHaps));
            }
            if (marker.end()!=-1) {
                out.print(";END=");
                out.print(marker.end());
            }
            if (isImputed) {
                out.print(";IMP");
            }
        }
    }

    private float r2(int allele) {
        float sum = sumAlProbs[allele];
        if (sum==0f) {
            return 0f;
        }
        else {
            float sum2 = sumAlProbs2[allele];
            float meanTerm = sum*sum/(nHaps);
            float num = (sum2 - meanTerm);
            float den = (sum - meanTerm);
            return num <= 0 ? 0f : num/den;
        }
    }

    private static String[] defaultHomRefFields() {
        String[] sa = new String[5];
        sa[1] = Const.tab + "0|0";
        sa[2] = Const.tab + "0|0:0";
        for (int j=3; j<sa.length; ++j) {
            sa[j] = sa[j-1] + ",0";
        }
        return sa;
    }

    private static String[] homRefFields(boolean ap, boolean gp) {
        String[] homRefField = new String[DEFAULT_HOM_REF_FIELDS.length];
        for (int nAl=1; nAl<homRefField.length; ++nAl) {
            StringBuilder sb = new StringBuilder(DEFAULT_HOM_REF_FIELDS[nAl]);
            if (ap) {
                for (int a=1; a<nAl; ++a) {
                    sb.append( (a==1) ? Const.colon : Const.comma );
                    sb.append(DS_VALS[0]);
                }
                for (int a=1; a<nAl; ++a) {
                    sb.append( (a==1) ? Const.colon : Const.comma );
                    sb.append(DS_VALS[0]);
                }
            }
            if (gp) {
                sb.append(Const.colon);
                sb.append(DS_VALS[100]);
                for (int i2=1; i2<nAl; ++i2) {
                    for (int i1=0; i1<=i2; ++i1) {
                        sb.append(Const.comma);
                        sb.append(DS_VALS[0]);
                    }
                }
            }
            homRefField[nAl] = sb.toString();
        }
        return homRefField;
    }
}
