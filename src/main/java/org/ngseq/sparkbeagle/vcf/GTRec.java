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

/**
 * <p>Interface {@code GTRec} represents represents genotype data for one
 * marker.
 * </p>
 * <p>All instances of {@code GTRec} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface GTRec extends DuplicatesGTRec {

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();

    /**
     * Returns {@code true} if the value returned by {@code this.gl()} is
     * determined by a called or missing genotype, and returns {@code false}
     * otherwise.
     * @return {@code true} if the value returned by {@code this.gl()} is
     * determined by a called or missing genotype
     *
     * @implSpec The default implementation returns {@code true}
     */
    boolean isGTData();

    /**
     * Returns the probability of the observed data for the specified sample
     * if the specified pair of ordered alleles is the true ordered genotype.
     * @param sample the sample index
     * @param allele1 the first allele index
     * @param allele2 the second allele index
     * @return the probability of the observed data for the specified sample
     * if the specified pair of ordered alleles is the true ordered genotype.
     *
     * @throws IndexOutOfBoundsException if
     * {@code samples < 0 || samples >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele1 < 0 || allele1 >= this.marker().nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code allele2 < 0 || allele2 >= this.marker().nAlleles()}
     *
     * @implSpec The default implementation returns {@code 1.0f} if the
     * corresponding genotype determined by the {@code isPhased()},
     * {@code allele1()}, and {@code allele2()} methods is consistent
     * with the specified ordered genotype, and returns {@code 0.0f} otherwise.
     */
    float gl(int sample, int allele1, int allele2);

    /**
     * Returns a VCF record corresponding to the specified {@code GTRec} object.
     * The returned VCF record will have missing QUAL and INFO fields,
     * will have "PASS" in the filter field, and will have a GT format field.
     * @param gtRec the genotype data
     * @return a VCF record corresponding to the specified {@code GTRec} object
     * @throws NullPointerException if {@code gtRec == null}
     */
    static String toVcfRec(GTRec gtRec) {
        StringBuilder sb = new StringBuilder(100);
        sb.append(gtRec.marker());
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);  // QUAL
        sb.append(Const.tab);
        sb.append("PASS");                   // FILTER
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);  // INFO
        sb.append(Const.tab);
        sb.append("GT");                     // FORMAT
        for (int j=0, n=gtRec.nSamples(); j<n; ++j) {
            int a1 = gtRec.allele1(j);
            int a2 = gtRec.allele2(j);
            sb.append(Const.tab);
            if (a1==-1) {
                sb.append(Const.MISSING_DATA_CHAR);
            }
            else {
                sb.append(a1);
            }
            sb.append(gtRec.isPhased(j) ? Const.phasedSep : Const.unphasedSep);
            if (a2==-1) {
                sb.append(Const.MISSING_DATA_CHAR);
            }
            else {
                sb.append(a2);
            }
        }
        return sb.toString();
    }
}