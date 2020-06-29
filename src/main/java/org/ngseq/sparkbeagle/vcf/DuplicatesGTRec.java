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

import org.ngseq.sparkbeagle.ints.IntArray;

/**
 * <p>Interface {@code DuplicatesGTRec} represents marker alleles for a
 * list of samples.  The samples in the list of samples are not
 * required to be unique.
 * </p>
 * All instances of {@code HapsMarkers} are required to be
 * immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface DuplicatesGTRec extends MarkerContainer, IntArray {

    /**
     * Returns the first allele for the specified sample or
     * -1 if the allele is missing.  The two alleles for a sample
     * are arbitrarily ordered if
     * {@code this.unphased(marker, sample) == false}.
     * @param sample a sample index
     * @return the first allele for the specified sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    int allele1(int sample);

    /**
     * Returns the second allele for the specified sample or
     * -1 if the allele is missing.  The two alleles for a sample
     * are arbitrarily ordered if
     * {@code this.unphased(marker, sample) == false}.
     * @param sample a sample index
     * @return the second allele for the specified sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    int allele2(int sample);

    /**
     * Returns the specified allele for the specified haplotype or
     * -1 if the allele is missing.  The two alleles for a sample
     * at a marker are arbitrarily ordered if
     * {@code this.unphased(marker, hap/2) == false}.
     * @param hap a haplotype index
     * @return the specified allele for the specified sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.size()}
     */
    @Override
    int get(int hap);

    /**
     * Returns an array of length {@code this.size()} whose {@code j}-th
     * element is equal to {@code this.allele(j}}
     * @return an array of length {@code this.size()} whose {@code j}-th
     * element is equal to {@code this.allele(j}}
     */
    int[] alleles();

    /**
     * Returns the number of haplotypes.  The returned value is equal to
     * {@code 2*this.nSamples()}.
     * @return the number of haplotypes
     */
    @Override
    int size();

    /**
     * Returns the number of samples.  The returned value is
     * equal to {@code this.size()/2}.
     * @return the number of samples
     */
    int nSamples();

    /**
     * Returns {@code true} if the genotype for the specified sample is
     * a phased, nonmissing genotype, and returns {@code false} otherwise.
     * @param sample a sample index
     * @return {@code true} if the genotype for the specified sample
     * is a phased, nonmissing genotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    boolean isPhased(int sample);

    /**
     * Returns {@code true} if every genotype for each sample is a phased,
     * non-missing genotype, and returns {@code false} otherwise.
     * @return {@code true} if the genotype for each sample is a phased,
     * non-missing genotype
     */
    boolean isPhased();
}
