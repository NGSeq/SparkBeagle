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
package org.ngseq.sparkbeagle.bref;

import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.vcf.RefGTRec;

import java.io.Closeable;

/**
 * <p>Interface {@code BrefWrites} writes phased, non-missing genotypes to a
 * binary reference format (bref) file.  The {@code close()} method must
 * be called after the last invocation of the {@code write()} method
 * in order to ensure that any buffered data are written to the output
 * binary reference file.
 * </p>
 * <p>Instances of class {@code BrefWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface BrefWriter extends Closeable {

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();

    /**
     * Writes the specified phased genotype data in binary reference format.
     * The Java virtual machine will exit with an error message if an I/O
     * error occurs during method execution, if {@code this.close()}
     * has previously been invoked, or if
     * {@code rec.samples().equals(this.samples()) == false}.
     *
     * @param rec phased genotype data
     *
     * @throws NullPointerException if {@code rec == null}
     */
    void write(RefGTRec rec);

    /**
     * Flushes any buffered output and releases any system resources that are
     * held by this {@code BrefWriter}.  The Java virtual machine will exit
     * with an error message if an I/O error occurs during method execution.
     */
    @Override
    void close();
}
