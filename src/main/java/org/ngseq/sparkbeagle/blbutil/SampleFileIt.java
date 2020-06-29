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
package org.ngseq.sparkbeagle.blbutil;

import org.ngseq.sparkbeagle.beagleutil.Samples;

/**
 * <p>An iterator for records in a file.  Each records contains
 * data for the same set of samples.
 *</p>
 * Instances of class {@code SampleFileIt} are not thread-safe.
 *
 * @param <E> the type of the elements returned by this iterator's
 * {@code next()} method.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface SampleFileIt<E> extends FileIt<E> {

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();
}
