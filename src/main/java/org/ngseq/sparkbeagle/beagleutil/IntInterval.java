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
package org.ngseq.sparkbeagle.beagleutil;

import java.util.Comparator;

/**
 * <p>Interface {@code IntInterval} represents an interval of
 * consecutive integers.
 * </p>
 * Instances of class {@code IntInterval} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface IntInterval {

    /**
     * Returns the start of the interval (inclusive).
     * @return the start of the interval (inclusive).
     */
    public int start();

    /**
     * Returns the end of the interval (inclusive).
     * @return the end of the interval (inclusive).
     */
    public int end();

    /**
     * Returns a {@code Comparator<IntInterval>} which orders
     * {@code IntInterval} objects in order of increasing {@code this.start()}
     * value and orders {@code IntInterval} objects with the same
     * {@code this.start()} value in order of increasing {@code this.end()}
     * value.
     * @return a {@code Comparator<IntInterval>} object
     */
    public static Comparator<IntInterval> incEndComp() {
        return (IntInterval t1, IntInterval t2) -> {
            if (t1.start() != t2.start()) {
                return (t1.start() < t2.start()) ? -1 : 1;
            }
            else if (t1.end() != t2.end()) {
                return (t1.end() < t2.end()) ? -1 : 1;
            }
            return 0;
        } ;
    }

    /**
     * Returns a {@code Comparator<IntInterval>} which orders
     * {@code IntInterval} objects in order of increasing {@code this.start()}
     * value and orders {@code IntInterval} objects with the same
     * {@code this.start()} value in order of decreasing {@code this.end()}
     * value.
     * @return a {@code Comparator<IntInterval>} object
     */
    public static Comparator<IntInterval> decEndComp() {
        return (IntInterval t1, IntInterval t2) -> {
            if (t1.start() != t2.start()) {
                return (t1.start() < t2.start()) ? -1 : 1;
            }
            else if (t1.end() != t2.end()) {
                return (t1.end() > t2.end()) ? -1 : 1;
            }
            return 0;
        } ;
    }
}
