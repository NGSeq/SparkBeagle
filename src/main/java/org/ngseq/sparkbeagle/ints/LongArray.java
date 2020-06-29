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
package org.ngseq.sparkbeagle.ints;

import java.util.Arrays;

/**
 * <p>Interface {@code LongArray} represents an immutable {@code long[]} array.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class LongArray {

    private final long[] values;

    /**
     * Constructs a {@code LongArray} instance from the specified values.
     * @param values a long array
     * @throws NullPointerException if {@code values == null}
     */
    public LongArray(long[] values) {
        this.values = values.clone();
    }

    /**
     * Returns the number of elements in this {@code LongArray}.
     * @return the number of elements in this {@code LongArray}
     */
    public int size() {
        return values.length;
    }

    /**
     * Returns the specified array element.
     * @param index an array index
     * @return the specified array element
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    public long get(int index) {
        return values[index];
    }

    /**
     * Returns a copy of the specified array.
     * @param la a list of longs
     * @return a copy of the specified array
     * @throws NullPointerException if {@code ia == null}
     */
    public static long[] toArray(LongArray la) {
        long[] copy = new long[la.size()];
        for (int j=0; j<copy.length; ++j) {
            copy[j] = la.get(j);
        }
        return copy;
    }

    /**
     * Returns a string representation of this {@code LongArray} by applying
     * {@code java.utils.Arrays.toString()} to an equivalent {@code int[]}
     * object.
     *
     * @param ia a list of longs
     * @return a string representation of this {@code LongArray}.
     * @throws NullPointerException if {@code ia == null}
     */
    public static String asString(LongArray ia) {
        return Arrays.toString(toArray(ia));
    }

    /**
     * Returns {@code true} if the specified {@code LongArray} objects
     * represent the same sequence of long values, and returns {@code false}
     * otherwise.
     * @param a a sequence of long values
     * @param b a sequence of long values
     * @return {@code true} if the specified {@code LongArray} objects
     * represent the same sequence of long values
     */
    public static boolean equals(LongArray a, LongArray b) {
        return Arrays.equals(a.values, b.values);
    }
}
