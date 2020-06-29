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

import java.util.Arrays;

/**
 * Class {@code FloatArray} represents an immutable list of float floating
 * point values.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class FloatArray {

    private final float[] values;

    /**
     * Constructs an {@code FloatArray} object with the specified
     * values.
     * @param values the list of floating point values
     * @throws NullPointerException if {@code values == null}
     */
    public FloatArray(float[] values) {
        this.values = values.clone();
    }

    /**
     * Returns the float at the specified position in this list.
     * @param index the index of the returned float
     * @return the float at the specified position in this list
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= size}
     */
    public float get(int index) {
        return values[index];
    }

    /**
     * Returns the number of elements in this list.
     * @return the number of elements in this list
     */
    public int size() {
        return values.length;
    }

    /**
     * Returns {@code true} if this list has no elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if this list has no elements, and returns
     * {@code false} otherwise
     */
    public boolean isEmpty() {
        return values.length==0;
    }

    /**
     * Returns an integer array containing the sequence of elements in this
     * list.
     * @return an integer array containing the sequence of elements in this
     * list
     */
    public float[] toArray() {
        return values.clone();
    }

    /**
     * Returns a string representation of this list that is
     * obtained by calling {@code java.util.Arrays.toString(this.toArray())}.
     *
     * @return a string representation of this list
     */
    @Override
    public String toString() {
        return Arrays.toString(values);
    }
}
