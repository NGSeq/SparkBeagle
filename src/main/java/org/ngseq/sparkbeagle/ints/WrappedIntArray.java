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
 * <p>Class {@code WrappedIntArray} represents an immutable
 * {@code int[]} array.
 * </p>
 * Instances of {@code WrappedIntArray} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class WrappedIntArray implements IntArray {

    private final int[] ia;

    /**
     * Constructs a new {@code WrappedIntArray} instance.
     * @param ia an array of integers
     * @throws NullPointerException if {@code ia == null}
     */
    public WrappedIntArray(int[] ia) {
        this.ia = ia.clone();
    }

    /**
     * Constructs a new {@code WrappedIntArray} instance.
     * @param ia an array of integers
     * @param valueSize the exclusive end of the range of non-negative
     * array values
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < ia.length)}
     * @throws NullPointerException if {@code ia == null}
     */
    public WrappedIntArray(int[] ia, int valueSize) {
        this.ia = new int[ia.length];
        for (int j=0; j<ia.length; ++j) {
            if (ia[j]<0 || ia[j]>=valueSize) {
                throw new IllegalArgumentException(String.valueOf(ia[j]));
            }
            this.ia[j] = ia[j];
        }
    }

    /**
     * Constructs a new {@code WrappedIntArray} instance.
     * @param il a list of integers
     * @throws NullPointerException if {@code il == null}
     */
    public WrappedIntArray(IntList il) {
        this.ia = il.toArray();
    }

    /**
     * Constructs a new {@code WrappedIntArray} instance.
     * @param il a list of integers
     * @param valueSize the exclusive end of the range of non-negative
     * array values
     * @throws IllegalArgumentException if
     * {@code (il[j] < 0 || il[j] > valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < il.length)}
     * @throws NullPointerException if {@code il == null}
     */
    public WrappedIntArray(IntList il, int valueSize) {
        this.ia = new int[il.size()];
        for (int j=0; j<ia.length; ++j) {
            int value = il.get(j);
            if (value<0 || value>=valueSize) {
                throw new IllegalArgumentException(String.valueOf(value));
            }
            ia[j] = value;
        }
    }

    @Override
    public int size() {
        return ia.length;
    }

    @Override
    public int get(int index) {
        return ia[index];
    }

    @Override
    public String toString() {
        return Arrays.toString(ia);
    }
}
