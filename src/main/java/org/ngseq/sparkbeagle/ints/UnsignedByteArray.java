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
 * <p>Class {@code ShiftedByteIndexArray} represents an immutable
 * array of integer values between 0 and 255 inclusive that is stored
 * as a {@code byte[]} array whose values have been translated by -128.
 * </p>
 * <p>
 * Instances of {@code UnsignedByteArray} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class UnsignedByteArray implements IntArray {

    private final byte[] ba;

    /**
     * Constructs a new {@code UnsignedByteArray} instance from
     * the specified data.
     * @param ba an array of bytes which are interpreted as unsigned byte
     * values between 0 and 255
     * @throws NullPointerException if {@code ba == null}
     */
    public UnsignedByteArray(byte[] ba) {
        this.ba = ba.clone();
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the
     * specified data.
     * @param ba an array of bytes which are interpreted as unsigned byte
     * values between 0 and 255
     * @param from the first element to be included (inclusive)
     * @param to the last element to be included (exclusive)
     * @throws IndexOutOfBoundsException if {@code (from < 0 || to > ia.length)}
     * @throws NegativeArraySizeException if {@code to > from}
     * @throws NullPointerException if {@code ba == null}
     */
    public UnsignedByteArray(byte[] ba, int from, int to) {
        this.ba = Arrays.copyOfRange(ba, from, to);
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from
     * the specified data.
     * @param ia an array of positive integer values whose lower order byte
     * will be stored
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > 255)} for any index {@code j}
     * satisfying {@code (j >= 0 && j < ia.length)}
     * @throws NullPointerException if {@code ia == null}
     */
    public UnsignedByteArray(int[] ia) {
        this(ia, 0, ia.length);
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the specified data.
     * @param il an list of integer values between 0 and 255 inclusive
     * @throws IllegalArgumentException if
     * {@code (il.get(j) < 0 || il.get(j) > 255)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < il.size())}
     * @throws NullPointerException if {@code il == null}
     */
    public UnsignedByteArray(IntList il) {
        this(il, 0, il.size());
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from
     * the specified data.
     * @param ia an array of nonnegative  integer values
     * @param valueSize the exclusive end of the range of non-negative
     * array values
     * @throws IllegalArgumentException if {@code valueSize < 1 || valueSize > 256}
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < ia.length)}
     * @throws NullPointerException if {@code ia == null}
     */
    public UnsignedByteArray(int[] ia, int valueSize) {
        if (valueSize < 1 || valueSize > 256) {
            throw new IllegalArgumentException(String.valueOf(valueSize));
        }
        this.ba = new byte[ia.length];
        for (int j=0; j<ia.length; ++j) {
            if (ia[j]<0 || ia[j]>=valueSize) {
                throw new IllegalArgumentException(String.valueOf(ia[j]));
            }
            this.ba[j] = (byte) ia[j];
        }
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from
     * the specified data.
     * @param il an list of nonnegative integers
     * @param valueSize the exclusive end of the range of non-negative
     * array values
     * @throws IllegalArgumentException if
     * {@code (valueSize < 1) || (valueSize > 256)}
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < ia.length)}
     * @throws NullPointerException if {@code il == null}
     */
    public UnsignedByteArray(IntList il, int valueSize) {
        if (valueSize < 1 || valueSize > 256) {
            throw new IllegalArgumentException(String.valueOf(valueSize));
        }
        this.ba = new byte[il.size()];
        for (int j=0; j<ba.length; ++j) {
            int value = il.get(j);
            if (value<0 || value>=valueSize) {
                throw new IllegalArgumentException(String.valueOf(value));
            }
            this.ba[j] = (byte) value;
        }
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the
     * specified data.
     * @param ia an array of integer values between 0 and 255 inclusive
     * @param from the first element to be included (inclusive)
     * @param to the last element to be included (exclusive)
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > 255)} for any index {@code j}
     * satisfying {@code (j >= from && j < to)}
     * @throws IndexOutOfBoundsException if {@code (from < 0 || to > ia.length)}
     * @throws NegativeArraySizeException if {@code to > from}
     * @throws NullPointerException if {@code ia == null}
     */
    public UnsignedByteArray(int[] ia, int from, int to) {
        this.ba = new byte[to - from];
        for (int j=from; j<to; ++j) {
            if (ia[j] < 0 || ia[j] > 255) {
                throw new IllegalArgumentException(String.valueOf(ia[j]));
            }
            ba[j - from] = (byte) ia[j];
        }
    }

    /**
     * Constructs a new {@code UnsignedByteArray} instance from the
     * specified data.
     * @param il an list of integer values between 0 and 255 inclusive
     * @param from the first element to be included (inclusive)
     * @param to the last element to be included (exclusive)
     * @throws IllegalArgumentException if
     * {@code (il.get(j) < 0 || il.get(j) > 255)} for any index {@code j}
     * satisfying  {@code (j >= from && j < to)}
     * @throws IndexOutOfBoundsException if {@code from < 0 || to > il.length}
     * @throws NegativeArraySizeException if {@code to > from}
     * @throws NullPointerException if {@code il == null}
     */
    public UnsignedByteArray(IntList il, int from, int to) {
        this.ba = new byte[to - from];
        for (int j=from; j<to; ++j) {
            int value = il.get(j);
            if (value < 0 || value > 255) {
                throw new IllegalArgumentException(String.valueOf(value));
            }
            ba[j - from] = (byte) value;
        }
    }

    @Override
    public int size() {
        return ba.length;
    }

    @Override
    public int get(int index) {
        return ba[index] & 0xff;
    }

    @Override
    public String toString() {
        return IntArray.asString(this);
    }
}
