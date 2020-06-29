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
import java.util.stream.IntStream;

/**
 * <p>Class {@code IntList} represents a list of integers.
 * Class {@code IntList} supports a {@code clear()} method, but it does not
 * support a {@code remove()} method.
 * </p>
 * Class {@code IntList} is not thread-safe.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IntList {

    /**
     * The default initial capacity of an {@code IntList}, which is 16.
     */
    public static final int DEFAULT_INIT_CAPACITY = 16;

    private int size;
    private int[] values;

    /**
     * Constructs an {@code IntList} object with the default
     * initial capacity.
     *
     * @see #DEFAULT_INIT_CAPACITY
     */
    public IntList() {
        this(DEFAULT_INIT_CAPACITY);
    }

    /**
     * Constructs an {@code IntList} object with the specified
     * initial capacity.
     *
     * @param initCapacity the initial capacity of this list
     * @throws IllegalArgumentException if {@code initCapacity < 0}
     */
    public IntList(int initCapacity) {
        if (initCapacity < 0) {
            throw new IllegalArgumentException(String.valueOf(initCapacity));
        }
        this.size = 0;
        this.values = new int[initCapacity];
    }

    /**
     * Constructs an {@code IntList} by cloning the specified array.
     *
     * @param ia a list of integer values
     * @throws NullPointerException if {@code ia == null}
     */
    public IntList(int[] ia) {
        this.size = ia.length;
        this.values = ia.clone();
    }

    /**
     * Constructs an {@code IntList} by copying the specified {@code IntList}.
     *
     * @param intList a list of integer values
     * @throws NullPointerException if {@code intList == null}
     */
    public IntList(IntList intList) {
        this.size = intList.size();
        this.values = Arrays.copyOf(intList.values, intList.size());
    }

    /**
     * Adds the specified integer to the end of this list.
     *
     * @param value the integer to be added to the end of this list
     */
    public void add(int value) {
        if (size==values.length) {
            int newCapacity = (values.length * 3)/2 + 1;
            this.values = Arrays.copyOf(this.values, newCapacity);
        }
        this.values[size++] = value;
    }

    /**
     * Removes and returns the last entry of this list.
     * @return the last entry of this list
     * @throws IndexOutOfBoundsException if {@code this.isEmpty() == true}
     */
    public int pop() {
        return this.values[--size];
    }

    /**
     * Returns the element at the specified position in this list.
     * @param index the index of the element to be returned
     * @return the element at the specified position in this list
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int get(int index) {
        if (index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return values[index];
    }

    /**
     * Replaces the element at the specified position in this list with the
     * specified element.
     * @param index the index of the element to be replaced
     * @param value the value to be stored at the specified position
     * in this list
     * @return the previous value at the specified position in this list
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int set(int index, int value) {
        if (index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int oldValue = values[index];
        values[index] = value;
        return oldValue;
    }

    /**
     * Increments by one the element at the specified position in this list.
     * @param index the index of the element to be incremented
     * @return the previous element at the specified position in this list
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int getAndIncrement(int index) {
        if (index < 0 || index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return values[index]++;
    }

    /**
     * Decrements by one the element at the specified position in this list.
     * @param index the index of the element to be decremented
     * @return the previous element at the specified position in this list
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int getAndDecrement(int index) {
        if (index < 0 || index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return values[index]--;
    }

    /**
     * Increments by one the element at the specified position in this list.
     * @param index the index of the element to be incremented
     * @return the updated element at the specified position in this list
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int incrementAndGet(int index) {
        if (index < 0 || index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return ++values[index];
    }

    /**
     * Decrements by one the element at the specified position in this list.
     * @param index the index of the element to be decremented
     * @return the updated element at the specified position in this list
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int decrementAndGet(int index) {
        if (index < 0 || index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return --values[index];
    }

    /**
     * Returns the number of elements in this list.
     * @return the number of elements in this list
     */
    public int size() {
        return size;
    }

    /**
     * Returns {@code true} if this list has no elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if this list has no elements
     */
    public boolean isEmpty() {
        return size==0;
    }

    /**
     * Sorts the elements of this list in increasing order
     */
    public void sort() {
        Arrays.sort(values, 0, size);
    }

    /**
     * Returns an index associated with the specified value, or
     * {@code (-insertionPoint - 1)} if the value is not found.  The returned
     * value is undefined if the list of integers is not sorted in
     * increasing order.
     * @param value the value to be searched for
     * @return an index associated with the specified value, or
     * {@code (-insertionPoint - 1)} if the value is not found
     */
    public int binarySearch(int value) {
        return Arrays.binarySearch(values, 0, size, value);
    }

    /**
     * Copies and returns this list of integers, truncating or padding with 0
     * as necessary so that the copy has the specified length.
     * @param newLength the length of the returned array
     * @return a copy of this list of integers, truncating or padding with 0
     * as necessary so that the copy has the specified length
     * @throws NegativeArraySizeException if {@code newLength < 0}
     */
    public int[] copyOf(int newLength) {
        int[] ia = Arrays.copyOf(values, newLength);
        if (newLength > size) {
            Arrays.fill(ia, size, newLength, 0);
        }
        return ia;
    }

    /**
     * Truncates this list of integer by removing all elements whose
     * index is greater than or equal to the specified size.  The list
     * is unchanged if the specified size is greater than the number of
     * elements in the list.
     * @param newSize the number of elements in the truncated list
     * @throws IllegalArgumentException if {@code newSize < 0}
     */
    public void truncate(int newSize) {
        if (newSize < 0) {
            throw new IllegalArgumentException(String.valueOf(newSize));
        }
        if (newSize < size) {
            size = newSize;
        }
    }

    /**
     * Returns an integer array containing the sequence of elements in this
     * list.
     * @return an integer array containing the sequence of elements in this
     * list
     */
    public int[] toArray() {
        return Arrays.copyOf(values, size);
    }

    /**
     * Returns an {@code IntStream} containing the sequence of elements in this
     * list.  The contract for this method is unspecified if this
     * list is modified during use of the returned stream.
     * @return an {@code IntStream} containing the sequence of elements in this
     * list
     */
    public IntStream stream() {
        return Arrays.stream(values, 0, size);
    }

    /**
     * Removes all elements from this list.
     */
    public void clear() {
        this.size = 0;
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.toArray())}
     *
     * @return {@code java.util.Arrays.toString(this.toArray())}
     */
    @Override
    public String toString() {
        return Arrays.toString(toArray());
    }
}
