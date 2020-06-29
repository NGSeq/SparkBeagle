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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code IntMap} represents a map with integer keys and generic type
 * values.
 * </p>
 * <p>Class {@code IntMap} is not thread-safe.
 * </p>
 *
 * @param <E> the type of values in this map
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IntMap<E> {

    private static final int NIL = -1;
    private static final float LOAD_FACTOR = 0.75f;

    private int size;
    private int nBuckets;

    private int[] next;
    private int[] data; // stores list index of keys and values
    private int[] keys;
    private List<E> values;
    private int firstFreeIndex;

    /**
     * Creates a new {@code IntMap} instance.
     *
     * @param capacity the initial capacity of this map
     * @throws IllegalArgumentException if
     * {@code capacity < 0 || (capacity > (1 << 30))}
     */
    public IntMap(int capacity) {
        if (capacity < 1 || capacity > (1<<30)) {
            throw new IllegalArgumentException(String.valueOf(capacity));
        }
        int numBuckets = (int) Math.ceil(capacity/LOAD_FACTOR) + 1;
        allocateArrays(capacity, numBuckets);
        values = new ArrayList<>(capacity);
        initializeFields(numBuckets);
    }

    private void allocateArrays(int capacity, int numBuckets) {
        this.next = new int[numBuckets + capacity];
        this.data = new int[numBuckets + capacity];
        this.keys = new int[capacity];
    }

    private void initializeFields(int numBuckets) {
        size = 0;
        nBuckets = numBuckets;
        firstFreeIndex = nBuckets;
        Arrays.fill(next, 0, nBuckets, NIL);
        for (int j=nBuckets; j<next.length; ++j) {
            next[j] = j+1;
        }
        values.clear();
    }

    /*
     * Increases the capacity of the internal hash table.
     */
    private void rehash(int newCapacity) {
        if (newCapacity > size) {
            int oldSize = size;
            int[] oldKeys = keys.clone();
            List<E> oldValues = new ArrayList<>(values);
            int newNumBuckets = (int) Math.ceil(newCapacity/LOAD_FACTOR);
            allocateArrays(newCapacity, newNumBuckets);
            initializeFields(newNumBuckets);
            for (int j=0; j<oldSize; ++j) {
                put(oldKeys[j], oldValues.get(j));
            }
        }
    }

    /**
     * Removes all keys from this map.
     */
    public void clear() {
        initializeFields(nBuckets);
    }

    /**
     * Returns {@code true} if the map contains the specified key,
     * and returns {@code false} otherwise.
     * @param key a key
     * @return {@code true} if the map contains the specified key
     */
    public boolean contains(int key) {
        return indexOf(key)>=0;
    }

    /**
     * Returns the current index of the specified key, or returns {@code -1}
     * if the key is not found in this map.
     * @param key a key index
     * @return the current index of the specified key
     */
    private int indexOf(int key) {
        int index = next[bucket(key)];
        while (index!=NIL && keys[data[index]]<key) {
            index = next[index];
        }
        return (index!=NIL && keys[data[index]]==key) ? index : -1;
    }

    /**
     * Adds the specified key and value to this map.  If there exists a value
     * in the map for the specified key, the old value is replaced with the
     * specified value. The indexing of keys immediately before and after this
     * method is invoked may differ if this map is changed by this operation.
     * @param key the key
     * @param value the value
     * @return the previous value associated with the key, or {@code null} if
     * there did not exist a value for the specified key
     */
    public E put(int key, E value) {
        int prevIndex = prevIndex(key);
        int nextIndex = next[prevIndex];
        if (nextIndex==NIL || keys[data[nextIndex]]!=key) {
            int index = firstFreeIndex;
            firstFreeIndex = next[firstFreeIndex];
            next[prevIndex] = index;
            data[index] = size;
            next[index] = nextIndex;
            keys[size] = key;
            values.add(value);
            ++size;
            if (size == keys.length) {
                int newCapacity = 3*keys.length/2 + 1;
                rehash(newCapacity);
            }
            return null;
        }
        else  {
            return values.set(data[nextIndex], value);
        }
    }

    /**
     * Removes the specified key from this map, and returns the value
     * previously associated with the key or {@code null} if no such
     * value exists. The indexing of keys immediately before and after
     * this method is invoked may differ if this map is changed by
     * this operation.
     *
     * @param key a key index
     * @return the previous value associated with the key, or {@code null} if
     * there did not exist a value for the specified key
     */
    public E remove(int key) {
        int prevIndex = prevIndex(key);
        int index = next[prevIndex];
        if (index==NIL || keys[data[index]]!=key) {
            return null;
        }
        else {
            E prevValue = values.get(data[index]);
            int oldListIndex = data[index];
            next[prevIndex] = next[index];
            next[index] = firstFreeIndex;
            firstFreeIndex = index;

            --size;
            if (oldListIndex!=size) {
                // overwrite removed key
                index = indexOf(keys[size]);
                data[index] = oldListIndex;
                keys[oldListIndex] = keys[size];
                values.set(oldListIndex, values.get(size));
                values.remove(size);
            }
            return prevValue;
        }
    }

    private int bucket(int key) {
        return Math.abs((71*key) % nBuckets);
    }

    private int prevIndex(int key) {
        int prevIndex = bucket(key);
        int index = next[prevIndex];
        while (index!=NIL && keys[data[index]]<key) {
            prevIndex = index;
            index = next[index];
        }
        return prevIndex;
    }

    /**
     * Returns the specified key.
     * @param index an index of a key in this map
     * @return the specified key
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int key(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return keys[index];
    }

    /**
     * Returns the value for the specified key or {@code null}
     * if the specified key is not present in this map.
     * @param key the key
     * @return the specified value
     */
    public E get(int key) {
        int index = indexOf(key);
        if (index == -1) {
            return null;
        }
        return values.get(data[index]);
    }

    /**
     * Returns the number of keys in this map.
     *
     * @return the number of keys in this map
     */
    public int size() {
        return size;
    }

    /**
     * Returns an array containing the keys in this map. The returned
     * array will satisfy:
     * {@code this.toArray()[j]==this.key(j)} for each
     * {@code j} satisfying {@code (0 <= j && j < this.size())}.
     * @return an array containing the keys in this map
     */
    public int[] keys() {
        return Arrays.copyOf(keys, size);
    }

    /**
     * Returns a list containing the values in this map. The returned
     * array will satisfy:
     * {@code this.get(this.keys()[j])==this.values.get(j)} for each
     * {@code j} satisfying {@code (0 <= j && j < this.size())}.
     * @return an array containing the values in this map
     */
    public List<E> values() {
        return new ArrayList<>(values);
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.keys())}.
     *
     * @return {@code java.util.Arrays.toString(this.keys())}
     */
    @Override
    public String toString() {
        return Arrays.toString(keys());
    }

    // xxx code for testing class
    private static void main(String[] args) {
        IntMap<Integer> map1 = new IntMap<>(4);
        java.util.Map<Integer, Integer> map2 = new java.util.HashMap<>(100);
        java.util.Random rand = new java.util.Random(0);
        int nTests = 5000;
        int maxKey = 100;
        for (int j=0; j<nTests; ++j) {
            int key = j % maxKey;
            int value = rand.nextInt();
            double d = rand.nextDouble();
            if (d < 0.005) {
                map1.clear();
                map2.clear();
            }
            else if (d < 0.4) {
                Integer i1 = map1.get(key);
                Integer i2 = map2.get(key);
                assert (i1==i2) || (i1!=null && i2!=null && i1.intValue()==i2);
            }
            else if (d < 0.7) {
                Integer i1 = map1.put(key, value);
                Integer i2 = map2.put(key, value);
                assert (i1==i2) || (i1!=null && i2!=null && i1.intValue()==i2);
            }
            else {
                Integer i1 = map1.remove(key);
                Integer i2 = map2.remove(key);
                assert (i1==i2) || (i1!=null && i2!=null && i1.intValue()==i2);
            }
            assert map1.size()==map2.size();
        }
    }
    // xx code for testing class
}
