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

import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CopyOnWriteArrayList;

/**
 * <p>Class {@code ChromIds} is a singleton class that represents a
 * list of chromosome identifiers.
 * </p>
 * The singleton instance of {@code ChromIds} is thread-safe.
 *
 * @author Brian L. Browning
 */
public final class ChromIds {

    private static final ChromIds chromIds = new ChromIds();

    private final ThreadSafeIndexer<String> indexer;

    private final ConcurrentMap<String, Integer> map;
    private volatile CopyOnWriteArrayList<String> ids;

    private ChromIds() {
        // private constructor to restrict instantiation.
        int initCapacity = 4;
        this.indexer = new ThreadSafeIndexer<>(initCapacity);

        this.map = new ConcurrentHashMap<>(initCapacity);
        this.ids = new CopyOnWriteArrayList<>(new String[0]);
    }

    /**
     * Returns the singleton {@code ChromIds} instance.
     * @return the singleton {@code ChromIds} instance
     */
    public static ChromIds instance() {
        return chromIds;
    }

    /**
     * Returns the index of the specified chromosome identifier.  If
     * the chromosome identifiers is not yet indexed, the chromosome identifier
     * will be indexed. Chromosome identifier indices are assigned in
     * consecutive order beginning with 0.
     * @param id a chromosome identifier
     * @return the index of the specified chromosome identifier
     * @throws IllegalArgumentException if {@code id.isEmpty()}
     * @throws NullPointerException if {@code id == null}
     */
    public int getIndex(String id) {
        if (id.isEmpty()) {
            throw new IllegalArgumentException("id.isEmpty()");
        }
        Integer index = map.get(id);
        if (index!=null) {
            return index;
        }
        else {
            int i = indexer.getIndex(id);
            map.putIfAbsent(id, i);
            return i;
        }
    }

    /**
     * Returns an array of chromosome identifier indices corresponding to the
     * specified array of chromosome identifiers.  If a sample identifier is
     * not yet indexed, the sample identifier will be indexed.  Sample
     * identifier indices are assigned in increasing order starting with 0.
     * @param ids an array of sample identifiers
     * @return an array of sample identifier indices
     * @throws IllegalArgumentException if there is a {@code j} satisfying
     * {@code (0 <= j && j < ids.length) && ids[j].isEmpty()}
     * @throws NullPointerException if {@code ids == null}
     * @throws NullPointerException if there is a {@code j} satisfying
     * {@code (0 <= j && j < ids.length) && (ids[j] == null)}
     */
    public int[] getIndices(String[] ids) {
        for (String id : ids) {
            if (id.isEmpty()) {
                throw new IllegalArgumentException("id.isEmpty()");
            }
        }
        return indexer.getIndices(ids);
    }

    /**
     * Returns the index of the specified chromosome identifier, or returns
     * {@code -1} if the specified chromosome identifier is not indexed.
     *
     * @param id a chromosome identifier.
     * @return the index of the specified chromosome identifier, or
     * {@code -1} if the specified chromosome identifier is not indexed.
     *
     * @throws IllegalArgumentException if {@code id.isEmpty()}
     * @throws NullPointerException if {@code id == null}
     */
    public int getIndexIfIndexed(String id) {
        if (id.isEmpty()) {
            throw new IllegalArgumentException("id.isEmpty()");
        }
        Integer index = map.get(id);
        if (index!=null) {
            return index;
        }
        else {
            int i = indexer.getIndexIfIndexed(id);
            if (i>=0) {
                map.putIfAbsent(id, i);
            }
            return i;
        }
    }

    /**
     * Returns the number of indexed chromosomes identifiers.
     * @return the number of indexed chromosomes identifiers
     */
    public int size() {
        return indexer.size();
    }

    /**
     * Returns the chromosome identifier with the specified index.
     * @param index a chromosome identifier index.
     * @return the specified chromosome identifier.
     * @throws IndexOutOfBoundsException if
     * {@code  index < 0 || index >= this.size()}
     */
    public String id(int index) {
        if (index >= ids.size()) {
            for (int j=ids.size(), n=indexer.size(); j<n; ++j) {
                ids.add(indexer.item(j));
            }
        }
        return ids.get(index);
    }

    /**
     * Returns the list of chromosome identifiers as an array.
     * The returned array will have length {@code this.size()}, and
     * it will satisfy {@code this.ids()[k].equals(this.id(k)) == true}
     * for {@code  0 <= k < this.size()}.
     *
     * @return an array of chromosome identifiers
     */
    public String[] ids() {
        return indexer.items().toArray(new String[0]);
    }

    /**
     * Returns  {@code java.util.Arrays.toString(this.ids())}.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return Arrays.toString(this.ids());
    }
}
