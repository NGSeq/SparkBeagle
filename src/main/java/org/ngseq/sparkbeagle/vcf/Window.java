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
package org.ngseq.sparkbeagle.vcf;

import java.util.ArrayList;
import java.util.List;

/**
 * Class {@code Window} represents a window of VCF recList.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 *
 * @param <E> the type of elements in this window
 */
public class Window<E extends GTRec> {

    private final List<E> recList;
    private final int overlapEnd;
    private final int overlapStart;
    private final boolean lastWindowOnChrom;
    private final boolean lastWindow;

    /**
     * Constructs a new {@code Window} instance from the specified data.
     * The contract for the constructed instance is undefined if any element
     * of the specified {@code recList} is {@code null}.
     * @param recList a list of marker recList
     * @param overlapEnd the index of the first marker after the overlap with
     * the preceding marker window
     * @param overlapStart the index of the first marker in the overlap with the
     * next marker window
     * @param lastWindowOnChrom  {@code true} if the sliding window of
     * VCF Records is the last window for its chromosome
     * @param lastWindow {@code true} if the sliding window of
     * VCF Records is the last window
     * @throws IllegalArgumentException if {@code recList.isEmpty()}
     * @throws IllegalArgumentException if
     * {@code lastWindowOnChrom && overlapStart != recList.size()}
     * @throws NullPointerException if {@code recList==null}
     */
    public Window(List<E> recList, int overlapEnd, int overlapStart,
            boolean lastWindowOnChrom, boolean lastWindow) {
        if (recList.isEmpty()) {
            throw new IllegalArgumentException(recList.toString());
        }
        if (lastWindowOnChrom && overlapStart!=recList.size()) {
            throw new IllegalArgumentException(String.valueOf(overlapStart));
        }
        this.recList = new ArrayList<>(recList);
        this.overlapEnd = overlapEnd;
        this.overlapStart = overlapStart;
        this.lastWindowOnChrom = lastWindowOnChrom;
        this.lastWindow = lastWindow;
    }

    /**
     * Returns the number of markers in this window.
     * @return the number of markers in this window
     */
    public int nMarkers() {
        return recList.size();
    }

    /**
     * Returns the list of recList in this window.
     * @return the list of recList in this window
     */
    public List<E> recList() {
        return new ArrayList<>(recList);
    }

    /**
     * Returns the specified record.
     * @param marker a marker index
     * @return the specified record
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nMarkers()}
     */
    public E rec(int marker) {
        return recList.get(marker);
    }

    /**
     * Adds the specified records to he specified list.
     * @param list the list to be added to
     * @param start the start record index (inclusive)
     * @param end the end record index (exclusive)
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end > this.nMarkers() || start > end}
     * @throws NullPointerException if {@code list == null}
     */
    public void addRecords(List<E> list, int start, int end) {
        list.addAll(recList.subList(start, end));
    }

    /**
     * Returns the index of the first marker after the overlap with the
     * preceding marker window. Returns 0 if the current window
     * is the first window.
     *
     * @return the index of the first marker after the overlap with the
     * preceding marker window
     */
    public int overlapEnd() {
        return overlapEnd;
    }

    /**
     * Returns the index of the first marker in the overlap with the
     * next marker window. Returns {@code this.size()} if the next marker
     * window does not exist or is from a different chromosome.
     * @return the first marker index in the overlap between this
     * marker window and the next marker window
     */
    public int overlapStart() {
        return overlapStart;
    }

    /**
     * Returns the chromosome index of the first maker in the window.
     * @return the chromosome index of the first maker in the window
     */
    public int chromIndex() {
        return recList.get(0).marker().chromIndex();
    }

    /**
     * Returns {@code true} if the sliding window of genotype records is the
     * last window for its chromosome and returns {@code false} otherwise.
     * @return {@code true} if the sliding window of genotype records is the
     * last window for its chromosome
     */
    public boolean lastWindowOnChrom() {
        return lastWindowOnChrom;
    }

    /**
     * Returns {@code true} if the sliding window of genotype records is the
     * last window and returns {@code false} otherwise.
     * @return {@code true} if the sliding window of genotype records is the
     * last window
     */
    public boolean lastWindow() {
        return lastWindow;
    }
}
