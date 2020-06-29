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

import org.apache.hadoop.fs.FSDataInputStream;
import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.blbutil.Pair;
import org.ngseq.sparkbeagle.blbutil.SampleFileIt;
import org.ngseq.sparkbeagle.blbutil.Utilities;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.function.Supplier;

/**
 * <p>Class {@code WindowIt} represents a sliding window of VCF recList.
 </p>
 * <p>Instances of class {@code WindowIt} are not thread-safe.
 * </p>
 * @param <E> the type of elements in the window corresponding to this
 * window iterator.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}

 */
public class WindowIt<E extends GTRec> implements SampleFileIt<Window<E>> {

    private final FSDataInputStream file;                    // immutable
    private final Samples samples;              // immutable
    private final BlockingQueue<Window<E>> q;   // thread-safe
    private final Reader reader;                // thread-safe

    private Window<E> currentWind;              // thread-confined

    /**
     * Constructs and returns a new {@code WindowIt} instance for
     * the specified data.
     * @param <E> the type of elements in windows for the returned
     * window iterator.
     * @param supplier the object which supplies the {@code SampleFileIt} which
     * reads input data
     * @param genMap the genetic map
     * @param windowCM the requested window length in cM
     * @param overlapCM the requested overlap in cM between consecutive
     * windows on the same chromosome
     * @return a new {@code WindowIt} instance for the specified data
     * @throws IllegalArgumentException if a format error is detected in the
     * input file
     * @throws IllegalArgumentException if
     * {@code overlap < 0 || Float.isFinite(overlap) == false}
     * @throws IllegalArgumentException if
     * {@code window <= overlap || Float.isFinite(window) == false}
     * @throws NullPointerException if {@code it == null || genMap == null}
     */
    public static <E extends GTRec> WindowIt<E> newInstance(
            Supplier<SampleFileIt<E>> supplier, GeneticMap genMap,
            float windowCM, float overlapCM) {

        BlockingQueue<Window<E>> q = new ArrayBlockingQueue<>(1);
        BlockingQueue<Pair<FSDataInputStream,Samples>> q1 = new ArrayBlockingQueue<>(1);

        Reader<E> reader = new Reader(supplier, q, q1, genMap, windowCM,
                overlapCM);
        Thread t = new Thread(reader);
        t.setDaemon(true);
        t.start();

        Pair<FSDataInputStream,Samples> pair = takeFromQ(q1);
        return new WindowIt<>(pair.first(), pair.second(), reader, q);
    }

    /**
     * Constructs a new {@code WindowIt} instance.
     * @param it an iterator that returns VCF recList
     * @param genMap the genetic map
     * @param windowCM the requested window length in cM
     * @param overlapCM the requested overlap in cM between consecutive
     * windows on the same chromosome
     * @throws IllegalArgumentException if {@code it.hasNext() == false}
     * @throws IllegalArgumentException if a format error is detected in
     * a VCF record
     * @throws IllegalArgumentException if
     * {@code overlap < 0 || Float.isFinite(overlap) == false}
     * @throws IllegalArgumentException if
     * {@code window <= overlap || Float.isFinite(window) == false}
     * @throws NullPointerException if {@code it == null || genMap == null}
     */
    private WindowIt(FSDataInputStream file, Samples samples, Reader<E> reader,
            BlockingQueue<Window<E>> q) {
        this.file = file;
        this.samples = samples;
        this.q = q;
        this.reader = reader;
    }

    /**
     * Returns the genetic map.
     * @return the genetic map
     */
    public GeneticMap genMap() {
        return reader.genMap;
    }

    /**
     * Returns {@code true} if the sliding window of VCF recList can advance
     * and returns {@code false} otherwise.
     * @return {@code true} if the sliding window of VCF recList can advance
     */
    @Override
    public boolean hasNext() {
        return currentWind==null || currentWind.lastWindow()==false;
    }

    /**
     * Advances the sliding window of VCF recList, and returns the advanced
     * window as a {@code RefGTRec[]} object.
     *
     * @return the advanced window of VCF recList
     *
     * @throws IllegalArgumentException if a format error in the input data
     * is detected
     * @throws IllegalArgumentException if
     * {@code windowCm <= this.overlap() || Float.isFinite(windowCm) == false}
     * @throws IllegalStateException if
     * {@code this.canAdvanceWindow() == false}
     */
    @Override
    public Window<E> next() {
        if (hasNext()==false) {
            throw new IllegalStateException("canAdvanceWindow()==false");
        }
        currentWind = takeFromQ(q);
        return currentWind;
    }

    private static <E> E takeFromQ(BlockingQueue<E> q) {
        E e = null;
        try {
            e = q.take();
        } catch (InterruptedException ex) {
            throw new RuntimeException(ex);
        }
        return e;
    }

    /**
     * Returns the file from which VCF recList are read, or returns
     * {@code null} if the source is standard input.
     * @return the file from which VCF recList are read, or
     * {@code null} if the source is standard input
     */
    @Override
    public FSDataInputStream file() {
        return file;
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    @Override
    public Samples samples() {
        return samples;
    }

    /**
     * Releases any I/O resources controlled by this object.
     */
    @Override
    public void close() {
        reader.terminate();
        while (hasNext()) {
            next();
        }
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(1100);
        sb.append(this.getClass().toString());
        return sb.toString();
    }

    private static class Reader<V extends GTRec> implements Runnable {

        private final Supplier<SampleFileIt<V>> supplier;
        private final GeneticMap genMap;
        private final float windowCM;
        private final float overlapCM;
        private final BlockingQueue<Window<V>> q;
        private final BlockingQueue<Pair<File, Samples>> q1;
        private volatile boolean finished;

        public Reader(Supplier<SampleFileIt<V>> supplier,
                BlockingQueue<Window<V>> q, BlockingQueue<Pair<File, Samples>> q1,
                GeneticMap genMap, float windowCM, float overlapCM) {
            if (genMap==null) {
                throw new NullPointerException(GeneticMap.class.toString());
            }
            if (overlapCM < 0f || Float.isFinite(overlapCM)==false) {
                throw new IllegalArgumentException(String.valueOf(overlapCM));
            }
            if (windowCM <= overlapCM || Float.isFinite(windowCM)==false) {
                throw new IllegalArgumentException(String.valueOf(windowCM));
            }
            this.supplier = supplier;
            this.genMap = genMap;
            this.windowCM = windowCM;
            this.overlapCM = overlapCM;
            this.q = q;
            this.q1 = q1;
            this.finished = false;
        }

        public void terminate() {
            finished = true;
        }

        @Override
        public void run() {
            try {
                try (SampleFileIt<? extends V> it = supplier.get()) {
                    q1.add(new Pair(it.file(), it.samples()));
                    if (it.hasNext()==false) {
                        throw new IllegalArgumentException("No VCF records after filtering");
                    }
                    List<V> nextList = new ArrayList<>(10000);
                    V next = it.next();
                    Window<V> prevW = null;
                    double endCm = Double.NaN;
                    while (next!=null) {
                        int chromIndex = next.marker().chromIndex();
                        if (prevW==null || prevW.chromIndex()!=chromIndex) {
                            endCm = genMap.genPos(next.marker()) + windowCM;
                        }
                        else {
                            endCm += (windowCM - overlapCM);
                            prevW.addRecords(nextList, prevW.overlapStart(), prevW.nMarkers());
                        }
                        int endPos = genMap.basePos(chromIndex, endCm);
                        int overlapEnd = nextList.size();
                        nextList.add(next);   // ensure at least on record is added
                        next = it.hasNext() ? it.next() : null;
                        while (next!=null
                                && next.marker().chromIndex()==chromIndex
                                && next.marker().pos() < endPos) {
                            nextList.add(next);
                            next = it.hasNext() ? it.next() : null;
                        }
                        if (finished) {
                            next = null;
                        }
                        prevW = addWindowToQ(chromIndex, overlapEnd, nextList, next);
                        nextList.clear();
                    }
                }
            }
            catch (Throwable e) {
                Utilities.exit(e);
            }
        }

        private Window<V> addWindowToQ(int chromIndex, int overlapEnd,
                List<V> nextList, V next) {
            boolean last = (next==null);
            boolean lastOnChrom = next==null
                    || next.marker().chromIndex()!=chromIndex;
            int overlapStart = overlapStart(lastOnChrom, nextList, overlapCM);

            Window<V> window = new Window<>(nextList,
                    overlapEnd, overlapStart, lastOnChrom, last);
            try {
                q.put(window);
            } catch (InterruptedException ex) {
                throw new RuntimeException(ex);
            }
            return window;
        }

        private int overlapStart(boolean lastOnChrom, List<V> wind,
                float overlapCM) {
            if (wind.isEmpty() || lastOnChrom) {
                return wind.size();
            }
            else {
                Marker m = wind.get(wind.size()-1).marker();
                double endGenPos = genMap.genPos(m);
                double startGenPos = endGenPos - overlapCM;
                int key = genMap.basePos(m.chromIndex(), startGenPos);
                int low = 0;
                int high = wind.size()-1;
                while (low <= high) {
                    int mid = (low + high) >>> 1;
                    int midPos = wind.get(mid).marker().pos();
                    if (midPos < key) {
                        low = mid + 1;
                    }
                    else if (midPos > key) {
                        high = mid - 1;
                    }
                    else {
                        return firstIndexWithPos(wind, mid);
                    }
                }
                assert high < low;
                return firstIndexWithPos(wind, high);
            }
        }

        private int firstIndexWithPos(List<V> wind, int index) {
            if (index<0) {
                return 0;
            }
            else {
                int pos = wind.get(index).marker().pos();
                while (index>0 && wind.get(index-1).marker().pos()==pos) {
                    --index;
                }
                return index;
            }
        }
    }
}
