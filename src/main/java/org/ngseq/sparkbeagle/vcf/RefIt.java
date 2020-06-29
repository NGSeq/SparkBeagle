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
import org.ngseq.sparkbeagle.blbutil.*;
import org.ngseq.sparkbeagle.bref.SeqCoder3;

import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.function.Function;
import java.util.stream.Collectors;

import static org.ngseq.sparkbeagle.bref.SeqCoder3.COMPRESS_FREQ_THRESHOLD;
import static org.ngseq.sparkbeagle.bref.SeqCoder3.MAX_NALLELES;

/**
 * <p>Class {@code RefIt} represents  an iterator whose {@code next()}
 * method returns an object storing data from a VCF record with
 * phased, non-missing genotypes.
 * </p>
 * <p>Instances of class {@code RefIt} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefIt implements SampleFileIt<RefGTRec> {

    /**
     * The default number of {@code GTRec} objects that are
     * stored in a buffer.
     */
    public static final int MAX_EM_BUFFER_SIZE = 10000;

    private final VcfHeader vcfHeader;
    private final FileIt<String> strIt;
    private final Function<String, RefGTRec> mapper;
    private final Filter<Marker> markerFilter;
    private final Thread fileReaderThread;
    private volatile boolean stopFileReadingThread = false;

    private final BlockingQueue<String[]> stringBuffers;
    private final Deque<RefGTRec> emBuffer;

    private final List<RefGTRec> buffer;
    private final SeqCoder3 seqCoder;

    /**
     * Create and returns a new {@code RefIt} instance from the specified
     * iterator.
     * @param strIt an iterator that returns lines of a VCF file
     * @return a new {@code RefIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if {@code strIt == null}
     */
    public static RefIt create(FileIt<String> strIt) {
        return RefIt.create(strIt, Filter.acceptAllFilter(),
                Filter.acceptAllFilter(), MAX_EM_BUFFER_SIZE);
    }

    /**
     * Create and returns a new {@code RefIt} instance from the specified
     * objects.
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param bufferSize the buffer size
     * @return a new {@code RefIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strItt}
     * @throws IllegalArgumentException if {@code bufferSize < 1}
     * @throws NullPointerException if {@code strIt == null}
     */
    public static RefIt create(FileIt<String> strIt,
            Filter<String> sampleFilter, Filter<Marker> markerFilter,
            int bufferSize) {
        RefIt refIt = new RefIt(strIt, sampleFilter, markerFilter, bufferSize);
        refIt.start();
        return refIt;
    }

    private RefIt(FileIt<String> strIt, Filter<String> sampleFilter,
            Filter<Marker> markerFilter, int bufferSize) {
        if (bufferSize < 1) {
            throw new IllegalArgumentException(String.valueOf(bufferSize));
        }
        if (markerFilter==null) {
            markerFilter = Filter.acceptAllFilter();
        }
        this.vcfHeader = new VcfHeader(strIt, sampleFilter);
        this.strIt = strIt;
        this.mapper = (String s) -> {
            return RefGTRec.alleleCodedInstance(new VcfRecGTParser(vcfHeader, s));
        } ;
        this.markerFilter = markerFilter;
        this.stringBuffers = new ArrayBlockingQueue<>(1);
        this.emBuffer = new ArrayDeque<>(MAX_EM_BUFFER_SIZE);
        this.buffer = new ArrayList<>();
        this.seqCoder = new SeqCoder3(vcfHeader.samples());
        this.fileReaderThread = fileReadingThread();
    }

    private void start() {
        this.fileReaderThread.setDaemon(true);
        this.fileReaderThread.start();
        fillEmissionBuffer();
        if (emBuffer.isEmpty()) {
            noRecordFoundError(strIt);
        }
    }

    private void noRecordFoundError(FileIt<String> it) {
        if (it.hasNext()==false) {
            StringBuilder sb = new StringBuilder(100);
            sb.append("No VCF records found (data source: ");
            sb.append(it.file()==null ? "stdin" : it.file());
            sb.append(")");
            sb.append(Const.nl);
            sb.append("Check that the chromosome identifiers are the same in each input VCF");
            sb.append(Const.nl);
            sb.append("file and in the \'chrom=\' command line argument (if \'chrom=\' is used).");
            throw new IllegalArgumentException(sb.toString());
        }
    }

    private Thread fileReadingThread() {
        Runnable runnable = () -> {
            try {
                String line = readLine(strIt);
                int bufferSize = stringBufferSize(line);
                while (line != null && stopFileReadingThread == false) {
                    String chromPlusTab = chromFieldPlusTab(line);
                    String[] sa = new String[bufferSize];
                    int size = 0;
                    while (line != null && size < bufferSize
                            && line.startsWith(chromPlusTab)) {
                        sa[size++] = line;
                        line = readLine(strIt);
                    }
                    if (size < bufferSize) {
                        sa = Arrays.copyOf(sa, size);
                    }
                    putInBlockingQueue(stringBuffers, sa);
                }
                if (stopFileReadingThread == false) {
                    putInBlockingQueue(stringBuffers, new String[0]);    // sentinel
                }
            }
            catch (Throwable e) {
                Utilities.exit(e);
            }
        };
        return new Thread(runnable);
    }

    private static int stringBufferSize(String line) {
        if (line == null) {
            return 0;
        }
        long nBytesPerLine = 2*line.length();
        Runtime rt = Runtime.getRuntime();
        long maxMem = rt.maxMemory();
        if (maxMem == Long.MAX_VALUE) {
            maxMem = 500 * (1 << 30);
        }
        long bufferSize = 1 + (maxMem / (20*nBytesPerLine));
        if (bufferSize > MAX_EM_BUFFER_SIZE) {
            bufferSize = MAX_EM_BUFFER_SIZE;
        }
        return (int) bufferSize;
    }

    private static <E> void putInBlockingQueue(BlockingQueue<E> q, E e) {
        try {
            q.put(e);
        } catch (InterruptedException ex) {
            Utilities.exit("Error: InterruptedException", ex);
        }
    }

    private static <E> E takeFromBlockingQueue(BlockingQueue<E> q) {
        try {
            return q.take();
        } catch (InterruptedException ex) {
            Utilities.exit("Error: InterruptedException", ex);
        }
        assert false;
        return null;
    }

    private static String chromFieldPlusTab(String vcfRecord) {
        int tabIndex = vcfRecord.indexOf(Const.tab);
        if (tabIndex == -1) {
            String s = Const.nl + "ERROR: Missing tab delimiter in VCV Record:"
                    + Const.nl + vcfRecord
                    + Const.nl + "Exiting Program";
            Utilities.exit(s);
        }
        return vcfRecord.substring(0, tabIndex + 1);
    }

    private void fillEmissionBuffer() {
        assert emBuffer.isEmpty();
        int lastLength = -1;
        while (lastLength != 0 && emBuffer.size() < MAX_EM_BUFFER_SIZE) {
            String[] stringBuffer = takeFromBlockingQueue(stringBuffers);
            lastLength = stringBuffer.length;
            if (stringBuffer.length>0) {
                List<RefGTRec> list = Arrays.stream(stringBuffer)
                        .parallel()
                        .map(mapper)
                        .filter(e -> markerFilter.accept(e.marker()))
                        .collect(Collectors.toList());
                for (int j=0, n=list.size(); j<n; ++j) {
                    if (buffer.size()==Integer.MAX_VALUE) {
                        flushBuffer();
                    }
                    RefGTRec e = list.get(j);
                    if (applySeqCoding(e)==false) {
                        buffer.add(e);
                    }
                    else {
                        boolean success = seqCoder.add(e);
                        if (success == false) {
                            flushBuffer();
                            success = seqCoder.add(e);
                            assert success;
                        }
                        buffer.add(null);
                    }
                }
            }
            else {
                // put sentinel element back
                putInBlockingQueue(stringBuffers, stringBuffer);
            }
        }
        if (lastLength==0) {
            flushBuffer();
        }
    }

    private void flushBuffer() {
        List<RefGTRec> list = seqCoder.getCompressedList();
        int index = 0;
        for (int j=0, n=buffer.size(); j<n; ++j) {
            GTRec ve = buffer.get(j);
            if (ve==null) {
                buffer.set(j, list.get(index++));
            }
        }
        emBuffer.addAll(buffer);
        buffer.clear();
    }

    private static String readLine(FileIt<String> it) {
        if (it.hasNext()==false) {
            return null;
        }
        String line = it.next();
        while (line.trim().isEmpty() && it.hasNext()) {
            line = it.next();
        }
        return line;
    }

    @Override
    public void close() {
        stopFileReadingThread = true;
        stringBuffers.poll();   // unblock file reading thread
        try {
            fileReaderThread.join();
        } catch (InterruptedException ex) {
            Utilities.exit("Error: InterruptedException", ex);
        }
        strIt.close();
        emBuffer.clear();
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return !emBuffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public RefGTRec next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        RefGTRec first = emBuffer.removeFirst();
        if (emBuffer.isEmpty()) {
            fillEmissionBuffer();
        }
        return first;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException(this.getClass().toString());
    }

    @Override
    public FSDataInputStream file() {
        return strIt.file();
    }

    @Override
    public Samples samples() {
        return vcfHeader.samples();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append("RefVcfIt from file: ");
        sb.append(strIt.file()==null ? "stdin" : strIt.file().toString());
        return sb.toString();
    }

    public static RefGTRec recodeIfLowFreq(RefGTRec rec) {
        int nHaps = 2*rec.nSamples();
        int[] alleleCounts = alleleCounts(rec);
        int majorCnt = max(alleleCounts);
        if (COMPRESS_FREQ_THRESHOLD < (1.0f + majorCnt)/nHaps) {
            int[][] hapIndices = hapIndices(rec, alleleCounts);
            return RefGTRec.alleleCodedInstance(rec.marker(), rec.samples(),
                    hapIndices);
        }
        else {
            return rec;
        }
    }

    private static int[][] hapIndices(GTRec ve, int[] alCnts) {
        int majorAllele = majorAllele(alCnts);
        int[][]  hapIndices = new int[alCnts.length][];
        for (int j=0; j<hapIndices.length; ++j) {
            hapIndices[j] = (j == majorAllele) ? null : new int[alCnts[j]];
        }
        int[] indices = new int[alCnts.length];
        for (int h=0, n=ve.size(); h<n; ++h) {
            int a = ve.get(h);
            if (a != majorAllele) {
                hapIndices[a][indices[a]++] = h;
            }
        }
        return hapIndices;
    }

    private static int majorAllele(int[] alleleCnts) {
        int major = 0;
        for (int j=1; j<alleleCnts.length; ++j) {
            if (alleleCnts[j] > alleleCnts[major]) {
                major = j;
            }
        }
        return major;
    }

    private static int max(int[] ia) {
        int maxIndex = 0;
        for (int j=1; j<ia.length; ++j) {
            if (ia[j] > ia[maxIndex]) {
                maxIndex = j;
            }
        }
        return ia[maxIndex];
    }

    private static int[] alleleCounts(GTRec ve) {
        int[] cnts = new int[ve.marker().nAlleles()];
        for (int h=0, n=ve.size(); h<n; ++h) {
            ++cnts[ve.get(h)];
        }
        return cnts;
    }

    private static boolean applySeqCoding(RefGTRec rec) {
        assert rec.isAlleleCoded();
        if (rec.marker().nAlleles() > MAX_NALLELES) {
            return false;
        }
        int nHaps = rec.size();
        int majAllele = rec.majorAllele();
        int majCnt = nHaps;
        for (int a=0, n=rec.nAlleles(); a<n; ++a) {
            if (a!=majAllele) {
                majCnt -= rec.alleleCount(a);
            }
        }
        return (COMPRESS_FREQ_THRESHOLD >=(1.0f + majCnt)/nHaps);
    }
}
