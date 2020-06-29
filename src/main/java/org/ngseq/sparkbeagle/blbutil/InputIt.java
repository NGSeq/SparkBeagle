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

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.PositionedReadable;
import org.apache.hadoop.fs.Seekable;
import org.apache.hadoop.io.compress.CompressionCodec;

import java.io.*;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;


/**
 * <p>Class {@code InputIt} is a buffered iterator whose {@code next()}
 * method returns lines of a text input stream.
 * </p>
 * <p>If an {@code IOException} is thrown when an {@code InputIt}
 * instance reads from the text input stream, the {@code IOException}
 * is trapped, an error message is written to standard out, and the
 * Java Virtual Machine is terminated.
 * </p>
 * Instances of class {@code InputIt} are not thread-safe.
 *
 * @see #DEFAULT_BUFFER_SIZE
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class InputIt implements FileIt<String> {

    /**
     * The default buffer size, which is 4,194,304 bytes.
     */
    public static final int DEFAULT_BUFFER_SIZE = 1<<26;

    private FSDataInputStream file = null;
    private final BufferedReader in;
    private String next = null;

    /**
     * Constructs a new {@code InputStreamIterator} with default buffer
     * size that will iterate through lines of the specified input stream.
     *
     * @param is input stream of text data
     *
     * @see #DEFAULT_BUFFER_SIZE
     */
    private InputIt(InputStream is, FSDataInputStream fis, int bufferSize) {
        BufferedReader br = null;
        try {
            InputStreamReader isr = new InputStreamReader(is);
            br = new BufferedReader(isr, bufferSize);
            next = br.readLine();
        }
        catch(IOException e) {
            Utilities.exit("Error reading " + is, e);
        }
        this.in = br;
        this.file = fis;
    }

    public InputIt(ByteArrayOutputStream bos, int defaultBufferSize) {


        byte[] barr = bos.toByteArray();
        ByteArrayInputStream sbai = new ByteArrayInputStream(barr);
        BlockCompressedInputStream is = new BlockCompressedInputStream(sbai);

        BufferedReader br = null;

        try {
            InputStreamReader isr = new InputStreamReader(is);
            br = new BufferedReader(isr, defaultBufferSize);
            next = br.readLine();
        }
        catch(IOException e) {
            Utilities.exit("Error reading " + is, e);
        }
        this.in = br;

        try {
            this.file =  new FSDataInputStream(new SeekableByteArrayInputStream(barr));
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    @Override
    public FSDataInputStream file() {
        return file;
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return (next != null);
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public String next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        String current = next;
        try {
            next = in.readLine();
        }
        catch (IOException e) {
            Utilities.exit("Error reading " + in, e);
        }
        return current;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        String s = this.getClass().toString() + ".remove()";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public void close() {
        try {
            in.close();
        }
        catch (IOException e) {
            Utilities.exit("Error closing " + in, e);
        }
        next=null;
    }

    /**
     * Returns a string representation of this iterator.  The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of this iterator
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(200);
        sb.append("[file= ");
        sb.append(file);
        sb.append("; next=\"");
        sb.append(next);
        sb.append("\"]");
        return sb.toString();
    }



    public static FileIt<String> fromGzipFile(FSDataInputStream targFile, CompressionCodec codec, boolean zipped) {

        try {

            if (zipped) {
                return new InputIt(codec.createInputStream(targFile), targFile, DEFAULT_BUFFER_SIZE);
            }
            else {
                return new InputIt(targFile, targFile, DEFAULT_BUFFER_SIZE);
            }
        }
        catch(Exception e) {
            Utilities.exit("Error reading " + targFile, e);
        }
        assert false;
        return null;
    }

    public static InputIt fromGzipFile(FSDataInputStream fis, boolean zipped) {
        try {

            if (zipped) {

                return new InputIt(new GZIPInputStream(fis), fis, DEFAULT_BUFFER_SIZE);
            }
            else {
                return new InputIt(fis, fis, DEFAULT_BUFFER_SIZE);
            }
        }
        catch(Exception e) {
            Utilities.exit("Error reading " + fis, e);
        }
        assert false;
        return null;
    }

    public static InputIt fromBGZFStream(ByteArrayOutputStream bos, boolean zipped) {
        try {

                return new InputIt(bos, DEFAULT_BUFFER_SIZE);

        }
        catch(Exception e) {
            Utilities.exit("Error reading " + bos, e);
        }
        assert false;
        return null;
    }

    public static InputIt fromIterator(Iterator<String> iterator, String header) {

        try {
            return new InputIt(iterator, header, DEFAULT_BUFFER_SIZE);
        } catch (IOException e) {
            e.printStackTrace();
        }
        assert false;
        return null;
    }

    private InputIt(Iterator<String> iterator, String header, int bufferSize) throws IOException {


        //Check wheater the partition has header, add if not
        ByteArrayOutputStream bao = new ByteArrayOutputStream(bufferSize);

        while(iterator.hasNext()) {
            try {
                bao.write(iterator.next().getBytes("UTF-8"));
                bao.write("\n".getBytes("UTF-8"));
            } catch (IOException e) {
                e.printStackTrace();
            }catch (OutOfMemoryError oe) {
                oe.printStackTrace();
            }
        }



        byte[] barr = new byte[bufferSize];
        InputStream is = new ByteArrayInputStream(barr);

        BufferedReader br = null;
        try {
            InputStreamReader isr = new InputStreamReader(is);
            br = new BufferedReader(isr, bufferSize);
            next = br.readLine();
        }
        catch(IOException e) {
            Utilities.exit("Error reading " + is, e);
        }
        this.in = br;

        this.file =  new FSDataInputStream(new SeekableByteArrayInputStream(barr));

    }

    static class SeekableByteArrayInputStream extends ByteArrayInputStream implements Seekable, PositionedReadable {

        public SeekableByteArrayInputStream(byte[] buf)
        {
            super(buf);
        }
        @Override
        public long getPos() throws IOException{
            return pos;
        }

        @Override
        public void seek(long pos) throws IOException {
            if (mark != 0)
                throw new IllegalStateException();

            reset();
            long skipped = skip(pos);

            if (skipped != pos)
                throw new IOException();
        }

        @Override
        public boolean seekToNewSource(long targetPos) throws IOException {
            return false;
        }

        @Override
        public int read(long position, byte[] buffer, int offset, int length) throws IOException {

            if (position >= buf.length)
                throw new IllegalArgumentException();
            if (position + length > buf.length)
                throw new IllegalArgumentException();
            if (length > buffer.length)
                throw new IllegalArgumentException();

            System.arraycopy(buf, (int) position, buffer, offset, length);
            return length;
        }

        @Override
        public void readFully(long position, byte[] buffer) throws IOException {
            read(position, buffer, 0, buffer.length);

        }

        @Override
        public void readFully(long position, byte[] buffer, int offset, int length) throws IOException {
            read(position, buffer, offset, length);
        }
    }



    private static boolean isBGZipFile(File file) throws IOException {
        try (InputStream is=new BufferedInputStream(new FileInputStream(file))) {
		return BlockCompressedInputStream.isValidFile(is);
	    }
    }

     /**
     * Constructs and returns an {@code InputIt} instance with the specified
     * buffer size that iterates through lines of the specified text file.
     *
     * @param file a text file
     * @param bufferSize the buffer size in bytes
     * @return a new {@code InputIt} instance that iterates through
     * lines of the specified text file
     *
     * @throws IllegalArgumentException if {@code bufferSize < 0}
     * @throws NullPointerException if {@code filename == null}
     */

}
