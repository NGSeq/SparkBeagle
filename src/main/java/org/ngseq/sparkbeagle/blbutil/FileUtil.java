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

import org.apache.hadoop.fs.FSDataInputStream;

import java.io.*;
import java.util.zip.GZIPOutputStream;

/**
 * Class {@code FileUtil} contains static methods for working with files.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class FileUtil {

    private FileUtil() {
        // private constructor prevents instantiation
    }

    /**
     * Creates and returns a new empty temporary file in the specified
     * directory.  The temporary file will be deleted when the JVM terminates
     * if the JVM terminates normally.
     * @param prefix a prefix of at least 3 characters for the temporary
     * filename
     * @param directory the directory in which the temporary file is to
     * be created.
     * @return a new empty temporary file in the specified directory.
     * @throws NullPointerException if
     * {@code prefix == null || directory == null}
     */
    public static File tempFile(String prefix, File directory) {
        if (prefix==null) {
            throw new NullPointerException(String.class.toString());
        }
        if (directory==null) {
            throw new NullPointerException(File.class.toString());
        }
        File tmpFile = null;
        try {
            tmpFile = File.createTempFile(prefix, null, directory);
            tmpFile.deleteOnExit();
        }
        catch (IOException ex) {
            Utilities.exit("Error opening file", ex);
        }
        return tmpFile;
    }

   /**
     * Returns a {@code java.io.RandomAccessFile} to read from or optionally
     * to write to the specified file.  If the input stream cannot
     * be opened, an error message will be printed and the java interpreter
     * will exit.
     * @param file a file
     * @param mode the access mode as described in the documentation for the
     * {@code java.io.RandomAccessFile} constructors
     * @return a {@code java.io.RandomAccessFile}
     * @throws NullPointerException if {@code file == null || mode == null}
     */
    public static RandomAccessFile randomAccessFile(File file, String mode) {
        RandomAccessFile raf = null;
        try {
            raf = new RandomAccessFile(file, mode);
        } catch (FileNotFoundException ex) {
            Utilities.exit("Error: file not found [" + file + "]", ex);
        }
        return raf;
    }

   /**
     * Returns an buffered {@code java.io.InputStream} with default buffer
     * size reading from the specified file.  If the input stream cannot
     * be opened, an error message will be printed and the java interpreter
     * will exit.
     * @param file a file
     * @return an buffered {@code java.io.InputStream}
     * @throws NullPointerException if {@code file == null}
     */
    public static InputStream bufferedInputStream(FSDataInputStream file) {
        InputStream is = null;
        try {
            is = new BufferedInputStream(file);
        } catch (Exception ex) {
            Utilities.exit("Error: FSDataStream not found [" + file + "]", ex);
        }
        return is;
    }

    /**
     * Returns an buffered {@code java.io.OutputStream} with default
     * buffer size writing to the specified file. If the file cannot
     * be opened for writing, an error message will be printed and the
     * Java Virtual Machine will exit.
     *
     * @param file a file
     * @return a buffered {@code java.io.OutputStream}
     * @throws NullPointerException if {@code file == null}
     */
    public static OutputStream bufferedOutputStream(File file) {
        OutputStream os = null;
        try {
            os = new BufferedOutputStream(new FileOutputStream(file));
        } catch (FileNotFoundException ex) {
            Utilities.exit("Error: file not found [" + file + "]", ex);
        }
        return os;
    }

    /**
     * Returns an buffered {@code java.io.OutputStream} with default
     * buffer size writing to the specified file. If the file cannot
     * be opened for writing, an error message will be printed and the
     * Java Virtual Machine will exit.  If the specified file exists and
     * {@code append} is {@code false}, bytes written by the returned
     * {@code java.io.OutputStream} will overwrite the previously existing file.
     *
     * @param file a file
     * @param append {@code true} if bytes will be appended to the end of
     * the file
     * @return an buffered {@code java.io.OutputStream}
     * @throws NullPointerException if {@code file == null}
     */
    public static OutputStream bufferedOutputStream(File file, boolean append) {
        OutputStream os = null;
        try {
            os = new BufferedOutputStream(new FileOutputStream(file, append));
        } catch (FileNotFoundException ex) {
            Utilities.exit("Error: file not found [" + file + "]", ex);
        }
        return os;
    }

    /**
     * Returns a {@code java.io.PrintWriter} that writes
     * to standard out.
     *
     * @return a {@code java.io.PrintWriter} that writes
     * to standard out
     */
    public static PrintWriter stdOutPrintWriter() {
        return new PrintWriter(
                new BufferedOutputStream(System.out));
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} writing to
     * the specified file.  The resulting file will be compressed using
     * the GZIP compression algorithm.  If the file cannot be opened, an
     * error message will be printed and the java interpreter will exit.
     * If the specified file exists, bytes written by the returned
     * {@code PrintWriter} will overwrite the previously existing file.
     * @param file a file
     * @return a {@code java.io.PrintWriter} writing to the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter gzipPrintWriter(File file) {
        PrintWriter out = null;
        try {
            FileOutputStream fos = new FileOutputStream(file);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            GZIPOutputStream gzos = new GZIPOutputStream(bos);
            out = new PrintWriter(gzos);
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} writing to
     * the specified file.  The resulting file will be compressed using
     * the GZIP compression algorithm.  If the file cannot be opened, an
     * error message will be printed and the java interpreter will exit.
     * If the specified file exists and {@code append} is {@code false}, bytes
     * written by the returned {@code PrintWriter} will overwrite the
     * previously existing file.
     * @param file a file
     * @param append {@code true} if bytes will be appended to the end of
     * the file
     * @return a {@code java.io.PrintWriter} writing to the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter gzipPrintWriter(File file, boolean append) {
        PrintWriter out = null;
        try {
            FileOutputStream fos = new FileOutputStream(file, append);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            GZIPOutputStream gzos = new GZIPOutputStream(bos);
            out = new PrintWriter(gzos);
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} that compresses
     * data using the BGZIP algorithm and writes the compressed data
     * to the specified file. The {@code close()} method of the returned
     * {@code PrintWriter} will write an empty BGZIP block to the end of the
     * output stream. If the file cannot be opened for writing, an error
     * message will be printed and the Java Virtual Machine will exit.
     * If the specified file exists, bytes written by the returned
     * {@code PrintWriter} will overwrite the previously existing file.
     *
     * @param file a file
     * @return a buffered {@code java.io.PrintWriter}
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter bgzipPrintWriter(File file) {
        boolean writeBuffer = true;
        OutputStream bos = FileUtil.bufferedOutputStream(file);
        return new PrintWriter(new BGZIPOutputStream(bos, writeBuffer));
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} that compresses
     * data using the BGZIP algorithm and writes the compressed data to
     * the specified file. The {@code close()} method of the returned
     * {@code PrintWriter} will write an empty BGZIP block to the end of the
     * output stream. If the file cannot be opened for writing, an error
     * message will be printed and the Java Virtual Machine will exit.
     * If the specified file exists and {@code append} is {@code false}, bytes
     * written by the returned {@code PrintWriter} will overwrite the
     * previously existing file.
     *
     * @param file a file
     * @param append {@code true} if bytes will be appended to the end of
     * the file
     * @return a buffered {@code java.io.PrintWriter}
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter bgzipPrintWriter(File file, boolean append) {
        boolean writeBuffer = true;
        OutputStream bos = bufferedOutputStream(file, append);
        return new PrintWriter(new BGZIPOutputStream(bos, writeBuffer));
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} writing to the
     * specified file.  If the file cannot be opened, an error message
     * will be printed and the Java Virtual Machine will exit.  If the specified
     * file exists, bytes written by  the returned {@code PrintWriter} will
     * overwrite the previously existing file.
     * @param file a file
     * @return a buffered {@code java.io.PrintWriter} writing to
     * the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter printWriter(File file) {
        return printWriter(file, false);
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} writing to
     * the specified file. If the file cannot be opened, an error message will
     * be printed and the Java Virtual Machine will exit.  If the specified
     * file exists and {@code append} is {@code false}, bytes written by the
     * returned {@code PrintWriter} will overwrite the previously existing file.
     *
     * @param file a file
     * @param append {@code true} if the data will be appended
     * to the end of any existing file
     * @return a buffered {@code java.io.PrintWriter} writing to
     * the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter printWriter(File file, boolean append) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(
                    new BufferedWriter(new FileWriter(file, append)));
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    /**
     * Returns an unbuffered {@code java.io.PrintWriter} writing to
     * the specified file. If the file cannot be opened, an error message will
     * be printed and the Java Virtual Machine will exit.  If the specified
     * file exists and {@code append} is {@code false}, bytes written by the
     * returned {@code PrintWriter} will overwrite the previously existing file.
     *
     * @param file a file
     * @param append {@code true} if the data will be appended
     * to the end of any existing file
     * @return a non-buffered {@code java.io.PrintWriter} writing to
     * the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter nonBufferedPrintWriter(File file, boolean append) {
        boolean autoflush = true;
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new FileWriter(file, append), autoflush);
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return pw;
    }

    /**
     * Returns a temporary {@code File} that will be deleted when
     * the Java virtual machine exits.
     *
     * @param prefix the filename prefix.
     *
     * @return a {@code File} a new empty file.
     *
     * @throws IllegalArgumentException if {@code prefix} contains fewer than
     * three characters
     */
    public static File tempFile(String prefix) {
        File tempFile = null;
        try {
            tempFile = File.createTempFile(prefix, null);
            tempFile.deleteOnExit();
        } catch (IOException e) {
            Utilities.exit("Exception thrown by createTempFile: ", e);
        }
        return tempFile;
    }
}
