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
package org.ngseq.sparkbeagle;

import org.apache.hadoop.fs.FSDataInputStream;
import org.ngseq.sparkbeagle.beagleutil.Samples;
import org.ngseq.sparkbeagle.blbutil.FileIt;
import org.ngseq.sparkbeagle.blbutil.InputIt;
import org.ngseq.sparkbeagle.blbutil.StringUtil;
import org.ngseq.sparkbeagle.ints.IntList;

import java.util.Arrays;

/**
 * <p>Class {@code Pedigree} stores parent-offspring relationships
 * in a list of samples.  In particular, class {@code Pedigree}
 * stores a list of the single individuals in the list of samples,
 * a list of the parent-offspring duos in the list of samples, and a list of
 * the parent-offspring trios in the list of samples. A single individual is
 * an  individuals without a parent or offspring in the list of samples.
 * </p>
 * <p>Instances of class {@code Pedigree} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Pedigree {

    private static final String NO_PAR = "0";   // NO_PARENT code

    private final Samples samples;

    private final int[] singles;
    private final int[] relateds;
    private final int[] duoOffspring;
    private final int[] trioOffspring;

    private final int[] mothers;
    private final int[] fathers;
    private final int[][] offspring;

    /**
     * Constructs a new {@code NuclearFamilies} instance.
     *
     * @param samples the list of samples.
     * @param pedFile a linkage-format pedigree file, or {@code null}
     * if no pedigree relationships are known.  A pedigree file must have
     * at least 4 white-space delimited columns.  The first column of the
     * pedigree file (family ID) is ignored.  The second, third, and fourth
     * columns are the individual's ID, the individual's father's ID, and
     * the individual's mother's ID respectively.
     *
     * @throws IllegalArgumentException if a pedigree file is specified,
     * and if the file has a non-blank line with less than 4 white-space
     * delimited fields
     * @throws IllegalArgumentException if a pedigree file is specified,
     * and if the file has duplicate individual identifiers in the
     * second white-space delimited column
     * @throws NullPointerException if {@code samples == null}
     */
    public Pedigree(Samples samples, FSDataInputStream pedFile) {
        int nSamples = samples.nSamples();
        this.samples = samples;
        this.fathers = new int[nSamples];
        this.mothers = new int[nSamples];
        this.offspring = new int[nSamples][];
        Arrays.fill(fathers, -1);
        Arrays.fill(mothers, -1);
        Arrays.fill(offspring, new int[0]);

        readPedFile(samples, pedFile, fathers, mothers, offspring);

        int[] cnts = counts(fathers, mothers, offspring);
        this.singles = new int[cnts[0]];
        this.duoOffspring = new int[cnts[1]];
        this.trioOffspring = new int[cnts[2]];
        this.relateds = new int[nSamples - cnts[0]];
        fillArrays(samples, fathers, mothers, offspring, singles, duoOffspring,
                trioOffspring, relateds);
    }

    private static void readPedFile(Samples samples, FSDataInputStream pedFile,
            int[] fathers, int[] mothers, int[][] offspring) {
        if (pedFile != null) {
            boolean[] processed = new boolean[samples.nSamples()];
            IntList[] children = new IntList[samples.nSamples()];
            try (FileIt<String> pedIt=InputIt.fromGzipFile(pedFile, false)) {
                while (pedIt.hasNext()) {
                    String line = pedIt.next().trim();
                    readPedLine(samples, processed, line, fathers, mothers,
                            children);
                }
            }
            storeOffspring(children, offspring);
        }
    }

    private static void readPedLine(Samples samples, boolean[] processed,
            String line, int[] fathers, int[] mothers, IntList[] children) {
        if (line.length() > 0) {
            String[] sa = getPedFields(line);
            String childId = sa[1];
            int child = samples.index(childId);
            if (child != -1) {
                markAsProcessed(processed, child, childId);
                int father = sa[2].equals(NO_PAR) ? -1 : samples.index(sa[2]);
                int mother = sa[3].equals(NO_PAR) ? -1 : samples.index(sa[3]);
                setParentOffspring(father, child, fathers, children);
                setParentOffspring(mother, child, mothers, children);
            }
        }
    }

    private static String[] getPedFields(String line) {
        String[] fields = StringUtil.getFields(line, 5);
        if (fields.length < 4) {
            String s = "invalid line in ped file: " + line;
            throw new IllegalArgumentException(s);
        }
        return fields;
    }

    private static void markAsProcessed(boolean[] processed, int index,
            String id) {
        if (processed[index]) {
            String s = "duplicate sample in pedigree file: "
                    + id;
            throw new IllegalArgumentException(s);
        }
        else {
            processed[index] = true;
        }
    }

    private static void setParentOffspring(int parent, int child,
            int[] sample2Parent,  IntList[] children) {
        if (parent != -1) {
            sample2Parent[child] = parent;
            if (children[parent]==null) {
                children[parent] = new IntList(3);
            }
            children[parent].add(child);
        }
    }

    private static void storeOffspring(IntList[] children, int[][] sample2Offspring) {
        for (int j=0, n=children.length; j<n; ++j) {
            IntList il = children[j];
            if (il != null) {
                int[] ia = il.toArray();
                Arrays.sort(ia);
                sample2Offspring[j] = ia;
            }
        }
    }

    private int[] counts(int[] fathers, int[] mothers, int[][] offspring) {
        assert fathers.length==mothers.length;
        assert fathers.length==offspring.length;
        int nSingles = 0;
        int nDuos = 0;
        int nTrios = 0;
        for (int s=0; s<offspring.length; ++s) {
            int nParents = nParents(s, fathers, mothers);
            switch (nParents) {
                case 0 :
                    if (offspring[s].length==0) {
                        ++nSingles;
                    }
                    break;
                case 1 :
                    ++nDuos;
                    break;
                case 2:
                    ++nTrios;
                    break;
                default:
                    assert false;
            }
        }
        return new int[] {nSingles, nDuos, nTrios};
    }

    private static void fillArrays(Samples samples,
            int[] father, int[] mother, int[][] offspring, int[] single,
            int[] duoOffspring, int[] trioOffspring, int[] relateds) {
        int singleIndex = 0;
        int duoIndex = 0;
        int trioIndex = 0;
        int relatedIndex = 0;
        for (int s=0, n=samples.nSamples(); s<n; ++s) {
            int nParents = nParents(s, father, mother);
            switch (nParents) {
                case 0:
                    if (offspring[s].length==0) {
                        single[singleIndex++] = s;
                    }
                    else {
                        relateds[relatedIndex++] = s;
                    }
                    break;
                case 1:
                    duoOffspring[duoIndex++] = s;
                    relateds[relatedIndex++] = s;
                    break;
                case 2:
                    trioOffspring[trioIndex++] = s;
                    relateds[relatedIndex++] = s;
                    break;
                default:
                    assert false;
            }
        }
        assert singleIndex==single.length;
        assert duoIndex==duoOffspring.length;
        assert trioIndex==trioOffspring.length;
        assert relatedIndex==relateds.length;
    }

    private static int nParents(int index, int[] father, int[] mother) {
        int cnt = 0;
        if (father[index]>=0) {
            ++cnt;
        }
        if (mother[index]>=0) {
            ++cnt;
        }
        return cnt;
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return samples.nSamples();
    }

    /**
     * Returns the number of single individuals in the list of samples.
     * A single individual has no parent or offspring in the list of samples.
     * @return the number of single individuals in the sample
     */
    public int nSingles() {
        return singles.length;
    }

    /**
     * Returns the number of parent-offspring duos in the list of samples.
     * The offspring of a parent-offspring duo has only one parent
     * in the sample.
     * @return the number of parent-offspring duos in the list of samples
     */
    public int nDuos() {
        return duoOffspring.length;
    }

    /**
     * Returns the number of parent-offspring trios in the list of samples.
     * The offspring of a parent-offspring trio has two parents
     * in the sample.
     * @return the number of parent-offspring trios in the list of samples
     */
    public int nTrios() {
        return trioOffspring.length;
    }

    /**
     * Returns an array of indices of samples with no parents or
     * children in the list of samples.
     * @return an array of indices of samples with no parents or
     * children in the list of samples
     */
    public int[] singles() {
        return singles.clone();
    }

    /**
     * Returns an array of indices of samples with at least one
     * parent or child in the list or samples.
     * @return an array of indices of samples with at least one
     * parent or child in the list or samples
     */
    public int[] relateds() {
        return relateds.clone();
    }

    /**
     * Returns the sample index of the specified single individual.
     * A single individual has no first-degree relative in the list of
     * samples.
     * @param index the index of a single individual
     * @return the sample index of the specified single individual
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nSingles()}
     */
    public int single(int index) {
        return singles[index];
    }

    /**
     * Returns the sample index of the parent of the specified
     * parent-offspring duo.
     * @param index the index of a parent-offspring duo
     * @return the sample index of the parent of the specified
     * parent-offspring duo
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nDuos()}
     */
    public int duoParent(int index) {
        int offspring = duoOffspring[index];
        if (fathers[offspring]>=0) {
            return fathers[offspring];
        }
        else {
            assert mothers[offspring]>=0;
            return mothers[offspring];
        }
    }

    /**
     * Returns the sample index of the offspring of the specified
     * parent-offspring duo.
     * @param index the index of a parent-offspring duo
     * @return the sample index of the offspring of the specified
     * parent-offspring duo
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nDuos()}
     */
    public int duoOffspring(int index) {
        return duoOffspring[index];
    }

    /**
     * Returns the sample index of the father of the specified
     * parent-offspring trio.
     * @param index the index of a parent-offspring trio
     * @return the sample index of the father of the specified
     * parent-offspring trio
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nTrios()}
     */
    public int trioFather(int index) {
        return fathers[trioOffspring[index]];
    }

    /**
     * Returns the sample index of the mother of the specified
     * parent-offspring trio.
     * @param index the index of a parent-offspring trio
     * @return the sample index of the mother of the specified
     * parent-offspring trio
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nTrios()}
     */
    public int trioMother(int index) {
        return mothers[trioOffspring[index]];
    }

    /**
     * Returns the sample index of the offspring of the specified
     * parent-offspring trio.
     * @param index the index of a parent-offspring trio
     * @return the sample index of the offspring of the specified
     * parent-offspring trio
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nTrios()}
     */
    public int trioOffspring(int index) {
        return trioOffspring[index];
    }

    /**
     * Returns the sample index of the father of the specified sample,
     * or returns {@code -1} if the father is unknown or is not present
     * in the list of samples.
     * @param sample a sample index
     * @return the sample index of the father of the specified sample,
     * or {@code -1} if the father is unknown or is not present in
     * the list of samples
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()()}
     */
    public int father(int sample) {
        return fathers[sample];
    }

    /**
     * Returns the sample index of the mother of the specified sample,
     * or returns {@code -1} if the mother is unknown or is not present
     * in the list of samples.
     * @param sample a sample index
     * @return the sample index of the mother of the specified sample,
     * or {@code -1} if the mother is unknown or is not present
     * in the list of samples
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()()}
     */
    public int mother(int sample) {
        return mothers[sample];
    }

    /**
     * Returns the number of offspring of the specified sample.
     * @param sample a sample index
     * @return the number of offspring of the specified sample
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()()}
     */
    public int nOffspring(int sample) {
        return offspring[sample].length;
    }

    /**
     * Returns the sample index of the offspring of the specified sample.
     * @param sample a sample index
     * @param index the offspring index
     * @return the sample index of the offspring of the specified sample.
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nOffspring(sample)}
     */
    public int offspring(int sample, int index) {
        return offspring[sample][index];
    }

    /**
     * Returns a string representation of {@code this}.  The exact details of
     * the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return this.getClass().toString();
    }
}
