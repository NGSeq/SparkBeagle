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
package org.ngseq.sparkbeagle.bref;

import org.ngseq.sparkbeagle.beagleutil.ChromIds;
import org.ngseq.sparkbeagle.blbutil.Const;

/**
 * <p>Class {@code BrefBlock} represents starting chromosome coordinates and
 * file offset for the start of a binary reference format (bref) data block.
 * </p>
 *
 * <p>Instances of class {@code BrefBlock} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class BrefBlock {

    private final int chromIndex;
    private final int pos;
    private final long offset;

    /**
     * Constructs a {@code BrefBlock} for the specified data.  It is the
     * caller's responsibility to ensure the consistency of the
     * constructor parameters.
     *
     * @param chromIndex the chromosome index
     * @param pos the starting chromosome position
     * @param offset the file offset in bytes for the bref data block
     */
    public BrefBlock(int chromIndex, int pos, long offset) {
        this.chromIndex = chromIndex;
        this.pos = pos;
        this.offset = offset;
    }

    /**
     * Returns the chromosome index of the first marker in this bref block.
     * @return the chromosome index of the first marker in this bref block
     */
    public int chromIndex() {
        return chromIndex;
    }

    /**
     * Returns the chromosome position of the first marker in this bref block.
     * @return the chromosome position of the first marker in this bref block
     */
    public int pos() {
        return pos;
    }

    /**
     * Returns the file offset of the first marker in this bref block.
     * @return the file offset of the first marker in this bref block
     */
    public long offset() {
        return offset;
    }

    /**
     * Returns a string description of {@code this}.  The exact details
     * of the representation are unspecified and subject to change.
     * @return  a string description of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(30);
        sb.append('[');
        sb.append(ChromIds.instance().id(chromIndex));
        sb.append(Const.tab);
        sb.append(pos);
        sb.append(Const.tab);
        sb.append(offset);
        sb.append(']');
        return sb.toString();
    }
}
