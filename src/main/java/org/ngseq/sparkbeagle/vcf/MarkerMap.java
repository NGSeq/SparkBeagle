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

/**
 * <p>Class {@code MarkerMap} represents genetic map positions for a
 * list of markers.
 * </p>
 * <p>Instances of class {@code MarkerMap} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MarkerMap implements GeneticMap {

    private final GeneticMap genMap;
    private final Markers markers;
    private final double[] genPos;

    /**
     * Returns a new {@code MarkerMap} instance that is constructed from the
     * specified data
     * @param genMap the genetic map
     * @param markers a list of markers
     * @return a returns new {@code MarkerMap} instance
     * @throws IllegalArgumentException if
     * {@code markers.marker(0).chromIndex()
     *        != markers.marker(markers.nMarkers()-1).chromIndex()}
     * @throws NullPointerException if {@code genMap == null || markers == null}
     */
    public static MarkerMap create(GeneticMap genMap, Markers markers) {
        if (markers.marker(0).chromIndex()
                != markers.marker(markers.nMarkers()-1).chromIndex()) {
            throw new IllegalArgumentException("multiple chromosomes");
        }
        return new MarkerMap(genMap, markers);
    }

    /**
     * Constructs a new {@code MarkerMap} instance from the specified data
     * @param genMap the genetic map
     * @param markers a list of markers
     * @throws IllegalArgumentException if {@code nHaps < 1}
     * @throws IllegalArgumentException if {@code ne < 1.0}
     * @throws NullPointerException if {@code genMap == null || markers == null}
     */
    private MarkerMap(GeneticMap genMap, Markers markers) {
        this.genMap = genMap;
        this.markers = markers;
        this.genPos = genMap.genPos(markers);
    }

    /**
     * Returns the chromosome corresponding to the list of markers.
     * @return the chromosome corresponding to the list of markers
     */
    public String chrom() {
        return markers.marker(0).chrom();
    }

    /**
     * Returns the index of the chromosome corresponding to the list of markers.
     * @return the index of the chromosome corresponding to the list of markers
     */
    public int chromIndex() {
        return markers.marker(0).chromIndex();
    }

    @Override
    public int basePos(int chrom, double geneticPosition) {
        return genMap.basePos(chrom, geneticPosition);
    }

    @Override
    public double index2GenPos(int chrom, int index) {
        return genMap.index2GenPos(chrom, index);
    }

    @Override
    public int nMapPositions(int chrom) {
        return genMap.nMapPositions(chrom);
    }

    @Override
    public double genPos(Marker marker) {
        return genMap.genPos(marker);
    }

    @Override
    public double genPos(int chrom, int basePosition) {
        return genMap.genPos(chrom, basePosition);
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return markers.nMarkers();
    }

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    public Markers markers() {
        return markers;
    }

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     */
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    /**
     * Returns the genetic map position of the specified marker.
     * @param marker a marker index
     * @return the genetic map position of the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public double genPos(int marker) {
        return genPos[marker];
    }

    /**
     * Returns the array of genetic map positions whose {@code k}-th element
     * equals {@code this.genPos(k)}.
     * @return the array of genetic map positions
     */
    public double[] genPos() {
        return genPos.clone();
    }
}
