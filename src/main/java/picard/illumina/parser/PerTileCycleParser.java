/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.illumina.parser;

import picard.PicardException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeSet;

/**
 * PerTilePerCycleParser is an abstract IlluminaParser that maintains a list of file parsers for the current tile (1 for each cycle)
 * and coordinates the construction/population of an IlluminaData object on a cycle by cycle basis.
 *
 * @param <ILLUMINA_DATA>
 */
abstract class PerTileCycleParser<ILLUMINA_DATA extends IlluminaData> implements IlluminaParser<ILLUMINA_DATA> {

    /** Location of illumina output files to be parsed.  Typically this is Data/Intensities/L00<lane> */
    private final File laneDirectory;

    /** The lane to iterate over */
    private final int lane;

    /** A parser for the current tile */
    private CycleFilesParser<ILLUMINA_DATA> cycleFilesParser;

    final OutputMapping outputMapping;

    /** The current tile number */
    protected int currentTile;

    /** Map of Cycle -> Tile -> List<File> */
    private final CycleIlluminaFileMap cyclesToTileFiles;

    private final TreeSet<Integer> tileOrder;

    /**
     * Construct a per tile parser
     *
     * @param directory         The directory containing the lane we are analyzing (i.e. the parent of the L00<lane> directory)
     * @param lane              The lane that is being iterated over
     * @param cyclesToTileFiles A map of tile to CycleFilesIterators whose iterators contain only the cycles we want to output
     * @param outputMapping     Data structure containing information on how we should output data
     */
    PerTileCycleParser(final File directory, final int lane, final CycleIlluminaFileMap cyclesToTileFiles, final OutputMapping outputMapping) {
        this.tileOrder = getTileOrder(cyclesToTileFiles);
        this.lane = lane;
        this.laneDirectory = new File(directory, IlluminaFileUtil.longLaneStr(this.lane));
        this.cyclesToTileFiles = cyclesToTileFiles;
        this.currentTile = tileOrder.first();
        this.outputMapping = outputMapping;
    }

    private TreeSet<Integer> getTileOrder(final CycleIlluminaFileMap cyclesToTileFiles) {
        final TreeSet<Integer> uniqueTiles = new TreeSet<Integer>();

        for (final IlluminaFileMap fileMap : cyclesToTileFiles.values()) {
            uniqueTiles.addAll(fileMap.keySet());
        }
        return uniqueTiles;
    }

    /**
     * For a given cycle, return a CycleFilesParser.  It will close the cycleFilesParser if not null.
     *
     * @param file             The file to parse
     * @param cycleFilesParser The previous cycle file parser, null otherwise.
     * @return A CycleFilesParser that will populate the correct position in the IlluminaData object with that cycle's data.
     */
    protected CycleFilesParser<ILLUMINA_DATA> makeCycleFileParser(final List<File> file, final CycleFilesParser<ILLUMINA_DATA> cycleFilesParser) {
        if (cycleFilesParser != null) cycleFilesParser.close();
        return makeCycleFileParser(file);
    }

    protected abstract CycleFilesParser<ILLUMINA_DATA> makeCycleFileParser(final List<File> file);

    public abstract void initialize();

    /**
     * CycleFileParsers iterate through the clusters of a file and populate an IlluminaData object with a single cycle's
     * value.
     *
     * @param <ILLUMINA_DATA>
     */
    protected interface CycleFilesParser<ILLUMINA_DATA> {
        public void close();

        public ILLUMINA_DATA next();

        public boolean hasNext();
    }

    /**
     * Clear the current set of cycleFileParsers and replace them with the ones for the tile indicated by oneBasedTileNumber
     *
     * @param tile requested tile
     */
    @Override
    public void seekToTile(final int tile) {
        currentTile = tile;

        int totalCycles = 0;
        final List<File> tileFiles = new ArrayList<File>();
        for (final Map.Entry<Integer, IlluminaFileMap> entry : cyclesToTileFiles.entrySet()) {
            tileFiles.add(entry.getValue().get(currentTile));
            ++totalCycles;
        }
        cycleFilesParser = makeCycleFileParser(tileFiles, cycleFilesParser);

        if (totalCycles != outputMapping.getTotalOutputCycles()) {
            throw new PicardException("Number of cycle OUTPUT files found (" + totalCycles + ") does not equal the number expected (" + outputMapping.getTotalOutputCycles() + ")");
        }
    }

    /**
     * Return the data for the next cluster by:
     * 1. Advancing tiles if we reached the end of the current tile.
     * 2. For each cycle, get the appropriate parser and have it populate it's data into the IlluminaData object.
     *
     * @return The IlluminaData object for the next cluster
     */
    @Override
    public ILLUMINA_DATA next() { //iterate over clusters
        if (!hasNext()) {
            throw new NoSuchElementException("IlluminaData is missing in lane " + lane + " at directory location " + laneDirectory.getAbsolutePath());
        }

        if (!cycleFilesParser.hasNext()) {
            seekToTile(tileOrder.higher(currentTile));
        }

        return cycleFilesParser.next();
    }

    @Override
    public boolean hasNext() {
        return cycleFilesParser.hasNext() || currentTile < tileOrder.last();
    }

    /**
     * Returns the tile of the next cluster that will be returned by PerTilePerCycleParser and therefore should be called before
     * next() if you want to know the tile for the data returned by next()
     *
     * @return The tile number of the next ILLUMINA_DATA object to be returned
     */
    public int getTileOfNextCluster() {
        //if the current parser still has more clusters, return the current tile number
        if (cycleFilesParser.hasNext()) {
            return currentTile;
        }

        //if the current parser is EMPTY, return the next tile number
        if (currentTile < tileOrder.last()) {
            return tileOrder.higher(currentTile);
        }

        //If we are at the end of clusters then this method should not be called, throw an exception
        throw new NoSuchElementException();
    }

    @Override
    public void verifyData(List<Integer> tiles, final int[] cycles) {
        if (tiles == null) {
            tiles = new ArrayList<Integer>(this.cyclesToTileFiles.keySet());
        }
        this.cyclesToTileFiles.assertValid(tiles, cycles);
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by " + this.getClass().getName());
    }

    @Override
    public void close() {
        cycleFilesParser.close();
    }
}
