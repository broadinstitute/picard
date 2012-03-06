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
package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;
import net.sf.samtools.util.CloserUtil;

import java.util.*;
import java.io.File;

/**
 * PerTilePerCycleParser is an abstract IlluminaParser that maintains a list of file parsers for the current tile (1 for each cycle)
 * and coordinates the construction/population of an IlluminaData object on a cycle by cycle basis.
 * @param <ILLUMINA_DATA>
 */
abstract class PerTilePerCycleParser<ILLUMINA_DATA extends IlluminaData> implements IlluminaParser<ILLUMINA_DATA> {

    /** Location of illumina output files to be parsed.  Typically this is Data/Intensities/L00<lane> */
    private final File laneDirectory;

    /** The lane to iterate over */
    private final int lane;

    /** A set of parsers for the current tile ordered by cycle */
    private final List<CycleFileParser<ILLUMINA_DATA>> cycleFileParsers;

    protected final OutputMapping outputMapping;

    /** The current tile number */
    private int tileNumber;

    /** Map of tiles -> CycledFilesIterator, the CycleFileIterators of this map should contain ONLY cycles to be output */
    protected final CycleIlluminaFileMap tilesToCycleFiles;

    /**
     * Construct a per tile parser
     * @param directory The directory containing the lane we are analyzing (i.e. the parent of the L00<lane> directory)
     * @param lane The lane that is being iterated over
     * @param tilesToCycleFiles A map of tile to CycleFilesIterators whose iterators contain only the cycles we want to output
     * @param outputMapping Data structure containing information on how we should output data
     */
    public PerTilePerCycleParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, final OutputMapping outputMapping) {
        this.lane = lane;
        this.laneDirectory = new File(directory, "L00" + this.lane);
        this.tilesToCycleFiles = tilesToCycleFiles;
        this.tileNumber = tilesToCycleFiles.firstKey();
        this.outputMapping = outputMapping;

        cycleFileParsers = new ArrayList<CycleFileParser<ILLUMINA_DATA>>(outputMapping.getTotalOutputCycles());

        seekToTile(tilesToCycleFiles.firstKey());
    }

    /**
     * Per cluster makeData will make the relevant IlluminaData object with the given outputLengths
     * @param outputLengths The expected lengths of the output data
     * @return An IlluminaData object that has been constructed with default(according to Java construction rules) data
     */
    protected abstract ILLUMINA_DATA makeData(final int [] outputLengths);

    /**
     * For a given cycle, return a CycleFileParser.
     * @param file The file to parse
     * @param cycle The cycle that file represents
     * @return A CycleFileParser that will populate the correct position in the IlluminaData object with that cycle's data.
     */
    protected abstract CycleFileParser<ILLUMINA_DATA> makeCycleFileParser(final File file, final int cycle);

    /**
     * CycleFileParsers iterate through the clusters of a file and populate an IlluminaData object with a single cycle's
     * value.
     * @param <ILLUMINA_DATA>
     */
    protected interface CycleFileParser<ILLUMINA_DATA>  {
        public void close();
        public void next(final ILLUMINA_DATA ild);
        public boolean hasNext();
    }

    /**
     * Clear the current set of cycleFileParsers and replace them with the ones for the tile indicated by oneBasedTileNumber
     * @param oneBasedTileNumber requested tile with indices beginning at 1
     */
    @Override
    public void seekToTile(final int oneBasedTileNumber) {
        tileNumber = oneBasedTileNumber;
        final CycleFilesIterator filesIterator = tilesToCycleFiles.get(tileNumber);
        filesIterator.reset();

        CloserUtil.close(cycleFileParsers);
        cycleFileParsers.clear();

        int totalCycles = 0;
        while(filesIterator.hasNext()) {
            final int nextCycle = filesIterator.getNextCycle();
            cycleFileParsers.add(makeCycleFileParser(filesIterator.next(), nextCycle));
            ++totalCycles;
        }

        if(totalCycles != outputMapping.getTotalOutputCycles()) {
            throw new PicardException("Number of cycle OUTPUT files found (" + totalCycles + ") does not equal the number expected (" + outputMapping.getTotalOutputCycles() +")");
        }
    }

    /**
     *  Return the data for the next cluster by:
     *  1. Advancing tiles if we reached the end of the current tile.
     *  2. For each cycle, get the appropriate parser and have it populate it's data into the IlluminaData object.
     * @return The IlluminaData object for the next cluster
     */
    @Override
    public ILLUMINA_DATA next() { //iterate over clusters
        if(!hasNext()) {
            throw new NoSuchElementException("IlluminaData is missing in lane " + lane + " at directory location " + laneDirectory.getAbsolutePath());
        }

        if(!cycleFileParsers.get(0).hasNext()) {
            seekToTile(tilesToCycleFiles.higherKey(tileNumber));
        }

        final ILLUMINA_DATA data = makeData(outputMapping.getOutputReadLengths());
        for(int i = 0; i < outputMapping.getTotalOutputCycles(); i++) {
            cycleFileParsers.get(i).next(data);
        }

        return data;
    }

    @Override
    public boolean hasNext() {
        if(cycleFileParsers.get(0).hasNext()) {
            return true;
        }

        return tileNumber < tilesToCycleFiles.lastKey();
    }

    /**
     * Returns the tile of the next cluster that will be returned by PerTilePerCycleParser and therefore should be called before
     * next() if you want to know the tile for the data returned by next()
     * @return The tile number of the next ILLUMINA_DATA object to be returned
     */
    public int getTileOfNextCluster() {
        //if the current parser still has more clusters, return the current tile number
        if(cycleFileParsers.get(0).hasNext()) {
            return tileNumber;
        }

        //if the current parser is EMPTY, return the next tile number
        if(tileNumber < tilesToCycleFiles.lastKey()) {
            return tilesToCycleFiles.higherKey(tileNumber);
        }

        //If we are at the end of clusters then this method should not be called, throw an exception
        throw new NoSuchElementException();
    }

    @Override
    public void verifyData(List<Integer> tiles, final int [] cycles) {
        if(tiles == null) {
            tiles = new ArrayList<Integer>(this.tilesToCycleFiles.keySet());
        }
        this.tilesToCycleFiles.assertValid(tiles, cycles);
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by " + this.getClass().getName());
    }
}
