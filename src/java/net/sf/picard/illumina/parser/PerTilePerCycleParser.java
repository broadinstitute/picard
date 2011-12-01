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

    /**
     * Computed on construction, since the output data spans multiple arrays this provides an index to the
     * correct array/element for each cycle
     */
    private final CompositeIndex [] cycleToIndex;

    /** The expected lengths of the outputs */
    private final int [] outputLengths;

    /** A sum of the output lengths == totalCycles */
    private final int totalCycles;

    /** The current tile number */
    private int tileNumber  = 0;

    /** Map of tiles -> CycledFilesIterator */
    protected final CycleIlluminaFileMap tilesToCycleFiles;

    public PerTilePerCycleParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles, int[] outputLengths) {
        this.lane = lane;
        this.laneDirectory = new File(directory, "L00" + this.lane);
        this.tilesToCycleFiles = tilesToCycleFiles;
        tileNumber = tilesToCycleFiles.firstKey();
        this.outputLengths = outputLengths;

        int cycles = 0;
        for(int i = 0; i < outputLengths.length; i++) { //TODO: Put this into as util
            cycles += outputLengths[i];
        }
        totalCycles = cycles;
        cycleFileParsers = new ArrayList<CycleFileParser<ILLUMINA_DATA>>(totalCycles);

        cycleToIndex = new CompositeIndex[cycles+1];
        cycleToIndex[0] = null;
        int arrIndex = 0;
        int elementIndex = 0;
        for(int i = 0; i < totalCycles; i++) {
            if(elementIndex >= outputLengths[arrIndex]) {
                elementIndex = 0;
                ++arrIndex;
            }
            cycleToIndex[i+1] = new CompositeIndex(arrIndex, elementIndex);
            ++elementIndex;
        }
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
     * @return A CycleFileParser that will populate the correct position in the IlluminaData object with that cycles data.
     */
    protected abstract CycleFileParser<ILLUMINA_DATA> makeCycleFileParser(final File file, final int cycle);

    /**
     * CycleFileParsers iterate through the clusters of a file and populate an IlluminaData object with a single cycle's
     * value.
     * @param <ILLUMINA_DATA>
     */
    protected interface CycleFileParser<ILLUMINA_DATA>  { //Alternatively since we seem to be using MemoryMappedByteBuffers and these files have a numClusters
        public void close();                              //we could just have CycleFileParser only have 1 method addData(ild, cluster)
        public void next(final ILLUMINA_DATA ild);
        public boolean hasNext();
    }

    /**
     * Return the precomputed output array and element indices of the given cycle
     * @param cycle requested cycle
     * @return A composite index with the correct array/element that will locate that cycle in output data arrays/multi-element objects
     */
    protected CompositeIndex getIndex(final int cycle) {
        return cycleToIndex[cycle];
    }

    /**
     * Clear the current set of cycleFileParsers and replace them with the ones for the tile indicated by oneBasedTileNumber
     * @param oneBasedTileNumber requested tile with indices beginning at 1
     */
    @Override
    public void seekToTile(final int oneBasedTileNumber) {
        tileNumber = oneBasedTileNumber;
        final CycleFilesIterator filesIterator = tilesToCycleFiles.get(tileNumber);

        CloserUtil.close(cycleFileParsers);
        cycleFileParsers.clear();

        int cycleIndex = 0;
        while(filesIterator.hasNext()) {
            cycleFileParsers.add(makeCycleFileParser(filesIterator.next(), ++cycleIndex));
        }

        if(cycleIndex != totalCycles) {
            throw new PicardException("Number of cycle files found (" + cycleIndex + ") does not equal the number expected (" + totalCycles +")");
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
            seekToTile(tileNumber+1);
        }

        final ILLUMINA_DATA data = makeData(outputLengths);
        for(int i = 0; i < totalCycles; i++) {
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

    @Override
    public void verifyData(final IlluminaRunConfiguration runConfig, List<Integer> tiles) {
        if(tiles == null) {
            tiles = new ArrayList<Integer>(this.tilesToCycleFiles.keySet());
        }
        this.tilesToCycleFiles.assertValid(tiles, runConfig.totalCycles);
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by " + this.getClass().getName());
    }
}
