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
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.*;

/**
 * A sorted map of tiles to Iterators, where each iterator traverses the cycle files for the given lane/tile.
 * @author jburke@broadinstitute.org
 */
class CycleIlluminaFileMap extends TreeMap<Integer, CycleFilesIterator> {
    /** Return a CycleIlluminaFileMap with only the tiles listed and all of the cycles provided, throws an exception if it does not contain any of the tiles listed
     * Important NOTE: this DOES NOT eliminate cycles from the cycles parameter passed in that are missing in the cyclesFileIterator of any given lane in the CycleIlluminaFileMap
     * */
    public CycleIlluminaFileMap keep(List<Integer> tilesToKeep, final int [] cycles) {
        if(tilesToKeep == null) {
            tilesToKeep = new ArrayList<Integer>(this.keySet());
        }

        final CycleIlluminaFileMap ciMap = new CycleIlluminaFileMap();
        for(final Integer tile : tilesToKeep) {
            final CycleFilesIterator template = this.get(tile);
            if(template != null) {
                ciMap.put(tile, new CycleFilesIterator(this.get(tile), cycles));
            }
        }
        return ciMap;
    }

    /**
     * Assert that this map has an iterator for all of the expectedTiles and each iterator has expectedCycles number
     * of files.  Also, assert that each cycle file for a given tile is the same size
     * @param expectedTiles A list of tiles that should be in this map
     * @param expectedCycles The total number of files(cycles) that should be in each CycledFilesIterator
     */
    public void assertValid(final List<Integer> expectedTiles, final int [] expectedCycles) {
        if(size() != expectedTiles.size()) {
            throw new PicardException("Expected CycledIlluminaFileMap to contain " + expectedTiles + " tiles but only " + size() + " were found!");
        }

        File curFile = null;

        for (final int tile : expectedTiles) {
            final CycleFilesIterator cycleFiles = new CycleFilesIterator(get(tile), null); //so as not to expend the iterator
            int total;
            for(total = 0; cycleFiles.hasNext(); total++) {
                if(cycleFiles.getNextCycle() != expectedCycles[total]) {
                    if(curFile == null) {
                        curFile = cycleFiles.next();
                    }
                    cycleFiles.reset();
                    throw new PicardException("Cycles in iterator(" + remainingCyclesToString(cycleFiles) +
                                              ") do not match those expected (" + StringUtil.intValuesToString(expectedCycles) +
                                              ") Last file checked (" + curFile.getAbsolutePath() + ")");
                }
                curFile = cycleFiles.next();
                if(!curFile.exists()) {
                    throw new PicardException("Missing cycle file " + curFile.getName() + " in CycledIlluminaFileMap");
                }
            }
            if(total != expectedCycles.length) {
                String message = "Expected tile " + tile + " of CycledIlluminaFileMap to contain " + expectedCycles + " cycles but " + total + " were found!";

                if(curFile != null) {
                    message += "Check to see if the following bcl exists: " + incrementCycleCount(curFile).getAbsolutePath();
                }
                throw new PicardException(message);
            }
        }
    }

    /**
     * All files in a CycleFileMap should have a cycle directory (at this point in time). Return the filePath with the cycle
     * in the cycle directory incremented by 1.
     * @param cycleFile The cycleFile for a given tile/cycle
     * @return A cycleFile with the same name but a cycleDir one greater than cycleFile
     */
    public static File incrementCycleCount(final File cycleFile) {
        final File cycleDir = cycleFile.getParentFile();
        final int cycle = Integer.parseInt(cycleDir.getName().substring(1, cycleDir.getName().lastIndexOf(".")));
        return new File(cycleDir.getParentFile(), "C" + cycle + ".1" + File.separator + cycleFile.getName());
    }

    public static String remainingCyclesToString(final CycleFilesIterator cfi) {
        String cycles = "";
        if(cfi.hasNext()) {
            cycles += cfi.getNextCycle();
            cfi.next();
        }

        while(cfi.hasNext()) {
            cycles += ", " + cfi.getNextCycle();
            cfi.next();
        }

        return cycles;
    }
}

/**
 * Given a lane directory, lane, tile number, and file extension get provide an iterator over files in the following order
 * <LaneDir>/<nextCycle>.1/s_<lane>_<tile>.<fileExt>
 * In other words iterate through the different cycle directory and find the file for the current cycle with the
 * given file extension while lane/tile stay the same.  Stop if the next file does not exist.
 */
 class CycleFilesIterator implements Iterator<File>, Iterable<File> {
    private final File parentDir;
    private final int lane;
    private final int tile;
    private final String fileExt;
    protected final int [] cycles;
    protected int nextCycleIndex;

    public CycleFilesIterator(final File laneDir, final int lane, final int tile, final int [] cycles, final String fileExt) {
        this.parentDir = laneDir;
        this.lane = lane;
        this.tile = tile;
        this.fileExt = fileExt;
        this.cycles  = cycles;
        this.nextCycleIndex = 0;
    }

    CycleFilesIterator(final CycleFilesIterator template, final int [] cycles) {
        this(template.parentDir, template.lane, template.tile, (cycles != null) ? cycles : template.cycles, template.fileExt);
    }

    public void reset() {
        this.nextCycleIndex = 0;
    }

    @Override
    public boolean hasNext() {
        return nextCycleIndex < cycles.length;
    }

    @Override
    public File next() {
        if (!hasNext()) {
            throw new NoSuchElementException(summarizeIterator());
        }

        final File cycleDir = new File(parentDir, "C" + cycles[nextCycleIndex] + ".1");
        final File curFile  = new File(cycleDir, "s_" + lane + "_" + tile + fileExt);

        nextCycleIndex++;
        return curFile;
    }

    public int getNextCycle() {
        return cycles[nextCycleIndex];
    }

    private String summarizeIterator() {
        String cyclesSummary = "";
        if(cycles.length > 0) {
            cyclesSummary = String.valueOf(cycles[0]);
            for(int i = 1; i < cycles.length; i++) {
                cyclesSummary += ", " + String.valueOf(cycles[i]);
            }
        }

        return " Parent dir (" + parentDir.getAbsolutePath() + ")" +
                " Lane (" + lane + ") Tile (" + tile + ") FileExt(" + fileExt + ")" +
                " CycleIndex (" + nextCycleIndex + ")" + "Cycles(" + cyclesSummary + ")";
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by " + CycleFilesIterator.class.getName());
    }

    @Override
    public Iterator<File> iterator() {
        return this;
    }
}
