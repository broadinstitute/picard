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

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.TreeMap;

/**
 * A sorted map of tiles to Iterators, where each iterator traverses the cycle files for the given lane/tile.
 * @author jburke@broadinstitute.org
 */
class CycleIlluminaFileMap extends TreeMap<Integer, CycleFilesIterator> {

    /** Return a CycleIlluminaFileMap with only the tiles listed , throws an exception if it does not contain any of the tiles listed */
    public CycleIlluminaFileMap keep(final List<Integer> tilesToKeep) {
        final CycleIlluminaFileMap ciMap = new CycleIlluminaFileMap();
        for(final Integer tile : tilesToKeep) {
            final CycleFilesIterator template = this.get(tile);
            if(template != null) {
                ciMap.put(tile, new CycleFilesIterator(this.get(tile)));
            }
        }
        return ciMap;
    }

    /**
     * Assert that this map has an iterator for all of the expectedTiles and each iterator has expectedCycles number
     * of files.
     * @param expectedTiles A list of tiles that should be in this map
     * @param expectedCycles The total number of files(cycles) that should be in each CycledFilesIterator
     */
    public void assertValid(final List<Integer> expectedTiles, final int expectedCycles) {
        if(size() != expectedTiles.size()) {
            throw new PicardException("Expected CycledIlluminaFileMap to contain " + expectedTiles + " tiles but only " + size() + " were found!");
        }

        for (final int tile : expectedTiles) {
            final CycleFilesIterator cycleFiles = new CycleFilesIterator(get(tile)); //so as not to expend the iterator
            int total;
            for(total = 0; cycleFiles.hasNext(); total++) {
                final File curFile = cycleFiles.next();
                if(!curFile.exists()) {
                    throw new PicardException("Missing cycle file " + curFile.getName() + " in CycledIlluminaFileMap");
                }
            }
            if(total != expectedCycles) {
                throw new PicardException("Expected tile " + tile + " of CycledIlluminaFileMap to contain " + expectedCycles + " cycles but " + total + " were found!");
            }
        }
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

    private File nextFile;
    private int nextCycle;

    public CycleFilesIterator(final File laneDir, final int lane, final int tile, final String fileExt) {
        this.parentDir = laneDir;
        this.lane = lane;
        this.tile = tile;
        this.fileExt = fileExt;
        this.nextCycle = 1;
        findNextFile();
    }

    CycleFilesIterator(final CycleFilesIterator template) {
        this(template.parentDir, template.lane, template.tile, template.fileExt);
    }

    private void findNextFile() {
        final File cycleDir = new File(parentDir, "C" + nextCycle + ".1");
        nextFile = new File(cycleDir, "s_" + lane + "_" + tile + fileExt);
    }

    public void reset() {
        this.nextCycle = 1;
        findNextFile();
    }

    @Override
    public boolean hasNext() {
        return nextFile.exists();
    }

    @Override
    public File next() {
        if (!hasNext()) {
            throw new NoSuchElementException( " Parent dir (" + parentDir.getAbsolutePath() + ")" +
                                              " Lane (" + lane + ") Tile (" + tile + ") FileExt(" + fileExt + ")" +
                                              " Cycle (" + nextCycle + ")");
        }
        final File curFile = nextFile;
        nextCycle++;
        findNextFile();
        return curFile;
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
