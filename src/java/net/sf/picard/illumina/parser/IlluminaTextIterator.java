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

import net.sf.picard.util.BasicInputParser;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import net.sf.picard.PicardException;
import net.sf.samtools.util.CloserUtil;

/**
 * Abstract class for parsing text-based whitespace-delimited Illumina output files, organized
 * by tile.  Concrete implementations must call setFiles() in order to provide the list of files
 * to be iterated over.
 * 
 * @author jburke@broadinstitute.org
 */
class IlluminaTextIterator implements Iterator<String[]> {
    // Location of illumina output files to be parsed
    private final int lane;
    private int currentTile = 0;

    // List of files of the given type, sorted by tile #
    private IlluminaFileMap files;

    private boolean treatGroupedDelimitersAsOne = true;
    private BasicInputParser parser;

    public IlluminaTextIterator(final int lane, final IlluminaFileMap files) {
        this.lane = lane;
        this.files = files;
        currentTile = files.firstKey();
    }

    public IlluminaTextIterator(final int lane, final IlluminaFileMap files,
                                final boolean treatGroupedDelimitersAsOne) {
        this.lane = lane;
        this.files = files;
        this.treatGroupedDelimitersAsOne = treatGroupedDelimitersAsOne;
        currentTile = files.firstKey();
    }
    
    /**
     * Jump so that the next record returned will be the first one from the specified tile.
     */
    public void seekToTile(final int oneBasedTileNumber) {
        CloserUtil.close(parser);
        currentTile = oneBasedTileNumber;
        initializeParser();
    }

    /**
     * Prepare to iterate.
     */
    private void initializeParser() {
        final List<File> fileSubset = files.getFilesStartingAt(currentTile);
        parser = new BasicInputParser(treatGroupedDelimitersAsOne, fileSubset.toArray(new File[fileSubset.size()]));
    }

    /**
     * Read the next record from the list of input files, and load into data argument.
     */
    @Override
    public String [] next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }

        return parser.next();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by IlluminaTextIterator");
    }

    public boolean hasNext() {
        if(parser == null) initializeParser();
        return parser.hasNext();
    }

    protected int getLane() {
        return lane;
    }

    public String getCurrentFilename() {
        if(parser == null) initializeParser();
        return parser.getFileName();
    }

    protected void validateLane(final int lane) {
        if (lane != getLane()) {
            throw new PicardException("Lane number mismatch: " + lane + " != " + getLane());
        }
    }
}
