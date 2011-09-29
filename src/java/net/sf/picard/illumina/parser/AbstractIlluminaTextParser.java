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
import java.util.NoSuchElementException;

import net.sf.picard.util.FormatUtil;
import net.sf.picard.PicardException;
import net.sf.picard.util.Log;

/**
 * Abstract class for parsing text-based whitespace-delimited Illumina output files, organized
 * by tile.  Concrete implementations must call setFiles() in order to provide the list of files
 * to be iterated over.
 * 
 * @author alecw@broadinstitute.org
 */
abstract class AbstractIlluminaTextParser implements IlluminaParser {
    private static final FormatUtil formatter = new FormatUtil();

    // Describe relationship between cycle numbers and read ends and barcode
    private final ReadConfiguration readConfiguration;

    // Location of illumina output files to be parsed
    private final File directory;
    private final int lane;

    // List of files of the given type, sorted by tile #
    private TiledIlluminaFile[] files;

    private boolean treatGroupedDelimitersAsOne = true;

    private BasicInputParser parser;
    private final Log log = Log.getInstance(AbstractIlluminaTextParser.class);

    public AbstractIlluminaTextParser(final ReadConfiguration readConfiguration, final int lane, final File directory) {
        this.readConfiguration = readConfiguration;
        this.lane = lane;
        this.directory = directory;
    }

    public AbstractIlluminaTextParser(final ReadConfiguration readConfiguration, final int lane, final File directory,
                                      final boolean treatGroupedDelimitersAsOne) {
        this.readConfiguration = readConfiguration;
        this.lane = lane;
        this.directory = directory;
        this.treatGroupedDelimitersAsOne = treatGroupedDelimitersAsOne;
    }

    /**
     * Jump so that the next record returned will be the first one from the specified tile.
     */
    @Override
    public void seekToTile(final int oneBasedTileNumber) {
        int i;
        for (i = 0; i < files.length; ++i) {
            if (files[i].tile == oneBasedTileNumber) {
                break;
            }
        }
        // If fall off the end of the for loop, parser should contain no elements.
        initializeParser(i);
    }

    /**
     * Prepare to iterate.
     * @param fileIndex 0-based index into this.files of tile file to start with.
     */
    protected void initializeParser(final int fileIndex) {
        final File[] files = new File[this.files.length - fileIndex];
        for (int i = fileIndex; i < this.files.length; ++i) {
            files[i-fileIndex] = this.files[i].file;
        }
        parser = new BasicInputParser(treatGroupedDelimitersAsOne, files);
    }

    /**
     * Read the next record from the list of input files, and load into data argument.
     * @param data Has already been set up appropriately for paired or single end read, barcoded or not.
     */
    @Override
    public void next(final ClusterData data) {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }

        try {
            final String[] fields = parser.next();
            processLine(data, fields);
        }
        catch (RuntimeException re) {
            log.error(re, "Error parsing " + parser.getFileName() + ":" + parser.getCurrentLineNumber() + " - " + parser.getCurrentLine());
            throw re;
        }
    }

    /**
     * Override this method to parse an input line in the appropriate way for the specific file type.
     * @param data Parsed input is stored in this object, which has already been set up appropriately
     * for paired or single end read, and for barcode or not.
     * @param fields Input line, split on whitespace.
     */
    protected abstract void processLine(ClusterData data, String[] fields);

    public boolean hasNext() {
        return parser.hasNext();
    }

    protected File getDirectory() {
        return directory;
    }

    protected static FormatUtil getFormatter() {
        return formatter;
    }

    protected int getLane() {
        return lane;
    }

    /**
     * Concrete class must call this method before calling initializeParser() in order to
     * inform this class what the input file are.
     * @param files Input files plus tile number, sorted by tile number.
     */
    protected void setFiles(final TiledIlluminaFile[] files) {
        this.files = files;
    }

    public ReadConfiguration getReadConfiguration() {
        return readConfiguration;
    }

    public String getCurrentFilename() {
        return parser.getFileName();
    }

    protected void validateLane(final int lane) {
        if (lane != getLane()) {
            throw new PicardException("Lane number mismatch: " + lane + " != " + getLane());
        }
    }
}
