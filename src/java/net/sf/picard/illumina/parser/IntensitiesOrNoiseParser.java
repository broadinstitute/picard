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

import java.io.File;
import java.util.List;
import java.util.Arrays;

/**
 * Abstract base class for parsing files of the form s_<lane>_<tile>_<type>.txt.  Typically these files store
 * raw intensities, processed intensities, or noise, which are file types "int", "sig2", and "nse" respectively.
 * Note that although the input files contain floating point values, these are truncated and stored as shorts.
 * These files contain data for all cycles, not just a single end.
 *
 * These files contain the following
 *     lane
 *     tile
 *     x-coordinate
 *     y-coordinate
 *     A value for cycle 1
 *     C value for cycle 1
 *     G value for cycle 1
 *     T value for cycle 1
 *     ... etc for each cycle
 *
 *
 *
 * @author alecw@broadinstitute.org
 */
public abstract class IntensitiesOrNoiseParser extends AbstractIlluminaTextParser {

    private static final int COLUMNS_PER_BASE = 4;
    // Index of first intensity column
    private static final int INTENSITIES_OFFSET = 4;

    /**
     * Prepare to parse an Illumina text file containing 4 channels of values for each cycle.
     * @param readConfiguration Defines correspondence btw cycle # and end.
     * @param directory Where to find the input files.
     * @param lane
     * @param fileType Typically "int", "sig2" or "nse".
     * @param tiles List of tiles to parse, in order, or null to indicate that all tiles should be parsed.
     */
    public IntensitiesOrNoiseParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                             final String fileType, final List<Integer> tiles) {
        super(readConfiguration, lane, directory);
        setFiles(IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(directory, fileType, lane, tiles));
        initializeParser(0);
    }

    /**
     * Process a single line of input.  Splits the input into first end, second end (if appropriate) and barcode
     * end (if appropriate), and then calls a method on the concrete class to stuff the parsed data into the correct
     * slot in IlluminaEndData.
     * @param data Parsed input is stored in this object, which has already been set up appropriately
     * for paired or single end read, and for barcode or not.
     * @param fields Input line, split on whitespace.
     */
    @Override
    protected void processLine(final IlluminaReadData data, final String[] fields) {
        final int lane = getFormatter().parseInt(fields[0]);
        validateLane(lane);
        final int tile = getFormatter().parseInt(fields[1]);
        final int x = getFormatter().parseInt(fields[2]);
        final int y = getFormatter().parseInt(fields[3]);
        data.setOrCheckLane(lane);
        data.setOrCheckTile(tile);
        data.setOrCheckX(x);
        data.setOrCheckY(y);

        // 0-based
        final int firstStart = convertCycleNumberToIntensityIndex(getReadConfiguration().getFirstStart());
        final int firstEnd = convertCycleNumberToIntensityIndex(getReadConfiguration().getFirstEnd() + 1);
        setValues(data.getFirstEnd(), IlluminaFileUtil.parseSig2Intensities(fields, firstStart, firstEnd));

        if (getReadConfiguration().isPairedEnd()) {
            // 0-based
            final int secondStart = convertCycleNumberToIntensityIndex(getReadConfiguration().getSecondStart());
            final int secondEnd = convertCycleNumberToIntensityIndex(getReadConfiguration().getSecondEnd() + 1);
            setValues(data.getSecondEnd(), IlluminaFileUtil.parseSig2Intensities(fields, secondStart, secondEnd));
        }

        if (getReadConfiguration().isBarcoded()) {
            // 0-based
            final int barcodeStart = convertCycleNumberToIntensityIndex(getReadConfiguration().getBarcodeStart());
            final int barcodeEnd = convertCycleNumberToIntensityIndex(getReadConfiguration().getBarcodeEnd() + 1);
            setValues(data.getBarcodeRead(), IlluminaFileUtil.parseSig2Intensities(fields, barcodeStart, barcodeEnd));
        }
    }

    /**
     * Convert btw 1-based cycle number to 0-based offset into array of fields in an input line.
     * @param cycleNumber 1-based cycle number.
     * @return 0-based index into array of input fields.
     */
    private int convertCycleNumberToIntensityIndex(final int cycleNumber) {
        return (cycleNumber - 1) * COLUMNS_PER_BASE + INTENSITIES_OFFSET;
    }

    /**
     * Override this method to store the data in the appropriate slot.
     * @param end Object on which to store the data.
     * @param values Data to store.
     */
    abstract protected void setValues(IlluminaEndData end, FourChannelIntensityData values);

    protected static boolean fileExists(final File directory, final int lane, final String fileType, final int tile) {
        try {
        return IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(directory, fileType, lane, Arrays.asList(tile)).length > 0;
        } catch (IlluminaFileNotFoundException ex) {
            return false;
        }
    }
}
