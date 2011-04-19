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
import net.sf.picard.util.IlluminaUtil;

import java.io.File;
import java.util.List;

/**
 * Parser for s_<lane>_<tile>_prb.txt(.gz)? files, produced by Illumina 1.1
 * These files contain, among other things, base quality scores.  All that is extracted here
 * is quality scores.  There is a single line for the entire read, so this class splits up the data
 * into first end, second end (if appropriate) and barcode end (if appropriate).
 * 
 * @author alecw@broadinstitute.org
 */
public class PrbParser extends AbstractIlluminaTextParser {

    /**
     * Preparse to parse prb files.
     * @param readConfiguration Tells how to convert btw cycle number and appropriate end and offset.
     * @param directory Where to find prb files.
     * @param lane
     * @param tiles List of tiles to process in order, or null to indicate that all tiles should be processed.
     */
    public PrbParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                     final List<Integer> tiles) {
        super(readConfiguration, lane, directory);
        setFiles(IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(directory, "prb", lane, tiles));
        initializeParser(0);
    }

    /**
     * Extract the quality values from a line of input, convert to phred binary, and store in the appropriate
     * IlluminaEndData object.  Each line of input contains 4 values for each cycle.  The highest value of the 4 values
     * for the cycle is considered the quality for the cycle.
     *
     * @param data Parsed input is stored in this object, which has already been set up appropriately
     * for paired or single end read, and for barcode or not.
     * @param fields Input line, split on whitespace.
     */
    @Override
    protected void processLine(final IlluminaReadData data, final String[] fields) {
        final int expectedLength = getReadConfiguration().getFirstLength() + getReadConfiguration().getSecondLength() +
                getReadConfiguration().getBarcodeLength();
        if (expectedLength  * 4 != fields.length) {
            throw new PicardException("Number of fields does not match expected in " + getCurrentFilename());
        }

        byte[] solexaQualityChars = new byte[getReadConfiguration().getFirstLength()];
        int startCycle = getReadConfiguration().getFirstStart();
        for (int i = 0; i < getReadConfiguration().getFirstLength(); ++i) {
            solexaQualityChars[i] = IlluminaUtil.getSolexaQualityCharFromFourQualities(fields, startCycle + i, getFormatter());
        }
        IlluminaUtil.convertSolexaQualityAscii_1_1_ToPhredBinary(solexaQualityChars);
        data.getFirstEnd().setQualities(solexaQualityChars);

        if (getReadConfiguration().isPairedEnd()) {
            solexaQualityChars = new byte[getReadConfiguration().getSecondLength()];
            startCycle = getReadConfiguration().getSecondStart();
            for (int i = 0; i < getReadConfiguration().getSecondLength(); ++i) {
                solexaQualityChars[i] = IlluminaUtil.getSolexaQualityCharFromFourQualities(fields, startCycle + i, getFormatter());
            }
            IlluminaUtil.convertSolexaQualityAscii_1_1_ToPhredBinary(solexaQualityChars);
            data.getSecondEnd().setQualities(solexaQualityChars);
        }
        if (getReadConfiguration().isBarcoded()) {
            solexaQualityChars = new byte[getReadConfiguration().getBarcodeLength()];
            startCycle = getReadConfiguration().getBarcodeStart();
            for (int i = 0; i < getReadConfiguration().getBarcodeLength(); ++i) {
                solexaQualityChars[i] = IlluminaUtil.getSolexaQualityCharFromFourQualities(fields, startCycle + i, getFormatter());
            }
            IlluminaUtil.convertSolexaQualityAscii_1_1_ToPhredBinary(solexaQualityChars);
            data.getBarcodeRead().setQualities(solexaQualityChars);
        }
    }
}
