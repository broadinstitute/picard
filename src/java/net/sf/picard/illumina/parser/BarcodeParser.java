/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

/**
 * Parse barcode file and store the matched barcode (not the actual sequence from the read),
 * in IlluminaReadData.matchedBarcode.  If read does not match a barcode, this attribute is not set,
 * and should be null.
 * 
 * @author alecw@broadinstitute.org
 */
public class BarcodeParser extends AbstractIlluminaTextParser {
    private static final int Y_OR_N_COLUMN = 1;
    private static final int BARCODE_COLUMN = 2;

    public BarcodeParser(final ReadConfiguration readConfiguration, final File basecallDirectory, final int lane,
                         final List<Integer> tiles) {
        super(readConfiguration, lane, basecallDirectory, false);
        setFiles(IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(basecallDirectory, "barcode", lane, tiles));
        initializeParser(0);
    }

    /**
     * Parse a line from the s_<lane>_barcode.txt file.
     *
     * @param data   Object on which matchedBarcode property will be set if this line has a barcode match.
     * @param fields Input line, split on whitespace.
     */
    @Override
    protected void processLine(final IlluminaReadData data, final String[] fields) {
        if (fields[Y_OR_N_COLUMN].equals("Y")) {
            data.setMatchedBarcode(fields[BARCODE_COLUMN]);
        }
    }
}
