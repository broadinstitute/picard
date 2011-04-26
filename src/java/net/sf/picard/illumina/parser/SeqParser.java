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
import net.sf.picard.PicardException;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.List;

/**
 * Parser for Illumina 1.1 s_<lane>_<tile>_seq.txt(.gz)? files.
 * Each record contains:
 *   lane
 *   tile
 *   x-coordinate
 *   y-coordinate
 *   bases for all cycles
 * 
 * @author alecw@broadinstitute.org
 */
public class SeqParser extends AbstractIlluminaTextParser {

    private static final int LANE_COLUMN = 0;
    private static final int TILE_COLUMN = 1;
    private static final int X_COLUMN = 2;
    private static final int Y_COLUMN = 3;
    private static final int BASES_COLUMN = 4;

    /**
     * Prepare to parse seq files.
     * @param readConfiguration Defines correspondence btw cycle number and end.
     * @param directory Where to find seq files.
     * @param lane
     * @param tiles The list of tiles to be processed, in order, or null to process all tiles.
     */
    public SeqParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                     final List<Integer> tiles) {
        super(readConfiguration, lane, directory);
        setFiles(IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(directory, "seq", lane, tiles));
        initializeParser(0);
    }

    /**
     * Process a line of seq data.  The bases are split as appropriate into first end, second end and barcode end,
     * and set into the bases property of the corresponding IlluminaEndData in the IlluminaReadData argument.
     * 
     * @param data Parsed input is stored in this object, which has already been set up appropriately
     * for paired or single end read, and for barcode or not.
     * @param fields Input line, split on whitespace.
     */
    @Override
    protected void processLine(final IlluminaReadData data, final String[] fields) {
        final int lane = getFormatter().parseInt(fields[LANE_COLUMN]);
        validateLane(lane);
        final int tile = getFormatter().parseInt(fields[TILE_COLUMN]);
        final int x = getFormatter().parseInt(fields[X_COLUMN]);
        final int y = getFormatter().parseInt(fields[Y_COLUMN]);

        final String baseString = fields[BASES_COLUMN];
        final int expectedLength = getReadConfiguration().getFirstLength() + getReadConfiguration().getSecondLength() +
                getReadConfiguration().getBarcodeLength();
        if (expectedLength != baseString.length()) {
            throw new PicardException("Length of bases does not match expected in " + getCurrentFilename());
        }

        data.setOrCheckLane(lane);
        data.setOrCheckTile(tile);
        data.setOrCheckX(x);
        data.setOrCheckY(y);

        final byte[] bases1 = StringUtil.stringToBytes(baseString, getReadConfiguration().getFirstStart() - 1,
                getReadConfiguration().getFirstLength());
        data.getFirstEnd().setBases(bases1);

        if (getReadConfiguration().isPairedEnd()) {
            final byte[] bases2 = StringUtil.stringToBytes(baseString, getReadConfiguration().getSecondStart() - 1,
                    getReadConfiguration().getSecondLength());
            data.getSecondEnd().setBases(bases2);
        }
        if (getReadConfiguration().isBarcoded()) {
            final byte[] barcodeBases = StringUtil.stringToBytes(baseString, getReadConfiguration().getBarcodeStart() - 1,
                    getReadConfiguration().getBarcodeLength());
            data.getBarcodeRead().setBases(barcodeBases);
        }
    }

    public static int getReadLength(final File seqFile) {
        final BasicInputParser parser = new BasicInputParser(true, seqFile);
        if (!parser.hasNext()) {
            throw new PicardException("Unexpected empty qseq file: " + seqFile);
        }
        final String[] fields = parser.next();
        return fields[BASES_COLUMN].length();
    }
}
