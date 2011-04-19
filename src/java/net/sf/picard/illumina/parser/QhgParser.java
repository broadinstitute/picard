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

/**
 * Parser for s_<lane>_<tile>_qhg.txt(.gz)? files, produced by Illumina 1.1
 * Currently the only data parsed from this file is Passing Filter score.  If the score is > PASSING_FILTER_THRESHOLD
 * then the read is considered to be PF.
 * 
 * @author alecw@broadinstitute.org
 */
public class QhgParser extends AbstractIlluminaTextParser {
    private static final double PASSING_FILTER_THRESHOLD = 0.6;
    private static final int LANE_COLUMN = 0;
    private static final int TILE_COLUMN = 1;
    private static final int X_COLUMN = 2;
    private static final int Y_COLUMN = 3;
    private static final int PF_SCORE_COLUMN = 5;

    /**
     * Prepare to parse a directory of qhg files.
     * @param readConfiguration Describes correcpondence btw cycle number and end, not needed by this implementation
     * because PF is global for the entire read.
     * @param directory Where to find qhg files.
     * @param lane
     * @param tiles List of tiles to process, in order, or null to process all tiles in the directory.
     */
    public QhgParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                     final List<Integer> tiles) {
        super(readConfiguration, lane, directory);
        setFiles(IlluminaFileUtil.getNonEndedIlluminaBasecallFiles(directory, "qhg", lane, tiles));
        initializeParser(0);
    }

    /**
     * Extract PF score for a read, and set pf attribute of IlluminaReadData according to whether
     * the score is >= PASSING_FILTER_THRESHOLD or not.
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
        final double pfScore = getFormatter().parseDouble(fields[PF_SCORE_COLUMN]);
        final boolean passingFilter = (pfScore >= PASSING_FILTER_THRESHOLD);

        data.setOrCheckLane(lane);
        data.setOrCheckTile(tile);
        data.setOrCheckX(x);
        data.setOrCheckY(y);
        data.setOrCheckPf(passingFilter);
    }
}
