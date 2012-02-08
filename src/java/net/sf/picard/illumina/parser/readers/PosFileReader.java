/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
package net.sf.picard.illumina.parser.readers;

import net.sf.picard.PicardException;
import net.sf.picard.util.BasicInputParser;
import net.sf.samtools.util.CloserUtil;

import java.io.File;

/**
 * The pos file format is one 3 Illumina formats(pos, locs, and clocs) that stores position data exclusively.
 * pos files store position data for successive clusters in tabbed delimited coordinated pairs, 1 per file row e.g.:
 *
 * xPos1\tyPos1
 * xPos2\tyPos2
 * 102.0\t303.3
 *     ...
 * xPosn-1\yPosn-1
 * xPosn\yPosn
 *
 * Where n = the total number of clusters (and therefore lines) in the file.
 */
public class PosFileReader extends AbstractIlluminaPositionFileReader {

    private final BasicInputParser parser;

    public PosFileReader(final File posFile) {
        super(posFile);
        this.parser  = new BasicInputParser(true, posFile);
    }

    /** Read a line of text and parse it into two float values, create a PositionInfo and return it */
    @Override
    protected PositionInfo unsafeNextInfo() {
        final String [] strVals = this.parser.next();
        if(strVals.length != 2) {
            throw new PicardException("Pos file number of values != 2, found (" + strVals.length +")" + makeExceptionMsg());
        }
        try {
            final float xVal = Float.parseFloat(strVals[0]);
            final float yVal = Float.parseFloat(strVals[1]);

            if(xVal <0 || yVal < 0) {
                throw new NumberFormatException("X and Y pos values cannot be negative!");
            }

            return new PositionInfo(xVal, yVal, getLane(), getTile());
        } catch(final NumberFormatException nfe) {
            throw new PicardException("Bad x or y value in " + makeExceptionMsg(), nfe);
        }
    }

    @Override
    protected String makeExceptionMsg() {
        return "pos file( "              + parser.getFileName()            +
               " ) on line number( "     + parser.getCurrentLineNumber()   +
               " ) with current line = " + parser.getCurrentLine();
    }

    public boolean hasNext() {
        return parser.hasNext();
    }

    public void close() {
        CloserUtil.close(parser);
    }
}
