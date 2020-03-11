/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package picard.illumina.parser;

import picard.illumina.parser.readers.FilterFileReader;

import java.io.File;
import java.util.Collections;
import java.util.List;

/**
 * Read filter file that contains multiple tiles in a single file.  A tile index is needed to parse this
 * file so that {tile number, cluster number} can be converted into absolute record number in file.
 */
public class MultiTileFilterParser extends MultiTileParser<PfData> {
    private final FilterFileReader reader;

    public MultiTileFilterParser(final TileIndex tileIndex, final List<Integer> requestedTiles, final File filterFile) {
        super(tileIndex, requestedTiles, Collections.singleton(IlluminaDataType.PF));
        reader = new FilterFileReader(filterFile);
    }

    @Override
    PfData readNext() {
        final boolean nextVal = reader.next();
        return new PfData() {
            @Override
            public boolean isPf() {
                return nextVal;
            }
        };
    }

    @Override
    void skipRecords(final int numToSkip) {
        reader.skipRecords(numToSkip);
    }

    @Override
    public void close() {
        //no-op
    }
}
