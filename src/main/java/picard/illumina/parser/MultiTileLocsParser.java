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

import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import picard.illumina.parser.readers.LocsFileReader;

import java.io.File;
import java.util.Collections;
import java.util.List;

/**
 * Read locs file that contains multiple tiles in a single file.  A tile index is needed to parse this
 * file so that {tile number, cluster number} can be converted into absolute record number in file.
 */
public class MultiTileLocsParser extends MultiTileParser<PositionalData> {
    private final LocsFileReader reader;
    private final int lane;

    public MultiTileLocsParser(final TileIndex tileIndex, final List<Integer> requestedTiles, final File locsFile, final int lane) {
        super(tileIndex, requestedTiles, Collections.singleton(IlluminaDataType.Position));
        final int tileNumber;
        if (requestedTiles.size() == 1) tileNumber = requestedTiles.get(0);
        else tileNumber = -1;
        this.reader = new LocsFileReader(locsFile, lane, tileNumber);
        this.lane = lane;
    }

    @Override
    PositionalData readNext() {
        final int tile = getTileOfNextCluster();
        final AbstractIlluminaPositionFileReader.PositionInfo nextVal = reader.next();
        return new PositionalData() {
            public int getXCoordinate() {
                return nextVal.xQseqCoord;
            }

            public int getYCoordinate() {
                return nextVal.yQseqCoord;
            }

            public int getLane() {
                return lane;
            }

            public int getTile() {
                return tile;
            }
        };
    }

    @Override
    void skipRecords(final int numToSkip) {
        reader.skipRecords(numToSkip);
    }

    @Override
    public void close() {
        reader.close();
    }
}
