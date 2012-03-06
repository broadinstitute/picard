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

import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * Interface for classes that parse information out of the Illumina Pipeline
 *
 * @author jburke@broadinstitute.org
 */
interface IlluminaParser<DATA_TYPE extends IlluminaData> extends Iterator<DATA_TYPE> {
    /** Jump so that the next record returned will be from the specified tile. */
    void seekToTile(int oneBasedTileNumber);

    /**
     * Read the next read's set of data and set it into the provided data object.  The object must have
     * the appropriate IlluminaEndData objects set into it for first end, second end, barcode.
     */
    DATA_TYPE next();

    /**
     * Is there a DATA_TYPE object for another cluster remaining.
     * @return TRUE if there is a DATA_TYPE object for the next cluster that can be provided by
     * next
     */
    boolean hasNext();

    /** Get the tile for the NEXT DATA_TYPE object that will be returned by this parser.  This should
     * be called BEFORE next if you want the tile for the value returned by next */
    public int getTileOfNextCluster();

    /**
     * Given the expected tiles and cycles for this run, make sure this parser can provide data for
     * all tiles/cycles or throw's an exception if it's missing any required data or data structures
     * it relies on do not disagree with the provided tiles/cycles
     * @param tiles The number of tiles in the current run
     * @param cycles The number of cycles in the current run
     */
    void verifyData(final List<Integer> tiles, final int [] cycles);

    /** The DataTypes that this parser can provide */
    Set<IlluminaDataType> supportedTypes();

}
