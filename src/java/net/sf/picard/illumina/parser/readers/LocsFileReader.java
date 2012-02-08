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
import net.sf.picard.util.UnsignedTypeUtil;

import java.io.File;
import java.nio.ByteBuffer;

/**
 * The locs file format is one 3 Illumina formats(pos, locs, and clocs) that stores position data exclusively.
 * locs files store position data for successive clusters in 4 byte float pairs, described as follows:
 *     bytes 1-4    : (int?) Version number (1)
 *     bytes 5-8    : 4 byte float equaling 1.0
 *     bytes 9-12   : unsigned int numClusters
 *     bytes 13-16: : X coordinate of first cluster (32-bit float)
 *     bytes 17-20: : Y coordinate of first cluster (32-bit float)
 *
 *     The remaining bytes of the file store the X and Y coordinates of the remaining clusters.
 */

public class LocsFileReader extends AbstractIlluminaPositionFileReader {
    /** Size of the opening file header, this is skipped by the iterator below*/
    private static final int HEADER_SIZE = 12;

    /** The first four bytes of a locs file should equal a little endian 1 */
    private static final int BYTES_1_TO_4 = 1;

    /** The expected version of locs files */
    private static final float VERSION = 1.0f;

    /** An iterator over all of the coordinate values in the file, remember next needs to be called
     * twice per coordinate pair */
    private BinaryFileIterator<Float> bbIterator;

    /** Total clusters in the file as read in the file header */
    private final long numClusters;

    /** The index of the next cluster to be returned */
    private int nextCluster;

    public LocsFileReader(final File file) {
        super(file);

        bbIterator = MMapBackedIteratorFactory.getFloatIterator(HEADER_SIZE, file);
        final ByteBuffer headerBuf = bbIterator.getHeaderBytes();

        final int firstValue = headerBuf.getInt();
        if(firstValue != BYTES_1_TO_4) {
            throw new PicardException("First header byte of locs files should be " + BYTES_1_TO_4 + " value found(" + firstValue + ")");
        }

        final float versionNumber = headerBuf.getFloat();
        if(versionNumber != VERSION) {
            throw new PicardException("First header byte of locs files should be " + VERSION + " value found(" + firstValue + ")");
        }

        numClusters = UnsignedTypeUtil.uIntToLong(headerBuf.getInt());
        bbIterator.assertTotalElementsEqual(numClusters * 2);
    }

    @Override
    protected PositionInfo unsafeNextInfo() {
        final float xVal = bbIterator.next();
        final float yVal = bbIterator.next();
        ++nextCluster;
        return new PositionInfo(xVal, yVal, getLane(), getTile());
    }

    @Override
    protected String makeExceptionMsg() {
        return "LocsFileReader(file=" + getFile().getAbsolutePath() + ", numClusters=" + numClusters + ") ";
    }

    @Override
    public boolean hasNext() {
        return nextCluster < numClusters;
    }

    public void close() {
        bbIterator = null;
    }
}
