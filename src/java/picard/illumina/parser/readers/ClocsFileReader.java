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
 * The clocs file format is one of 3 Illumina formats(pos, locs, and clocs) that stores position data exclusively.
 * clocs files store position data for successive clusters, compressed in bins as follows:
 *     Byte 0   : unused
 *     Byte 1-4 : unsigned int numBins
 *     The rest of the file consists of bins/blocks, where a bin consists of an integer
 *     indicating number of blocks, followed by that number of blocks and a block consists
 *     of an x-y coordinate pair.  In otherwords:
 *
 *     for each bin
 *         byte 1: Unsigned int numBlocks
 *         for each block:
 *             byte 1 : byte xRelativeCoordinate
 *             byte 2 : byte yRelativeCoordinate
 *
 *     Actual x and y values are computed using the following algorithm
 *
 *     xOffset = yOffset = 0
 *     imageWidth = 2048
 *     blockSize = 25
 *     maxXbins:Int = Math.Ceiling((double)ImageWidth/(double)blockSize)
 *     for each bin:
 *         for each location:
 *             x = convert.ToSingle(xRelativeCoordinate/10f + xoffset)
 *             y = convert.toSingle(yRelativeCoordinate/10f + yoffset)
 *         if (binIndex > 0 && ((binIndex + 1) % maxXbins == 0)) {
 *            xOffset = 0; yOffset += blockSize
 *         } else xOffset += blockSize
 */
public class ClocsFileReader extends AbstractIlluminaPositionFileReader {

    private static final int HEADER_SIZE = 5;

    private static final int IMAGE_WIDTH = 2048;
    private static final int BLOCK_SIZE = 25;
    private static final int NUM_BINS_IN_ROW = (int)Math.ceil((double) IMAGE_WIDTH / (double) BLOCK_SIZE);

    /** Total number of bins */
    private final long numBins;

    /** An iterator through clocsFile's bytes */
    private final BinaryFileIterator<Byte> byteIterator;

    //mutable vars
    private float xOffset;
    private float yOffset;
    private long currentBin;
    private int numClustersInBin;   //MAX 255
    private long currentClusterInBin;

    public ClocsFileReader(final File clocsFile) {
        super(clocsFile);

        byteIterator = MMapBackedIteratorFactory.getByteIterator(HEADER_SIZE, clocsFile);

        final ByteBuffer hbs = byteIterator.getHeaderBytes();
        hbs.get(); //unusedByte
        numBins = UnsignedTypeUtil.uIntToLong(hbs.getInt());

        xOffset = 0;
        yOffset = 0;
        currentBin = 0;
        startBlock();

        checkAndAdvanceBin();
    }

    /**
     * Grab the next set off offset values, decompress them and then
     * @return
     */
    @Override
    protected PositionInfo unsafeNextInfo() {
        final byte xByte = byteIterator.next();
        final byte yByte = byteIterator.next();

        final float xPos = UnsignedTypeUtil.uByteToInt(xByte)/10f + xOffset;
        final float yPos = UnsignedTypeUtil.uByteToInt(yByte)/10f + yOffset;
        ++currentClusterInBin;
        checkAndAdvanceBin();

        return new PositionInfo(xPos, yPos, getLane(), getTile());
    }

    /** Compute offset for next bin and then increment the bin number and reset block information*/
    private void checkAndAdvanceBin() {
        while(currentClusterInBin >= numClustersInBin && currentBin < numBins) { //While rather than if statement to skip empty blocks
            if ((currentBin + 1) % NUM_BINS_IN_ROW == 0) {
                xOffset = 0;
                yOffset += BLOCK_SIZE;
            } else{
                xOffset += BLOCK_SIZE;
            }

            currentBin += 1;
            if(currentBin < numBins) {
                startBlock();
            }
        }
    }

    /** Start the next block by reading it's numBlocks byte and setting the currentBlock index to 0 */
    private void startBlock() {
        numClustersInBin = UnsignedTypeUtil.uByteToInt(byteIterator.next());
        currentClusterInBin = 0;
    }
    
    @Override
    protected String makeExceptionMsg() {
        return "ClocsFileReader(file=" + getFile().getName() + ", lane=" + getLane() +  ", tile=" + getTile() +
                               ", currentBin="   + currentBin +  ", numBins=" + numBins +  ", xOffset=" + xOffset  + ", yOffset" + yOffset +
                               ", currentBlock=" + currentClusterInBin + ", numBlocks=" + numClustersInBin;
    }

    @Override
    public boolean hasNext() { 
        boolean valuesRemain = currentClusterInBin < numClustersInBin || currentBin < (numBins -1);
        if(!valuesRemain && byteIterator.hasNext()) {
            throw new PicardException("Read the number of expected bins( " + numBins + ") but still had more elements in file( " + byteIterator.getFile().getAbsolutePath() + ") ");
        }
        return valuesRemain;
    }

    public void close() {
    }
}
