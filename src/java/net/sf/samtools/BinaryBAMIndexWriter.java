/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

package net.sf.samtools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.Arrays;
import java.util.List;

/**
 * Class for writing binary BAM index files
 */
public class BinaryBAMIndexWriter extends AbstractBAMIndexWriter {

    private final int bufferSize = 1000000;
    private final ByteBuffer bb;
    private final FileChannel fileChannel;
    private final FileOutputStream stream;
    private final boolean sortBins;

    /**
     * constructor
     *
     * @param n_ref    Number of reference sequences
     * @param output   BAM Index output file
     * @param sortBins Whether to sort the bins - useful for comparison to c-generated index
     */
    public BinaryBAMIndexWriter(final int n_ref, final File output, boolean sortBins) {
        super(output, n_ref);

        this.sortBins = sortBins;

        try {
            stream = new FileOutputStream(output, true);
            fileChannel = stream.getChannel();
            bb = ByteBuffer.allocateDirect(bufferSize);
            bb.order(ByteOrder.LITTLE_ENDIAN);
        } catch (FileNotFoundException e) {
            throw new SAMException("Can't find output file " + output, e);
        }
    }


    public void writeHeader() {
        // magic string
        final byte[] magic = BAMFileConstants.BAM_INDEX_MAGIC;
        bb.put(magic);
        // n_ref
        bb.putInt(n_ref);
    }

    /**
     * Write this content as binary output
     */
    public void writeReference(final BAMIndexContent content, int reference) {

        if (content == null) {
            writeNullContent(bb);
            return;
        }
        final List<Bin> bins = content.getBins();
        final LinearIndex linearIndex = content.getLinearIndex();

        if (bins == null || bins.size() == 0) {
            writeNullContent(bb);
            return;
        }

        final int size = bins.size();
        bb.putInt(size);

        // copy bins into an array so that it can be sorted
        final Bin[] binArray = new Bin[size];
        if (size != 0) {
            bins.toArray(binArray);
        }
        if (sortBins) Arrays.sort(binArray);  // Sort for easy text comparisons

        // todo - don't copy the array when no sorting
        // and instead loop using
        // for (Bin bin: bins) {

        for (int j = 0; j < size; j++) {
            Bin bin = binArray[j];

            if (bin.getBinNumber() == BAMIndex.MAX_BINS)  break;
            
            bb.putInt(bin.getBinNumber()); // todo uint32_t vs int32_t in spec?
            if (bin.getChunkList() == null){
                bb.putInt(0);
                continue;
            }
            final List<Chunk> chunkList = bin.getChunkList();
            final int n_chunk = chunkList.size();
            bb.putInt(n_chunk);
            for (final Chunk c : chunkList) {
                bb.putLong(c.getChunkStart());   // todo uint32_t vs int32_t in spec?
                bb.putLong(c.getChunkEnd());     // todo uint32_t vs int32_t in spec?
            }
        }
        writeChunkMetaData(content.getMetaDataChunks());

        final long[] entries = linearIndex == null ? null : linearIndex.getIndexEntries();
        final int indexStart = linearIndex == null ? 0 : linearIndex.getIndexStart();
        final int n_intv = entries == null ? indexStart : entries.length + indexStart;
        bb.putInt(n_intv);
        if (entries == null) {
            return;
        }

        for (int i = 0; i < indexStart; i++) {
            bb.putLong(0);          // todo uint32_t vs int32_t in spec?
        }
        for (int k = 0; k < entries.length; k++) {
            bb.putLong(entries[k]); // todo uint32_t vs int32_t in spec?
        }
        // write out data and reset the buffer for each reference
        bb.flip();
        try {
            fileChannel.write(bb);
            stream.flush();
        } catch (IOException e) {
            throw new SAMException("IOException in BinaryBAMIndexWriter reference " + reference, e);
        }
        bb.position(0);
        bb.limit(bufferSize);
    }

    /**
     * Write the meta data represented by the chunkLists associated with bin MAX_BINS 37450
     *
     * @param chunkList contains metadata describing numAligned records, numUnAligned, etc
     */
    private void writeChunkMetaData(List<Chunk> chunkList) {
        bb.putInt(BAMIndex.MAX_BINS);
        final int n_chunk = chunkList.size();   // should be 2
        if (n_chunk != 2){
            System.err.println("Undexpect # chunks of meta data= " + n_chunk); // throw new SAMException
        }
        bb.putInt(n_chunk);
        for (final Chunk c : chunkList) {
            bb.putLong(c.getChunkStart());   // todo uint32_t vs int32_t in spec?
            bb.putLong(c.getChunkEnd());     // todo uint32_t vs int32_t in spec?
        }
   }


    private static void writeNullContent(ByteBuffer bb) {
        bb.putInt(0);  // 0 bins
        bb.putInt(0);  // 0 intv
    }

    public void close(Long noCoordinateCount) {
        bb.putLong(noCoordinateCount == null ? 0 : noCoordinateCount);
        bb.flip();
        try {
            fileChannel.write(bb);
            fileChannel.close();
            stream.flush();
            stream.close();
        } catch (IOException e) {
            throw new SAMException("IOException in BinaryBAMIndexWriter ", e);
        }
    }
}