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

import net.sf.samtools.util.BinaryCodec;

import java.io.*;
import java.util.List;

/**
 * Class for writing binary BAM index files
 */
class BinaryBAMIndexWriter implements BAMIndexWriter {

    protected final int nRef;
    private final BinaryCodec codec;
    private int count = 0;

    /**
     * constructor
     *
     * @param nRef    Number of reference sequences
     * @param output  BAM Index output file
     */
    public BinaryBAMIndexWriter(final int nRef, final File output) {

        this.nRef = nRef;

        try {
            codec = new BinaryCodec(output, true);
            writeHeader();
        } catch (Exception e) {
            throw new SAMException("Exception opening output file " + output, e);
        }
    }

    /**
     *
     * @param nRef Number of reference sequences.
     * @param output BAM index output stream.  This stream will be closed when BinaryBAMIndexWriter.close() is called.
     */
    public BinaryBAMIndexWriter(final int nRef, final OutputStream output) {

        this.nRef = nRef;

        try {
            codec = new BinaryCodec(output);
            writeHeader();
        } catch (Exception e) {
            throw new SAMException("Exception opening output stream", e);
        }
    }

    /**
     * Write this content as binary output
     */
    public void writeReference(final BAMIndexContent content) {

        if (content == null) {
            writeNullContent();
            count++;
            return;
        }

        if (content.getReferenceSequence() != count){
            throw new SAMException("Unexpectedly writing reference " + content.getReferenceSequence() +
                ", expecting reference " + count);
        }
        count ++;

        // write bins

        final BAMIndexContent.BinList bins = content.getBins();
        final int size = bins == null ? 0 : content.getNumberOfNonNullBins();

        if (size == 0) {
            writeNullContent();
            return;
        }

        //final List<Chunk> chunks = content.getMetaData() == null ? null
        //        : content.getMetaData().getMetaDataChunks();
        BAMIndexMetaData metaData = content.getMetaData();

        codec.writeInt(size + ((metaData != null)? 1 : 0 ));
        // codec.writeInt(size);
        for (Bin bin : bins) {   // note, bins will always be sorted
            if (bin.getBinNumber() == AbstractBAMFileIndex.MAX_BINS)
                continue;
            writeBin(bin);
        }

        // write metadata "bin" and chunks        
        if (metaData != null)
            writeChunkMetaData(metaData);

        // write linear index

        final LinearIndex linearIndex = content.getLinearIndex();
        final long[] entries = linearIndex == null ? null : linearIndex.getIndexEntries();
        final int indexStart = linearIndex == null ? 0 : linearIndex.getIndexStart();
        final int n_intv = entries == null ? indexStart : entries.length + indexStart;
        codec.writeInt(n_intv);
        if (entries == null) {
            return;
        }
        // since indexStart is usually 0, this is usually a no-op
        for (int i = 0; i < indexStart; i++) {
            codec.writeLong(0);
        }
        for (int k = 0; k < entries.length; k++) {
            codec.writeLong(entries[k]);
        }
        try {
            codec.getOutputStream().flush();
        } catch (IOException e) {
            throw new SAMException("IOException in BinaryBAMIndexWriter reference " + content.getReferenceSequence(), e);
        }
    }

    /**
     * Writes out the count of records without coordinates
     *
     * @param count
     */
    public void writeNoCoordinateRecordCount(final Long count) {
        codec.writeLong(count == null ? 0 : count);
    }

    /**
     * Any necessary processing at the end of the file
     */
    public void close() {
        codec.close();
    }

    private void writeBin(Bin bin) {
        final int binNumber = bin.getBinNumber();
        if (binNumber >= AbstractBAMFileIndex.MAX_BINS){
            throw new SAMException("Unexpected bin number when writing bam index " + binNumber);
        }
        
        codec.writeInt(binNumber);
        if (bin.getChunkList() == null){
            codec.writeInt(0);
            return;
        }
        final List<Chunk> chunkList = bin.getChunkList();
        final int n_chunk = chunkList.size();
        codec.writeInt(n_chunk);
        for (final Chunk c : chunkList) {
            codec.writeLong(c.getChunkStart());
            codec.writeLong(c.getChunkEnd());
        }
    }

    /**
     * Write the meta data represented by the chunkLists associated with bin MAX_BINS 37450
     *
     * @param metaData information describing numAligned records, numUnAligned, etc
     */
    private void writeChunkMetaData(BAMIndexMetaData metaData) {
        codec.writeInt(AbstractBAMFileIndex.MAX_BINS);
        final int nChunk = 2;
        codec.writeInt(nChunk);
        codec.writeLong(metaData.getFirstOffset());
        codec.writeLong(metaData.getLastOffset());
        codec.writeLong(metaData.getAlignedRecordCount());
        codec.writeLong(metaData.getUnalignedRecordCount());

    }

    private void writeHeader() {
        // magic string
        final byte[] magic = BAMFileConstants.BAM_INDEX_MAGIC;
        codec.writeBytes(magic);
        codec.writeInt(nRef);
    }

    private void writeNullContent() {
        codec.writeLong(0);  // 0 bins , 0 intv
    }
}
