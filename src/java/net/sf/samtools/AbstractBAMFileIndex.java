/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import net.sf.samtools.util.RuntimeIOException;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.MappedByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

/**
 * Provides basic, generic capabilities to be used reading BAM index files.  Users can
 * subclass this class to create new BAM index functionality for adding querying facilities,
 * changing caching behavior, etc.
 *
 * Of particular note: the AbstractBAMFileIndex is, by design, the only class aware of the
 * details of the BAM index file format.  Anyone wanting to implement a reader for a differing
 * or extended BAM index format should implement BAMIndex directly.
 */
abstract class AbstractBAMFileIndex implements BAMIndex {
    /**
     * Reports the maximum number of bins that can appear in a BAM file.
     */
    private static final int MAX_BINS = 37450; // =(8^6-1)/7+1

    /**
     * Reports the total amount of genomic data that any bin can index.
     */
    protected static final int BIN_SPAN = 512*1024*1024;

    /**
     * What is the starting bin for each level?
     */
    private static final int[] LEVEL_STARTS = {0,1,9,73,585,4681};

    private final File mFile;
    private MappedByteBuffer mFileBuffer;

    protected AbstractBAMFileIndex(final File file) {
        mFile = file;
        if (file != null){
            open();
        }
    }

    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    public int getNumIndexLevels() {
        return LEVEL_STARTS.length;
    }

    /**
     * Gets the first bin in the given level.
     * @param levelNumber Level number.  0-based.
     * @return The first bin in this level.
     */
    public int getFirstBinInLevel(final int levelNumber) {
        return LEVEL_STARTS[levelNumber];
    }

    /**
     * Gets the number of bins in the given level.
     * @param levelNumber Level number.  0-based.
     * @return The size (number of possible bins) of the given level.
     */
    public int getLevelSize(final int levelNumber) {
        if(levelNumber == getNumIndexLevels())
            return MAX_BINS+1-LEVEL_STARTS[levelNumber];
        else
            return LEVEL_STARTS[levelNumber+1]-LEVEL_STARTS[levelNumber];
    }

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin  for which to determine the level.
     * @return the level associated with the given bin number.
     */
    public int getLevelForBin(final Bin bin) {
        if(bin.getBinNumber() >= MAX_BINS)
            throw new SAMException("Tried to get level for invalid bin.");
        for(int i = getNumIndexLevels()-1; i >= 0; i--) {
            if(bin.getBinNumber() >= LEVEL_STARTS[i])
                return i;
        }
        throw new SAMException("Unable to find correct bin for bin "+bin);
    }

    /**
     * Gets the first locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getFirstLocusInBin(final Bin bin) {
        final int level = getLevelForBin(bin);
        final int levelStart = LEVEL_STARTS[level];
        final int levelSize = ((level==getNumIndexLevels()-1) ? MAX_BINS-1 : LEVEL_STARTS[level+1]) - levelStart;
        return (bin.getBinNumber() - levelStart)*(BIN_SPAN/levelSize)+1;
    }

    /**
     * Gets the last locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getLastLocusInBin(final Bin bin) {
        final int level = getLevelForBin(bin);
        final int levelStart = LEVEL_STARTS[level];
        final int levelSize = ((level==getNumIndexLevels()-1) ? MAX_BINS-1 : LEVEL_STARTS[level+1]) - levelStart;
        return (bin.getBinNumber()-levelStart+1)*(BIN_SPAN/levelSize);
    }

    public void open() {
        // Open the file stream.
        if (mFileBuffer != null) {
            return;
        }
        try {
            FileInputStream fileStream = new FileInputStream(mFile);
            FileChannel fileChannel = fileStream.getChannel();
            mFileBuffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, 0L, fileChannel.size());
            mFileBuffer.order(ByteOrder.LITTLE_ENDIAN);

            fileChannel.close();
            fileStream.close();
        }
        catch (IOException exc) {
            throw new RuntimeIOException(exc.getMessage(), exc);
        }

        // Verify the magic number.
        seek(0);
        final byte[] buffer = new byte[4];
        readBytes(buffer);
        if (!Arrays.equals(buffer, BAMFileConstants.BAM_INDEX_MAGIC)) {
            close();
            throw new RuntimeException("Invalid file header in BAM index " + mFile +
                                       ": " + new String(buffer));
        }
    }

    public int getNumberOfReferences() {
        if(mFileBuffer == null)
            throw new SAMException("Cannot query a closed index file");
        seek(4);
        return readInteger();
    }

    public void close() {
        mFileBuffer = null;
    }

    /**
     * Write textual bam index file
     */
    public void writeText(final int n_ref, final File OUTPUT, final boolean sortBins) throws Exception {

        final PrintWriter pw = new PrintWriter(OUTPUT);
        pw.println("n_ref=" + n_ref);
        for (int i = 0; i < n_ref; i++) {
            try {
                if (getQueryResults(i) == null) {
                    BAMIndexContent.writeNullTextContent(pw, i);
                    continue;
                }
                getQueryResults(i).writeText(pw, sortBins);
            } catch (Exception e) {
                System.err.println(e.getMessage() + " Exception writing text for reference " + i);
            }
        }
        pw.close();
    }

    /**
     * Write binary bam index file
     *
     * @param
     */
    public void writeBinary(final int n_ref, final File OUTPUT, final boolean sortBins, final long bamFileSize) throws Exception {

        final int bufferSize; //  = 1000000; // 1M  works, but doesn't need to be this big
        final int defaultBufferSize = 1000000;  // 1M
        if (bamFileSize < defaultBufferSize  && bamFileSize != 0) {
            bufferSize = (int) bamFileSize;
        } else {
            bufferSize = defaultBufferSize;
        }
        // log.info("ByteBuffer size is " + bufferSize);

        final FileOutputStream stream = new FileOutputStream(OUTPUT, true);
        final FileChannel fileChannel = stream.getChannel();
        final ByteBuffer bb = ByteBuffer.allocateDirect(bufferSize);
        bb.order(ByteOrder.LITTLE_ENDIAN);

        // magic string
        final byte[] magic = BAMFileConstants.BAM_INDEX_MAGIC;
        bb.put(magic);
        // n_ref
        bb.putInt(n_ref);
        for (int i = 0; i < n_ref; i++) {
            if (getQueryResults(i) == null){
                BAMIndexContent.writeNullBinaryContent(bb);
                continue;
            }
            getQueryResults(i).writeBinary(bb, sortBins);
            //  write out data and reset the buffer for each reference
            bb.flip();
            fileChannel.write(bb, 0);
            // stream.flush();    // todo will flushing the stream at every reference help memory?
            bb.position(0);
            bb.limit(bufferSize);
        }
        bb.flip();
        fileChannel.write(bb, 0);
        fileChannel.close();
        stream.close();
    }

    /**
     * Use to get close to the unmapped reads at the end of a BAM file.
     * @return The file offset of the first record in the last linear bin, or -1
     * if there are no elements in linear bins (i.e. no mapped reads).
     */
    public long getStartOfLastLinearBin() {
        if(mFileBuffer == null)
            throw new SAMException("Cannot query a closed index file");        
        seek(4);

        final int sequenceCount = readInteger();
        // Because no reads may align to the last sequence in the sequence dictionary,
        // grab the last element of the linear index for each sequence, and return
        // the last one from the last sequence that has one.
        long lastLinearIndexPointer = -1;
        for (int i = 0; i < sequenceCount; i++) {
            // System.out.println("# Sequence TID: " + i);
            final int nBins = readInteger();
            // System.out.println("# nBins: " + nBins);
            for (int j1 = 0; j1 < nBins; j1++) {
                // Skip bin #
                skipBytes(4);
                final int nChunks = readInteger();
                // Skip chunks
                skipBytes(16 * nChunks);
            }
            final int nLinearBins = readInteger();
            if (nLinearBins > 0) {
                // Skip to last element of list of linear bins
                skipBytes(8 * (nLinearBins - 1));
                lastLinearIndexPointer = readLong();
            }
        }

        return lastLinearIndexPointer;
    }

    protected BAMIndexContent query(final int referenceSequence, final int startPos, final int endPos) {
        if(mFileBuffer == null)
            throw new SAMException("Cannot query a closed index file");
        seek(4);

        List<Bin> bins = null;
        SortedMap<Bin, List<Chunk>> binToChunks = new TreeMap<Bin,List<Chunk>>();
        LinearIndex linearIndex = null;

        final int sequenceCount = readInteger();

        if (referenceSequence >= sequenceCount) {
            return null;
        }

        final BitSet regionBins = regionToBins(startPos, endPos);
        if (regionBins == null) {
            return null;
        }

        skipToSequence(referenceSequence);

        final int binCount = readInteger();
        bins = new ArrayList<Bin>(binCount);
        for(int binNumber = 0; binNumber < binCount; binNumber++) {
            List<Chunk> chunks = new ArrayList<Chunk>();
            final int indexBin = readInteger();
            final int nChunks = readInteger();
            // System.out.println("# bin[" + i + "] = " + indexBin + ", nChunks = " + nChunks);
            if (regionBins.get(indexBin)) {
                for (int ci = 0; ci < nChunks; ci++) {
                    final long chunkBegin = readLong();
                    final long chunkEnd = readLong();
                    chunks.add(new Chunk(chunkBegin, chunkEnd));
                }
            } else {
                skipBytes(16 * nChunks);
            }
            Bin bin = new Bin(referenceSequence,indexBin);
            bins.add(bin);
            binToChunks.put(bin,chunks);
        }
        // Reorder the bins in binNumber order.
        Collections.sort(bins);

        final int nLinearBins = readInteger();

        final int regionLinearBinStart = LinearIndex.convertToLinearIndexOffset(startPos);
        final int regionLinearBinStop = LinearIndex.convertToLinearIndexOffset(endPos)>0 ? LinearIndex.convertToLinearIndexOffset(endPos) : nLinearBins-1; 

        long[] linearIndexEntries = new long[0];
        if (regionLinearBinStart < nLinearBins) {
            linearIndexEntries = new long[regionLinearBinStop-regionLinearBinStart+1];
            skipBytes(8 * regionLinearBinStart);
            for(int linearBin = regionLinearBinStart; linearBin <= regionLinearBinStop; linearBin++)
                linearIndexEntries[linearBin-regionLinearBinStart] = readLong();
        }

        linearIndex = new LinearIndex(referenceSequence,regionLinearBinStart,linearIndexEntries);

        return new BAMIndexContent(referenceSequence, bins, binToChunks, linearIndex);
    }

    abstract protected BAMIndexContent getQueryResults(int reference);

    /**
     * Gets the possible number of bins for a given reference sequence.
     * @return How many bins could possibly be used according to this indexing scheme to index a single contig.
     */
    protected int getMaxAddressibleGenomicLocation() {
        return BIN_SPAN;
    }

    /**
     * Get candidate bins for the specified region
     * @param startPos 1-based start of target region, inclusive.
     * @param endPos 1-based end of target region, inclusive.
     * @return bit set for each bin that may contain SAMRecords in the target region.
     */
    protected BitSet regionToBins(final int startPos, final int endPos) {
        final int maxPos = 0x1FFFFFFF;
        final int start = (startPos <= 0) ? 0 : (startPos-1) & maxPos;
        final int end = (endPos <= 0) ? maxPos : (endPos-1) & maxPos;
        if (start > end) {
            return null;
        }
        int k;
        final BitSet bitSet = new BitSet(MAX_BINS);
        bitSet.set(0);
        for (k =    1 + (start>>26); k <=    1 + (end>>26); ++k) bitSet.set(k);
        for (k =    9 + (start>>23); k <=    9 + (end>>23); ++k) bitSet.set(k);
        for (k =   73 + (start>>20); k <=   73 + (end>>20); ++k) bitSet.set(k);
        for (k =  585 + (start>>17); k <=  585 + (end>>17); ++k) bitSet.set(k);
        for (k = 4681 + (start>>14); k <= 4681 + (end>>14); ++k) bitSet.set(k);
        return bitSet;
    }    

    protected List<Chunk> optimizeChunkList(final List<Chunk> chunks, final long minimumOffset) {
        Chunk lastChunk = null;
        Collections.sort(chunks);
        final List<Chunk> result = new ArrayList<Chunk>();
        for (final Chunk chunk : chunks) {
            if (chunk.getChunkEnd() <= minimumOffset) {
                continue;
            }
            if (result.isEmpty()) {
                result.add(chunk);
                lastChunk = chunk;
                continue;
            }
            // Coalesce chunks that are in adjacent file blocks.
            // This is a performance optimization.
            final long lastFileBlock = getFileBlock(lastChunk.getChunkEnd());
            final long chunkFileBlock = getFileBlock(chunk.getChunkStart());
            if (chunkFileBlock - lastFileBlock > 1) {
                result.add(chunk);
                lastChunk = chunk;
            } else {
                if (chunk.getChunkEnd() > lastChunk.getChunkEnd()) {
                    lastChunk.setChunkEnd(chunk.getChunkEnd());
                }
            }
        }
        return result;
    }

    private void skipToSequence(final int sequenceIndex) {
        for (int i = 0; i < sequenceIndex; i++) {
            // System.out.println("# Sequence TID: " + i);
            final int nBins = readInteger();
            // System.out.println("# nBins: " + nBins);
            for (int j = 0; j < nBins; j++) {
                final int bin = readInteger();
                final int nChunks = readInteger();
                // System.out.println("# bin[" + j + "] = " + bin + ", nChunks = " + nChunks);
                skipBytes(16 * nChunks);
            }
            final int nLinearBins = readInteger();
            // System.out.println("# nLinearBins: " + nLinearBins);
            skipBytes(8 * nLinearBins);
        }
    }

    private long getFileBlock(final long bgzfOffset) {
        return ((bgzfOffset >> 16L) & 0xFFFFFFFFFFFFL);
    }

    private void readBytes(final byte[] bytes) {
        mFileBuffer.get(bytes);
    }

    private int readInteger() {
        return mFileBuffer.getInt();
    }

    private long readLong() {
        return mFileBuffer.getLong();
    }

    private void skipBytes(final int count) {
        mFileBuffer.position(mFileBuffer.position() + count);
    }

    private void seek(final int position) {
        mFileBuffer.position(position);
    }
}