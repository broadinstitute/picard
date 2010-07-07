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

import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.util.*;

/** Represents the contents of a bam index file for one reference.
 * A BAM index (.bai) file contains information for all references in the bam file.
 * This class describes the data present in the index file for one of these references;
 * including the bins, chunks, and linear index.
 */
class BAMIndexContent {
    /**
     * The reference sequence for the data currently loaded.
     */
    private final int mReferenceSequence;

    /**
     * A list of all bins in the above reference sequence.
     */
    private final List<Bin> mBins;

    /**
     * A mapping from bin to the chunks contained in that bin.
     */
    private final SortedMap<Bin, List<Chunk>> mBinToChunks;

    /**
     * The linear index for the reference sequence above.
     */
    private final LinearIndex mLinearIndex;

    public BAMIndexContent(final int referenceSequence, final List<Bin> bins, final SortedMap<Bin,List<Chunk>> binToChunks, final LinearIndex linearIndex) {
        this.mReferenceSequence = referenceSequence;
        this.mBins = bins;
        this.mBinToChunks = binToChunks;
        this.mLinearIndex = linearIndex;
    }

    public int getReferenceSequence() {
        return mReferenceSequence;
    }

    public boolean containsBin(final Bin bin) {
        return Collections.binarySearch(mBins,bin) >= 0;
    }

    public List<Bin> getBins() {
        return Collections.unmodifiableList(mBins);
    }

    public List<Chunk> getChunksForBin(final Bin bin) {
        if(!mBinToChunks.containsKey(bin))
            throw new SAMException("No chunks found for the given bin.");
        return Collections.unmodifiableList(mBinToChunks.get(bin));
    }

    public List<Chunk> getAllChunks() {
        List<Chunk> allChunks = new ArrayList<Chunk>();
        for(List <Chunk> moreChunks: mBinToChunks.values())
            allChunks.addAll(moreChunks);
        return Collections.unmodifiableList(allChunks);
    }

    public LinearIndex getLinearIndex() {
        return mLinearIndex;
    }

    // todo - remove - only used by Old BuildBamIndex
      /**
     * Write this content as human-readable text
     *
     * @param pw PrintWriter for text output file
     */
    public void writeText(final PrintWriter pw, boolean sortBins) {

        if (mBins == null || mBins.size() == 0) {
            writeNullTextContent(pw, mReferenceSequence);
            return;
        }

        final int size = mBins.size();
        pw.println("Reference " + mReferenceSequence + " has n_bin= " + size);

        // copy into an array so that it can be sorted
        final Bin[] bins = new Bin[size];
        if (size != 0) {
            getBins().toArray(bins);
        }
        if (sortBins) Arrays.sort(bins);  // Sort for easy text comparisons
        for (int j = 0; j < size; j++) {
            if (mBinToChunks.get(bins[j]) == null) {
                pw.println("  Ref " + mReferenceSequence + " bin " + bins[j].getBinNumber() + " has no mBinToChunks"); // remove?
                continue;
            }
            final List<Chunk> chunkList = mBinToChunks.get(bins[j]);
            if (chunkList == null) {
                pw.println("  Ref " + mReferenceSequence + " bin " + bins[j].getBinNumber() + " has no chunkList");
                continue;
            }
            if (chunkList.size() == 0){
                // continue; // todo
                pw.println();
            }
            pw.print("  Ref " + mReferenceSequence + " bin " + bins[j].getBinNumber() + " has n_chunk= " + chunkList.size());
            for (final Chunk c : chunkList) {
                pw.println("     Chunk: " + c.toString() +
                        " start: " + Long.toString(c.getChunkStart(), 16) +
                        " end: " + Long.toString(c.getChunkEnd(), 16));
            }
        }
        if (mLinearIndex == null || mLinearIndex.getIndexEntries() == null) {
            pw.println("Reference " + mReferenceSequence + " has n_intv= 0");
            return;
        }
        final long[] entries = mLinearIndex.getIndexEntries();
        final int indexStart = mLinearIndex.getIndexStart();
        // System.out.println("index start is " + indexStart);
        final int n_intv = entries.length + indexStart;
        pw.println("Reference " + mReferenceSequence + " has n_intv= " + n_intv);
        for (int k = 0; k < entries.length; k++) {
            if (entries[k] != 0) {
                pw.println("  Ref " + mReferenceSequence + " ioffset for " + (k + indexStart) + " is " + Long.toString(entries[k]));
            }
        }
    }

    /**
     * Write this content as binary output
     *
     * @param bb ByteBuffer for output file must exist and be writable
     */
    public void writeBinary(ByteBuffer bb, boolean sortBins) {
        if (mBins == null || mBins.size() == 0) {
            writeNullBinaryContent(bb);
            return;
        }

        final int size = mBins.size();
        bb.putInt(size);

        // copy into an array so that it can be sorted
        final Bin[] bins = new Bin[size];
        if (size != 0) {
            getBins().toArray(bins);
        }
        if (sortBins) Arrays.sort(bins);  // Sort for easy text comparisons
        for (int j = 0; j < size; j++) {

            bb.putInt(bins[j].getBinNumber()); // todo uint32_t vs int32_t in spec?
            if (mBinToChunks.get(bins[j]) == null) {
                bb.putInt(0);
                continue;
            }
            final List<Chunk> chunkList = mBinToChunks.get(bins[j]);
            final int n_chunk = chunkList.size();
            bb.putInt(n_chunk);
            for (final Chunk c : chunkList) {
                bb.putLong(c.getChunkStart());   // todo uint32_t vs int32_t in spec?
                bb.putLong(c.getChunkEnd());     // todo uint32_t vs int32_t in spec?
            }
        }
        final long[] entries = mLinearIndex == null ? null : mLinearIndex.getIndexEntries();
        final int indexStart = mLinearIndex == null ? 0 : mLinearIndex.getIndexStart();
        final int n_intv = entries == null ? indexStart : entries.length + indexStart; // +1;
        bb.putInt(n_intv);
        if (entries == null){
            return;
        }
        // System.out.println("index start is " + indexStart);
        for (int i = 0; i < indexStart; i++) {
            bb.putLong(0);          // todo uint32_t vs int32_t in spec?
        }
        for (int k = 0; k < entries.length; k++) {
            bb.putLong(entries[k]); // todo uint32_t vs int32_t in spec?
        }
    }

    static void writeNullTextContent(PrintWriter pw, int reference) {
        pw.println("Reference " + reference + " has n_bin=0");
        pw.println("Reference " + reference + " has n_intv=0");
    }

    static void writeNullBinaryContent(ByteBuffer bb) {
        bb.putInt(0);  // 0 bins
        bb.putInt(0);  // 0 intv
    }
}
