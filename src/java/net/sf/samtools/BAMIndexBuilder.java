/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sub-license, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.samtools;

import static net.sf.samtools.util.BlockCompressedInputStream.getFileBlock;
import static net.sf.samtools.BAMIndex.MAX_BINS;

import java.util.*;

/**
 * Class for constructing BAM index files.
 * Only used by BAMIndexer
 */
class BAMIndexBuilder {

    // the bins for the current reference
    private Bin[] bins = new Bin[MAX_BINS +1]; // one extra for metaData
    private int binsSeen = 0;

    // linear index for the current reference
    private long[] index = new long[LinearIndex.MAX_LINEAR_INDEX_SIZE];
    private int largestIndexSeen = -1;

    // information in meta data
    private long firstOffset = -1;
    private long lastOffset = 0;
    private int alignedRecords = 0;
    private int unalignedRecords = 0;
    private long noCoordinateRecords = 0;

    /**
     * Record any index information for a given BAM record
     *
     * @param rec The BAM record
     */
    public void processAlignment(final int reference, final SAMRecord rec) {

        final int alignmentStart = rec.getAlignmentStart();

        if (alignmentStart == SAMRecord.NO_ALIGNMENT_START) {
            noCoordinateRecords++;
            return; // do nothing for records without coordinates, but count them
        }

        if (rec.getReadUnmappedFlag()) {
            unalignedRecords++;
        } else {
            alignedRecords++;
        }

        // process bins

        final Integer binNumber = rec.getIndexingBin();
        final int binNum = binNumber == null ? rec.computeIndexingBin(): binNumber;

        // is there a bin already represented for this index?  if not, add one
        final Bin bin;
        if (bins[binNum] != null) {
            bin = bins[binNum];
        } else {
            bin = new Bin(reference, binNum);
            bins[binNum] = bin;
            binsSeen++;
        }

        // process chunks

        final SAMFileSource source = rec.getFileSource();
        if (source == null) {
            throw new SAMException("No source for BAM Record " + rec);
        }
        final BAMFileSpan newSpan = (BAMFileSpan) source.getFilePointer();
        Chunk newChunk = newSpan.getChunks().get(0); // only 1 chunk in a single source span

        final long chunkStart = newChunk.getChunkStart();
        if (chunkStart < firstOffset || firstOffset == -1) {
            firstOffset = chunkStart;
        }
        final long chunkEnd = newChunk.getChunkEnd();
        if (chunkEnd > lastOffset) {
            lastOffset = chunkEnd;
        }

        List<Chunk> oldChunks = bin.getChunkList();
        if (oldChunks == null) {
            oldChunks = new ArrayList<Chunk>();
            bin.setChunkList(oldChunks);
            bin.setLastChunk(newChunk);
            oldChunks.add(newChunk);

        } else {
            final Chunk lastChunk = bin.getLastChunk();

            // Coalesce chunks that are in adjacent file blocks.
            final long lastFileBlock = getFileBlock(lastChunk.getChunkEnd());
            final long chunkFileBlock = getFileBlock(chunkStart);
            if (chunkFileBlock - lastFileBlock <= 1) {
                lastChunk.setChunkEnd(chunkEnd);  // coalesced
            } else {
                oldChunks.add(newChunk);
                bin.setLastChunk(newChunk);
            }
        }

        // process linear index

        // the smallest file offset that appears in the 16k window for this bin
        final int alignmentEnd = rec.getAlignmentEnd();
        int startWindow = LinearIndex.convertToLinearIndexOffset(alignmentStart); // the 16k window
        final int endWindow;

        if (alignmentEnd == SAMRecord.NO_ALIGNMENT_START) {   // assume alignment uses one position
            // Next line for C (samtools index) compatibility. Differs only when on a window boundary
            startWindow = LinearIndex.convertToLinearIndexOffset(alignmentStart - 1);
            endWindow = startWindow;
        } else {
            endWindow = LinearIndex.convertToLinearIndexOffset(alignmentEnd);
        }

        if (endWindow > largestIndexSeen) {
            largestIndexSeen = endWindow;
        }

        // set linear index at every 16K window that this alignment overlaps
        for (int win = startWindow; win <= endWindow; win++) {
            if (index[win] == 0 || chunkStart < index[win]) {
                index[win] = chunkStart;
            }
        }
    }

    /**
     * Processing after all alignments of a reference have already been processed
     */
    public BAMIndexContent processReference(int reference) {

        // process bins
        if (binsSeen == 0) return null;  // no bins for this reference
        
        // process chunks, adding an extra bin for the meta data so that the count of n_bins comes out right
        final Bin metaBin = new Bin(reference, MAX_BINS);
        final List<Chunk> metaData = getMetaDataChunks();
        metaBin.setChunkList(metaData);
        bins[MAX_BINS]= metaBin;

        // process linear index
        // linear index will only be as long as the largest index seen
        final long[] newIndex = new long[largestIndexSeen + 1]; // in java1.6 Arrays.copyOf(index, largestIndexSeen + 1);

        // C (samtools index) also fills in intermediate 0's with values.  This seems unnecessary, but safe
        long lastNonZeroOffset = 0;
        for (int i = 0; i <= largestIndexSeen; i++) {
            if (index[i] == 0) {
                index[i] = lastNonZeroOffset; // not necessary, but C (samtools index) does this
            } else {
                lastNonZeroOffset = index[i];
            }
            newIndex[i] = index[i];
        }

        final LinearIndex linearIndex = new LinearIndex(reference, 0, newIndex);

        return new BAMIndexContent(reference, bins, binsSeen + 1, metaData, linearIndex);
    }

    /** @return the count of records with no coordinate positions */
    public long finish(){
        return noCoordinateRecords;
    }

    /**  reinitialize all data structures when the reference changes */
    void startNewReference() {
        if (binsSeen > 0){
            Arrays.fill(bins, null);
            Arrays.fill (index, 0);
        }
        binsSeen = 0;
        largestIndexSeen = -1;
        firstOffset = -1;
        lastOffset = 0;
        alignedRecords = 0;
        unalignedRecords = 0;
    }


    private List<Chunk> getMetaDataChunks() {
        // An extra bin #37450  (MAX_BINS) with 2 chunks of extra meta information
        //       offset_begin: offset_end
        //       n_mapped: n_unmapped
        final List<Chunk> metaChunkList = new ArrayList<Chunk>(1);
        metaChunkList.add(new Chunk(firstOffset == -1 ? 0 : firstOffset, lastOffset));
        metaChunkList.add(new Chunk(alignedRecords, unalignedRecords));
        return metaChunkList;
    }

    private void verbose(String message) {
        boolean verbose = true;
        if (verbose) {
            System.out.println("BAMIndexBuilder: " + message);
        }
    }
}