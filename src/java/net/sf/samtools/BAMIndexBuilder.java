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
 * Class for constructing BAM index files
 */
public class BAMIndexBuilder {

    // the bins for the current reference
    private Bin[] bins = new Bin[MAX_BINS];
    private boolean binsSeen = false;

    // linear index for the current bin
    private long[] index = new long[MAX_BINS];
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
    public void processAlignment(final SAMRecord rec) {

        if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
            noCoordinateRecords++;
            return; //  null;  // do nothing for records without coordinates, but count them
        }

        if (rec.getReadUnmappedFlag()) {
            unalignedRecords++;
        } else {
            alignedRecords++;
        }

        // process bins
        final int binNumber;
        if (rec.getIndexingBin() == null) {
            binNumber = rec.computeIndexingBin();
        } else {
            binNumber = rec.getIndexingBin();
        }

        // is there a bin already represented for this index?  if not, add it
        final Bin bin;
        if (bins[binNumber] != null) {
            bin = bins[binNumber];
        } else {
            final int reference = rec.getReferenceIndex();
            bin = new Bin(reference, binNumber);
            bins[binNumber] = bin;
            binsSeen = true;
        }

        // add chunk information from this record to the bin

        final BAMFileSpan newSpan;
        if (rec.getFileSource() == null) {
            throw new SAMException("No source for BAM Record " + rec);
        } else {
            newSpan = (BAMFileSpan) rec.getFileSource().getFilePointer();
        }

        List<Chunk> oldChunks = bin.getChunkList();
        boolean firstChunkList = false;
        if (oldChunks == null) {
            oldChunks = new ArrayList<Chunk>();
            bin.setChunkList(oldChunks);
            firstChunkList = true;
        }
        for (Chunk newChunk : newSpan.getChunks()) {
            // note, there should only be one newChunk
            // optimize chunkList as we go
            final long newFirstOffset = newChunk.getChunkStart();
            if (newFirstOffset < firstOffset || firstOffset == -1){
                firstOffset = newFirstOffset;
            }
            final long newLastOffset = newChunk.getChunkEnd();
            if (newLastOffset > lastOffset){
                lastOffset = newLastOffset;
            }
            simpleOptimizeChunkList(firstChunkList, oldChunks, newChunk, binNumber);
        }

        // add linear index information
        // the smallest file offset that appears in the 16k window for this bin
        final long iOffset = newSpan.getFirstOffset();
        final int alignmentStart = rec.getAlignmentStart();
        int alignmentEnd = rec.getAlignmentEnd();
        int startWindow = LinearIndex.convertToLinearIndexOffset(alignmentStart); // the 16k window
        int endWindow = LinearIndex.convertToLinearIndexOffset(alignmentEnd);

        if (alignmentEnd == 0) {
            if (alignmentStart != 0){
                startWindow = LinearIndex.convertToLinearIndexOffset(alignmentStart - 1);
            }
            endWindow = startWindow;  // assume alignment uses one position
        }

        /* almost, but not quite right
          if (alignmentEnd == 0) {
              alignmentEnd = alignmentStart + 1;
              endWindow = LinearIndex.convertToLinearIndexOffset(alignmentEnd);
          }
        */

        /* diagnostic
        if (reference == 1){
            final long endPos = newSpan.toCoordinateArray()[1];
            verbose("Ref " + reference + " " + rec.getReadName() + " bin=" + binNumber + " startPos=" +
              iOffset + "(" + Long.toString(iOffset, 16) + "x)" +
              "endPos=" + endPos + "(" + Long.toString(endPos, 16) + "x)" +
              " sA=" + alignmentStart + "(" + Long.toString(alignmentStart,16) + "x)" +
              " eA=" + alignmentEnd + "(" + Long.toString(alignmentEnd,16) + "x)" +
              "     sWin " + startWindow + " eWin " + endWindow );
        }
        */

        // set linear index at every 16K window that this alignment overlaps
        for (int win = startWindow; win <= endWindow; win++) {
            if (index[win] != 0) {
                if (iOffset < index[win]) {
                    index[win] = iOffset;
                }
            } else {
                if (win > largestIndexSeen){
                    largestIndexSeen = win;
                }
                index[win] = iOffset;
            }
        }
    }

    /**
     * Processing after all alignments of a reference have already been processed
     */
    public BAMIndexContent processReference(int reference) {

        // process bins
        if (!binsSeen) return null;  // no bins for this reference

        List<Bin> binList = new ArrayList <Bin> ();

        // process chunks
        for (Bin bin : bins) {
            if (bin == null) continue;
            binList.add(bin);
        }

        // for c compatibility, add an extra bin 37450 with 2 chunks of extra meta information
        //       offset_begin: offset_end
        //       n_mapped: n_unmapped
        final Bin metaBin = new Bin(reference, MAX_BINS);
        List<Chunk> metaChunkList = new ArrayList<Chunk>(1);
        metaBin.setChunkList(metaChunkList);
        metaChunkList.add(new Chunk(firstOffset, lastOffset));
        metaChunkList.add(new Chunk(alignedRecords, unalignedRecords));
        binList.add(metaBin);

        // process linear index
        final LinearIndex linearIndex = computeLinearIndex(reference);

        return new BAMIndexContent(reference, binList, linearIndex);
    }

    /**  reinitialize all data structures when the reference changes */
    public void startNewReference() {
        if (binsSeen){
            Arrays.fill(bins, null);
            Arrays.fill (index, 0);
        }
        binsSeen = false;
        largestIndexSeen = -1;
        firstOffset = -1;
        lastOffset = 0;
        alignedRecords = 0;
        unalignedRecords = 0;
    }

    /** @return the count of records with no coordinate positions */
    public long finish(){
        return noCoordinateRecords;
    }

    /**
     * This is a simpler version of optimizeChunkList found in AbstractBamFileIndex.
     */
    private void simpleOptimizeChunkList(final boolean firstChunkList, final List<Chunk> chunks, final Chunk newChunk, final int binNumber) {
        if (firstChunkList|| binNumber == MAX_BINS) {
            chunks.add(newChunk);
            return;
        }
        Chunk lastChunk = null;
        for (final Chunk chunk : chunks) {  // todo - performance maybe each bin should keep its lastChunk ?
            lastChunk = chunk;
        }
        // Coalesce chunks that are in adjacent file blocks.
        final long lastFileBlock = getFileBlock(lastChunk.getChunkEnd());
        final long chunkFileBlock = getFileBlock(newChunk.getChunkStart());
        if (chunkFileBlock - lastFileBlock <= 1) {     // probably always true
            lastChunk.setChunkEnd(newChunk.getChunkEnd());
        } else {
            chunks.add(newChunk);
        }
    }

    /**
     * The linear index stored can be smaller than full array of all possible windows
     */
    private LinearIndex computeLinearIndex(int reference) {
        final LinearIndex linearIndex;
        if (largestIndexSeen == -1) {
            // skip linear index for this reference
            linearIndex = null;
        } else {
            // linear index will be as long as the largest index seen
            long[] newIndex = new long[largestIndexSeen + 1]; // in java1.6 Arrays.copyOf(index, largestIndexSeen + 1);
            for (int i = 0; i <= largestIndexSeen; i++) {
                newIndex[i] = index[i];
            }

            // c samtools index also fills in intermediate 0's with values.  This seems unnecessary, but safe
            long lastOffset = 0;
            for (int i=0; i <= largestIndexSeen; i++){
                if (newIndex[i] == 0){
                    newIndex[i] = lastOffset;
                } else {
                    lastOffset = newIndex[i];
                }
            }
            linearIndex = new LinearIndex(reference, 0, newIndex);
        }
        return linearIndex;
    }

    private void verbose(String message) {
        boolean verbose = true;
        if (verbose) {
            System.out.println("BAMIndexBuilder: " + message);
        }
    }
}