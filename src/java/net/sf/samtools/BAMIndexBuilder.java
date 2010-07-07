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

import java.util.*;

import static net.sf.samtools.BAMIndex.MAX_BINS;

/**
 * Class for constructing BAM index files
 */
public class BAMIndexBuilder {

    private int currentReference = 0;

    // the bins for the current reference
    private BitSet binsSeen = new BitSet(MAX_BINS);
    private Bin[] bins = new Bin[MAX_BINS];

    // linear index for the current bin
    private BitSet indexSeen = new BitSet(MAX_BINS);
    private long[] index = new long[MAX_BINS];

    // information in meta data
    private int alignedRecords = 0;
    private int unalignedRecords = 0;
    private int noCoordinateRecords = 0;

    /**
     * Record any index information for a given BAM record
     *
     * @param rec The BAM record
     * @return the last finished reference, if done with it, otherwise null;
     */
    public BAMIndexContent processAlignment(final SAMRecord rec) {

        BAMIndexContent result = null;  // non-null only when we've finished a reference

        final BAMFileSpan newSpan;
        if (rec.getFileSource() == null) {
            throw new SAMException("No source for BAM Record " + rec);
        } else {
            newSpan = (BAMFileSpan) rec.getFileSource().getFilePointer();
        }

        if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
            noCoordinateRecords++;
            return result;  // do nothing for records without coordinates, but count them
        }

        final int binNumber;
        if (rec.getIndexingBin() == null) {
            binNumber = rec.computeIndexingBin();
        } else {
            binNumber = rec.getIndexingBin();
        }

        final int reference = rec.getReferenceIndex();
        if (reference != currentReference) {
            result = processReference(currentReference);
            currentReference = reference;
            startNewReference();
        }

        if (rec.getReadUnmappedFlag()) {
            unalignedRecords++;
        } else {
            alignedRecords++;
        }

        // is there a bin already represented for this index?  if not, add it
        final Bin bin;
        if (binsSeen.get(binNumber)) {
            bin = bins[binNumber];
        } else {
            binsSeen.set(binNumber);
            bin = new Bin(reference, binNumber);
            bins[binNumber] = bin;
        }

        // add chunk information from this record to the bin

        List<Chunk> chunks = bin.getChunkList();
        if (chunks == null) {
            chunks = new ArrayList<Chunk>();
            bin.setChunkList(chunks);
        }
        for (Chunk c : newSpan.getChunks()) {
            // there should only be 1 chunk c in this loop
            // chunks should have 0 or 1 chunk in its list

            /* diagnostic
           if (reference == 1){
               verbose(++ totalRecords + " Chunks for bin " + Integer.toString(binNumber) + " " +
                       Long.toString(c.getChunkStart()) + "(" + Long.toString(c.getChunkStart(),16) + "x)" + " : " +
                       Long.toString(c.getChunkEnd()) + "(" + Long.toString(c.getChunkEnd(),16) + "x)" +
                       "SA=" + (rec.getAlignmentStart() - 1) );
           }  */
            // optimize chunkList as we go
            simpleOptimizeChunkList(chunks, c, binNumber);
        }

        // add linear index information
        // the smallest file offset that appears in the 16k window for this bin
        final long iOffset = newSpan.toCoordinateArray()[0];
        final int alignmentStart = rec.getAlignmentStart();
        int alignmentEnd = rec.getAlignmentEnd();
        int startWindow = LinearIndex.convertToLinearIndexOffset(alignmentStart); // the 16k window
        int endWindow = LinearIndex.convertToLinearIndexOffset(alignmentEnd);

        if (alignmentEnd == 0) {
            if (alignmentStart != 0)
                startWindow = LinearIndex.convertToLinearIndexOffset(alignmentStart - 1);
            /* if (alignmentStart == 1){
                verbose("Forcing endWin to 18! *********"); // todo just an experiment! Fixes test case
                endWindow=18;
            } else  */
            endWindow = startWindow;  // assume alignment uses one position
        }

        /* almost, but not quite
          if (alignmentEnd == 0) {
              alignmentEnd = alignmentStart + 1;
              endWindow = LinearIndex.convertToLinearIndexOffset(alignmentEnd);
          }
        */

        /* diagnostic if (startWindow != endWindow && (endWindow -startWindow != 1)){
        if (reference == 1){
            // final Cigar cigar = rec.getCigar();
            final long endPos = newSpan.toCoordinateArray()[1];
            verbose("Ref " + reference + " " + rec.getReadName() + " bin=" + binNumber + " startPos=" +
              iOffset + "(" + Long.toString(iOffset, 16) + "x)" +
                  "endPos=" + endPos + "(" + Long.toString(endPos, 16) + "x)" + " sA=" +
              alignmentStart + "(" + Long.toString(alignmentStart,16) + "x)" +
              " eA=" + alignmentEnd + "(" + Long.toString(alignmentEnd,16) + "x)" +
              "     sWin " + startWindow + " eWin " + endWindow );
            //  " rec.getReadUnmappedFlag=" + rec.getReadUnmappedFlag());
            // " getCigar().getReferenceLength() " + cigar.getReferenceLength() + " cigar=" + cigar );
        }
        */

        // set linear index at every 16K window that this alignment overlaps
        for (int win = startWindow; win <= endWindow; win++) { // set linear index at every 16K window

            // Has this window already been seen? If so, won't replace existing index,
            // since first ref encountered is always earliest (since coordinate sorted)
            if (indexSeen.get(win)) {
                if (iOffset < index[win]) {
                    index[win] = iOffset;
                }
            } else {
                indexSeen.set(win);
                index[win] = iOffset;
            }
        }
        return result;
    }

    /**
     * Processing after all alignments of a reference have already been processed
     */
    public BAMIndexContent processReference(int reference) {

        // process bins
        if (binsSeen.isEmpty()) return null;  // no bins for this reference

        List<Bin> binList = new ArrayList <Bin> ();

        // process chunks
        final SortedMap<Bin, List<Chunk>> binToChunks = new TreeMap<Bin, List<Chunk>>();
        long firstOffset = -1;
        long lastOffset = 0;
        for (Bin bin : bins) {
            if (bin == null) continue;
            List<Chunk> chunkList = bin.getChunkList();
            if (chunkList != null){
                // calculate first and last offset
                long newFirstOffset = chunkList.get(0).getChunkStart();
                if (firstOffset == -1  || newFirstOffset < firstOffset){
                    firstOffset = newFirstOffset;
                }
                lastOffset = chunkList.get(chunkList.size()-1).getChunkEnd();

                bin.setChunkList(chunkList);
                binToChunks.put(bin, chunkList);
            }
            binList.add(bin);
        }
        // for c compatibility, add an extra bin 37450 with 2 chunks of extra meta information
        //       offset_begin: offset_end
        //       n_mapped: n_unmapped
        final Bin metaBin = new Bin(reference, MAX_BINS);
        List<Chunk> metaChunkList = new ArrayList<Chunk>(1);
        metaBin.setChunkList(metaChunkList);
        binToChunks.put(metaBin, metaChunkList);
        metaChunkList.add(new Chunk(firstOffset == -1 ? 0 : firstOffset, lastOffset)) ;
        metaChunkList.add(new Chunk(alignedRecords, unalignedRecords)) ;
        binList.add(metaBin);

        // process linear index
        final LinearIndex linearIndex = computeLinearIndex(reference);

        return new BAMIndexContent(reference, binList, binToChunks, linearIndex);
    }


    /**  reinitialize all data structures when the reference changes */
    public void startNewReference() {
        binsSeen.clear();
        bins = new Bin[MAX_BINS];
        indexSeen.clear();
        index = new long[MAX_BINS];
        alignedRecords = 0;
        unalignedRecords = 0;
    }

    // This is a simpler version of optimizeChunkList found in AbstractBamFileIndex
    // It can be simpler because we are processing chunks in bam file, coordinate sorted order
    private void simpleOptimizeChunkList(final List<Chunk> chunks, final Chunk newChunk, final int binNumber) {
        if (chunks.size() == 0  || binNumber == MAX_BINS) {
            chunks.add(newChunk);
            return;
        }
        Chunk lastChunk = null;
        for (final Chunk chunk : chunks) {
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
     * transform the sparse representation of linear index to an array representation
     */
    private LinearIndex computeLinearIndex(int reference) {
        final LinearIndex linearIndex;
        if (indexSeen.isEmpty()) {
            // skip linear index for this reference
            linearIndex = null;
        } else {
            // get the max window in the list; linear index will be this long
             // get the max window in the list; linear index will be this long

            int maxWindow = indexSeen.length();
            long[] newIndex = Arrays.copyOf(index, maxWindow);

            // c index also fills in intermediate 0's with values.  This seems incorrect todo
            long lastOffset = 0;
            for (int i=0; i < maxWindow; i++){
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