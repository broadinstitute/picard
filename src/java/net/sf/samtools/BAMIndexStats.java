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

import net.sf.samtools.util.BlockCompressedFilePointerUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Metadata about the bam index contained within the bam index.
 * One instance created per index file.
 */
public class BAMIndexStats {

    // information for the entire index.
    // stored at the end of the index
    private long noCoordinateRecords = 0;

    // information for each reference.
    // stored in two chunks in bin # MAX_BINS
    private long firstOffset = -1;
    private long lastOffset = 0;
    private int alignedRecords = 0;
    private int unAlignedRecords = 0;  // unmapped, but associated with this reference

    /** construct one instance for each index generated */
    BAMIndexStats(){
       noCoordinateRecords = 0;
        newReference();
    }

    /** Call for each new reference sequence encountered */
    void newReference() {
        firstOffset = -1;
        lastOffset = 0;
        alignedRecords = 0;
        unAlignedRecords = 0;
    }

    /** Metadata available for each reference */
     List<Chunk> getMetaDataChunks() {
        // An extra bin #37450  (MAX_BINS) with 2 chunks of extra meta information
        //       offset_begin: offset_end
        //       n_mapped: n_unmapped
        final List<Chunk> metaChunkList = new ArrayList<Chunk>(1);
        metaChunkList.add(new Chunk(firstOffset == -1 ? 0 : firstOffset, lastOffset));
        metaChunkList.add(new Chunk(alignedRecords, unAlignedRecords));
        return metaChunkList;
     }

    /** Call whenever a reference with no coordinate information is encountered in the bam file */
    void  incrementNoCoordinateRecordCount(){
        noCoordinateRecords++;
    }

    /** @return the count of records with no coordinate information in the bam file */
    long getNoCoordinateRecordCount(){
        return noCoordinateRecords;
    }

    /** Call whenever an aligned record is encountered in the bam file to record the smallest file offset */
    void setFirstOffsetIfSmaller(long newFirstOffset){
        if (BlockCompressedFilePointerUtil.compare(newFirstOffset, firstOffset) < 1  || firstOffset == -1) {
            this.firstOffset = newFirstOffset;
        }
    }

    /** Call whenever an aligned record is encountered in the bam file to record the largest file offset */
     void setLastOffsetIfLarger(long newLastOffset){
        if (BlockCompressedFilePointerUtil.compare(lastOffset, newLastOffset) < 1) {
            this.lastOffset = newLastOffset;
        }
    }

    /** Call once per record whenever an aligned record is encountered in the bam file */
    void  incrementAlignedRecordCount(){
        alignedRecords++;
    }

    /** Call once per record whenever a reference that is unaligned is encountered in the bam file */
    void  incrementUnAlignedRecordCount(){
        unAlignedRecords++;
    }
}
