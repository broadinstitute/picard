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

import net.sf.samtools.seekablestream.SeekableStream;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * A class for reading BAM file indices, hitting the disk once per query.
 */
class DiskBasedBAMFileIndex extends AbstractBAMFileIndex
{
    DiskBasedBAMFileIndex(final File file, SAMSequenceDictionary dictionary) {
        super(file, dictionary);
    }

    DiskBasedBAMFileIndex(final SeekableStream stream, SAMSequenceDictionary dictionary) {
        super(stream, dictionary);
    }

    DiskBasedBAMFileIndex(final File file, SAMSequenceDictionary dictionary, boolean useMemoryMapping) {
        super(file, dictionary, useMemoryMapping);
    }

    /**
     * Get list of regions of BAM file that may contain SAMRecords for the given range
     * @param referenceIndex sequence of desired SAMRecords
     * @param startPos 1-based start of the desired interval, inclusive
     * @param endPos 1-based end of the desired interval, inclusive
     * @return array of pairs of virtual file positions.  Each pair is the first and last
     * virtual file position in a range that can be scanned to find SAMRecords that overlap the given
     * positions. The last position in each pair is a virtual file pointer to the first SAMRecord beyond
     * the range that may contain the indicated SAMRecords.
     */
    public BAMFileSpan getSpanOverlapping(final int referenceIndex, final int startPos, final int endPos) {
        BAMIndexContent queryResults = query(referenceIndex,startPos,endPos);

        if(queryResults == null)
            return null;

        List<Chunk> chunkList = new ArrayList<Chunk>();
        for(Chunk chunk: queryResults.getAllChunks())
            chunkList.add(chunk.clone());
        chunkList = optimizeChunkList(chunkList,queryResults.getLinearIndex().getMinimumOffset(startPos));
        return new BAMFileSpan(chunkList);
    }

     protected BAMIndexContent getQueryResults(int reference){
         throw new UnsupportedOperationException();
         // todo: there ought to be a way to support this using the first startPos for the reference and the last
         // return query(reference, 1, -1);
         // If this were implemented, BAMIndexer.createAndWriteIndex could extend DiskBasedBAMFileIndex -or- CachingBAMFileIndex
    }
}
