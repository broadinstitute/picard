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
import java.util.*;

/**
 * Class for reading BAM file indices, caching each contig as it's loaded and
 * dropping values when the next contig is loaded.
 */
class CachingBAMFileIndex extends AbstractBAMFileIndex implements BrowseableBAMIndex
{
    private Integer mLastReferenceRetrieved = null;
    private WeakHashMap<Integer,BAMIndexContent> mQueriesByReference = new WeakHashMap<Integer,BAMIndexContent>();

    public CachingBAMFileIndex(final File file, SAMSequenceDictionary dictionary) {
        super(file, dictionary);
    }

    public CachingBAMFileIndex(final SeekableStream stream, SAMSequenceDictionary dictionary) {
        super(stream, dictionary);
    }

    public CachingBAMFileIndex(final File file, SAMSequenceDictionary dictionary, boolean useMemoryMapping) {
        super(file, dictionary, useMemoryMapping);
    }

    /**
     * Get list of regions of BAM file that may contain SAMRecords for the given range
     * @param referenceIndex sequence of desired SAMRecords
     * @param startPos 1-based start of the desired interval, inclusive
     * @param endPos 1-based end of the desired interval, inclusive
     * @return the virtual file position.  Each pair is the first and last virtual file position
     *         in a range that can be scanned to find SAMRecords that overlap the given positions.
     */
    public BAMFileSpan getSpanOverlapping(final int referenceIndex, final int startPos, final int endPos) {
        BAMIndexContent queryResults = getQueryResults(referenceIndex);

        if(queryResults == null)
            return null;

        BinList overlappingBins = getBinsOverlapping(referenceIndex,startPos,endPos);

        // System.out.println("# Sequence target TID: " + referenceIndex);
        List<Bin> bins = new ArrayList<Bin>();
        for(Bin bin: queryResults.getBins()) {
            if (overlappingBins.getBins().get(bin.getBinNumber()))
                bins.add(bin);
        }

        if (bins.isEmpty()) {
            return null;
        }

        List<Chunk> chunkList = new ArrayList<Chunk>();
        for(Bin bin: bins) {
            for(Chunk chunk: bin.getChunkList())
                chunkList.add(chunk.clone());
        }

        if (chunkList.isEmpty()) {
            return null;
        }

        chunkList = optimizeChunkList(chunkList,queryResults.getLinearIndex().getMinimumOffset(startPos));
        return new BAMFileSpan(chunkList);
    }

    /**
     * Get a list of bins in the BAM file that may contain SAMRecords for the given range.
     * @param referenceIndex sequence of desired SAMRecords
     * @param startPos 1-based start of the desired interval, inclusive
     * @param endPos 1-based end of the desired interval, inclusive
     * @return a list of bins that contain relevant data.
     */
    public BinList getBinsOverlapping(final int referenceIndex, final int startPos, final int endPos) {
        final BitSet regionBins = regionToBins(startPos,endPos);
        if (regionBins == null) {
            return null;
        }
        return new BinList(referenceIndex,regionBins);        
    }

    /**
     * Perform an overlapping query of all bins bounding the given location.
     * @param bin The bin over which to perform an overlapping query.
     * @return The file pointers
     */
    public BAMFileSpan getSpanOverlapping(final Bin bin) {
        if(bin == null)
            return null;

        final int referenceSequence = bin.getReferenceSequence();
        BAMIndexContent indexQuery = getQueryResults(referenceSequence);

        if(indexQuery == null)
            return null;

        final int binLevel = getLevelForBin(bin);
        final int firstLocusInBin = getFirstLocusInBin(bin);

        // Add the specified bin to the tree if it exists.
        List<Bin> binTree = new ArrayList<Bin>();
        if(indexQuery.containsBin(bin))
            binTree.add(indexQuery.getBins().getBin(bin.getBinNumber()));

        int currentBinLevel = binLevel;
        while(--currentBinLevel >= 0) {
            final int binStart = getFirstBinInLevel(currentBinLevel);
            final int binWidth = getMaxAddressibleGenomicLocation()/getLevelSize(currentBinLevel);
            final int binNumber = firstLocusInBin/binWidth + binStart;
            Bin parentBin = indexQuery.getBins().getBin(binNumber);
            if(parentBin != null && indexQuery.containsBin(parentBin))
                binTree.add(parentBin);
        }

        List<Chunk> chunkList = new ArrayList<Chunk>();
        for(Bin coveringBin: binTree) {
            for(Chunk chunk: coveringBin.getChunkList())
                chunkList.add(chunk.clone());
        }

        final int start = getFirstLocusInBin(bin);
        chunkList = optimizeChunkList(chunkList,indexQuery.getLinearIndex().getMinimumOffset(start));
        return new BAMFileSpan(chunkList);
    }

    /**
     * Looks up the cached BAM query results if they're still in the cache and not expired.  Otherwise,
     * retrieves the cache results from disk.
     * @param referenceIndex The reference to load.  CachingBAMFileIndex only stores index data for entire references. 
     * @return The index information for this reference.
     */
    protected BAMIndexContent getQueryResults(final int referenceIndex) {
        // WeakHashMap is a bit weird in that its lookups are done via equals() equality, but expirations must be
        // handled via == equality.  This implementation jumps through a few hoops to make sure that == equality still
        // holds even in the context of boxing/unboxing.

        // If this query is for the same reference index as the last query, return it.
        if(mLastReferenceRetrieved!=null && mLastReferenceRetrieved == referenceIndex)
            return mQueriesByReference.get(referenceIndex);

        // If not, check to see whether it's available in the cache.
        BAMIndexContent queryResults = mQueriesByReference.get(referenceIndex);
        if(queryResults != null) {
            mLastReferenceRetrieved = referenceIndex;
            mQueriesByReference.put(referenceIndex,queryResults);
            return queryResults;
        }

        // If not in the cache, attempt to load it from disk.
        queryResults = query(referenceIndex,1,-1);
        if(queryResults != null) {
            mLastReferenceRetrieved = referenceIndex;
            mQueriesByReference.put(referenceIndex,queryResults);
            return queryResults;
        }

        // Not even available on disk.
        return null;
    }
}
