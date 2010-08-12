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

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.io.Serializable;

/**
 * A interface representing a collection of (possibly) discontinuous segments in the
 * BAM file, possibly representing the results of an index query.
 */
public interface SAMFileSpan extends Cloneable {
    /**
     * Gets a pointer over the data immediately following this span.
     * @return The a pointer to data immediately following this span.
     */
    public SAMFileSpan getContentsFollowing();

    /**
     * Remove all pointers in this file span before the given file span starts.
     * @param fileSpan The filespan before which to eliminate.
     * @return The portion of the chunk list after the given chunk.
     */
    public SAMFileSpan removeContentsBefore(final SAMFileSpan fileSpan);

    /**
     * Does this file span point to any data, or is it completely empty?
     * @return True if the file span is empty, false otherwise.
     */
    public boolean isEmpty();
}

/**
 * An ordered list of chunks, capable of representing a set of discontiguous
 * regions in the BAM file.  FileSpans are mutable within the package, but perceived
 * as immutable outside the package.
 *
 * Some operations on FileSpans assume that the spans are sorted.  In these cases,
 * sort order will be validated.
 *
 * @author mhanna
 * @version 0.1
 */
class BAMFileSpan implements SAMFileSpan, Serializable {
    private static final long serialVersionUID = 1L;    

    /**
     * The constituent chunks of this list.
     */
    private final List<Chunk> chunks;

    /**
     * Create a new empty list of chunks.
     */
    protected BAMFileSpan() {
        this.chunks = new ArrayList<Chunk>();
    }

    /**
     * Convenience constructor to construct a BAM file span from
     * a single chunk.
     * @param chunk Chunk to use as the sole region in this span.
     */
    protected BAMFileSpan(final Chunk chunk) {
        this.chunks = new ArrayList<Chunk>();
        chunks.add(chunk);
    }

    /**
     * Create a new chunk list from the given list of chunks.
     * @param chunks Constituent chunks.
     */
    protected BAMFileSpan(final List<Chunk> chunks) {
        this.chunks = new ArrayList<Chunk>(chunks);
    }

    /**
     * Does this chunk list map to any position within the BAM file?
     * @return True iff the ChunkList points to any data within the BAM.
     */
    public boolean isEmpty() {
        return chunks.isEmpty();    
    }

    /**
     * Deep clone the given chunk list.
     * @return A copy of the chunk list.
     */
    public BAMFileSpan clone() {
        BAMFileSpan clone = new BAMFileSpan();
        for(Chunk chunk: chunks)
            clone.chunks.add(chunk.clone());
        return clone;
    }

    /**
     * Remove all chunks in this file span before the given file span starts.
     * If a chunk in the chunk list starts before and ends after the given
     * chunk, the first portion of the chunk will be deleted.
     * @param fileSpan The filespan before which to eliminate.
     * @return The portion of the chunk list after the given chunk.
     */
    public SAMFileSpan removeContentsBefore(final SAMFileSpan fileSpan) {
        if(fileSpan == null)
            return clone();

        if(!(fileSpan instanceof BAMFileSpan))
            throw new SAMException("Unable to compare ");

        BAMFileSpan bamFileSpan = (BAMFileSpan)fileSpan;

        if(bamFileSpan.isEmpty())
            return clone();

        validateSorted();

        BAMFileSpan trimmedChunkList = new BAMFileSpan();
        for(Chunk chunkToTrim: chunks) {
            if(chunkToTrim.getChunkEnd() > chunkToTrim.getChunkStart()) {
                if(chunkToTrim.getChunkStart() >= bamFileSpan.chunks.get(0).getChunkStart()) {
                    // This chunk from the list is completely beyond the start of the filtering chunk.
                    trimmedChunkList.add(chunkToTrim.clone());
                }
                else {
                    // This chunk from the list partially overlaps the filtering chunk and must be trimmed.                    
                    trimmedChunkList.add(new Chunk(bamFileSpan.chunks.get(0).getChunkStart(),chunkToTrim.getChunkEnd()));
                }
            }
        }
        return trimmedChunkList;
    }

    /**
     * Gets a file span over the data immediately following this span.
     * @return The a pointer to data immediately following this span.
     */
    public SAMFileSpan getContentsFollowing() {
        if(chunks.isEmpty())
            throw new SAMException("Unable to get the file pointer following this one: no data present.");
        validateSorted();
        return new BAMFileSpan(new Chunk(chunks.get(chunks.size()-1).getChunkEnd(),Long.MAX_VALUE));
    }

    /**
     * Merge one span into another
     *
     * @param span - span with chunks to add to this one
     */
    public void add(final BAMFileSpan span) {
        for (Chunk c : span.chunks) {
            chunks.add(c);
        }
    }

    /**
     * Adds a new chunk to this list.  Visible only within the BAm.
     * @param chunk Chunk to add.
     */
    protected void add(final Chunk chunk) {
        chunks.add(chunk);
    }
    
    /**
     * Convert the chunk list to an array of offsets, paired in [start,end) format.
     * @return Array of offsets.
     */
    protected long[] toCoordinateArray() {
        final int count = chunks.size() * 2;
        if (count == 0) {
            return null;
        }
        int index = 0;
        final long[] result = new long[count];
        for (final Chunk chunk : chunks) {
            result[index++] = chunk.getChunkStart();
            result[index++] = chunk.getChunkEnd();
        }
        return result;
    }

    /**
     * Find the first offset in the chunk list
     * @return The first offset in the span
     */
    protected long getFirstOffset() {
        long result = 0;
        if (chunks == null){
            return result;
        }
        for (final Chunk chunk : chunks) {
            return chunk.getChunkStart();
        }
        return result;
    }

    /**
     * Gets the constituent chunks stored in this span.
     * @return An unmodifiable list of chunks.
     */
    protected List<Chunk> getChunks() {
        return Collections.unmodifiableList(chunks);
    }

    /**
     * Checks that there is only a single chunk for this span and returns it.
     * @return The single chunk stored in this span
     */
    protected Chunk getSingleChunk() {
        if (chunks.size() != 1){
            throw new SAMException("Expecting a single chunk for span. Found " + chunks.size());
        }
        return chunks.get(0);
    }

    /**
     * The list of chunks is often represented as an array of
     * longs where every even-numbered index is a start coordinate
     * and every odd-numbered index is a stop coordinate.  Convert
     * from that format back to a list of chunks.
     * @param coordinateArray List of chunks to convert.
     * @return A list of chunks.
     */
    protected static SAMFileSpan toChunkList(long[] coordinateArray) {
        if(coordinateArray.length % 2 != 0)
            throw new SAMException("Data supplied does not appear to be in coordinate array format.");

        BAMFileSpan chunkList = new BAMFileSpan();
        for(int i = 0; i < coordinateArray.length; i += 2)
            chunkList.add(new Chunk(coordinateArray[i],coordinateArray[i+1]));

        chunkList.validateSorted();

        return chunkList;
    }

    /**
     * Validates the list of chunks to ensure that they appear in sorted order.
     */
    private void validateSorted() {
        for(int i = 1; i < chunks.size(); i++) {
            if(chunks.get(i).getChunkStart() < chunks.get(i-1).getChunkEnd())
                throw new SAMException(String.format("Chunk list is unsorted; chunk %s is before chunk %s",chunks.get(i-1),chunks.get(i)));
        }
    }

    /**
     * Creates a string representation of this chunk list.
     */
    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        boolean first = true;
        for(Chunk chunk: chunks) {
            if(!first) {
                builder.append(';');
                first = false;
            }
            builder.append(chunk);
        }
        return builder.toString();
    }
}
