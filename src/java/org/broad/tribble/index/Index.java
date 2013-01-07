/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble.index;

import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.IOException;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

/**
 * Interface for all index implementations.  
 */
public interface Index {
    /**
     * 
     * @param chr the chromosome
     * @param start the start position
     * @param end the end position
     * @return a list of blocks that contain the specified interval.
     * @throws IllegalArgumentException of chr isn't part of this index
     */
    List<Block> getBlocks(String chr, int start, int end);

    /**
     * does this index represent an up-to-date version
     * @return true if the index is up to date, false otherwise
     */
    public boolean isCurrentVersion();

    /**
     * get a list of the sequence names we've seen during indexing, in order
     * @return a LinkedHashSet, which guarantees the ordering
     */
    LinkedHashSet<String> getSequenceNames();

    /**
     * do we have an entry for the target chromosome?
     * @param chr the chromosome (or contig) name
     * @return true if we have an entry; false otherwise
     */
    public boolean containsChromosome(final String chr);

    /**
     * read in the index
     * @param stream an input stream to read from
     * @throws IOException if we have problems reading the index from the stream
     */
    public void read(LittleEndianInputStream stream)  throws IOException;

    /**
     * all indexes are writable to disk
     * @param stream the stream to write the index to
     * @throws IOException if the index is unable to write to the specified location
     */
    public void write(LittleEndianOutputStream stream) throws IOException;

    /**
     * this method allows properties to added to the index; warning: if you don't write out the index
     * to disk you'll lose these changes.
     * @param key the key
     * @param value the value, stored as a string, though it may represent an different underlying type
     */
    public void addProperty(String key, String value);

    /**
     * To be called after the index has been created and is ready to be used.  Filling in final metadata or
     * otherwise optimizes the index given that no more records will be added
     */
    public void finalizeIndex();

/**
     * this method allows properties to added to the index
     * @return get the list of properties for this index
     */
    public Map<String,String> getProperties();

    /**
     * Returns true if this and obj are 'effectively' equivalent indices.  Ignores the
     * time stamp on the file, as this may not be the same for even identical indices
     * @param obj
     * @return
     */
    public boolean equalsIgnoreProperties(Object obj);
}
