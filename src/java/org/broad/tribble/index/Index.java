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

import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

/**
 * Interface for all index implementations.
 * An index file is used for efficient lookup of features from a feature file;
 * and Index represents that index file.
 */
public interface Index {
    /**
     * Query the index.
     * @param chr the chromosome
     * @param start the start position
     * @param end the end position
     * @return a list of blocks that contain the specified interval.  Can never return null
     * @throws IllegalArgumentException of chr isn't part of this index
     */
    List<Block> getBlocks(String chr, int start, int end);

    /**
     * @return true if the index is up to date, false otherwise
     */
    public boolean isCurrentVersion();

    /**
     * @return a list of the sequence names we've seen during indexing, in order
     */
    List<String> getSequenceNames();

    /**
     * @param chr the chromosome (or contig) name
     * @return true if we have an entry; false otherwise
     */
    public boolean containsChromosome(final String chr);

    /**
     * all indexes are writable to disk
     * @param stream the stream to write the index to.  Caller must close after invocation.
     * @throws IOException if the index is unable to write to the specified location
     */
    public void write(LittleEndianOutputStream stream) throws IOException;

    /**
     * Write an appropriately named and located Index file based on the name and location of the featureFile.
     * If featureFile is not a normal file, the index will silently not be written.
     * @param featureFile
     */
    public void writeBasedOnFeatureFile(File featureFile) throws IOException;

    /**
     * @return get the list of properties for this index.  Returns null if no properties.
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
