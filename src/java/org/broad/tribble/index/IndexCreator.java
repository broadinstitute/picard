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

import org.broad.tribble.Feature;

import java.io.File;

/**
 *
 * An interface for creating indexes
 *
 * @author jrobinso
 */                                                                           
public interface IndexCreator {
    /**
     * Initialize the index creator with the input file and the bin size.  Be warned, the bin size
     * is HIGHLY dependent on the index implementation; in one implementation 100 may result in excessively
     * large files, and other this may be too small for effective discernment between bins.  It's recommended to
     * use the defaultBinSize() function to get an appropriately sized bin.
     *
     * @param inputFile the input file
     * @param binSize the bin size
     */
    public void initialize(File inputFile, int binSize);

    /**
     * Add a feature to the index
     * @param feature the feature, of which start, end, and contig must be filled in
     * @param filePosition the current file position, at the beginning of the specified feature
     */
    public void addFeature(Feature feature, long filePosition);

    /**
     * Create the index, given the stream of features passed in to this point
     * @param finalFilePosition the final file position, for indexes that have to close out with the final position
     * @return an index object
     */
    public Index finalizeIndex(long finalFilePosition);

    /**
     * The default bin size for this index type; use this unless you're aware of the nuances of the particular index type.
     * @return the default bin size appropriate for this index type
     */
    public int defaultBinSize();

    /**
     * @eturn the bin size of associated with the index with are creating
     * @return the index bin size
     */
    public int getBinSize();
}


