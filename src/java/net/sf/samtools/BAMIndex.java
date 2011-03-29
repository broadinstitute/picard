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

/**
 * A basic interface for querying BAM indices.
 *
 * @author mhanna
 * @version 0.1
 */
public interface BAMIndex {

    public static final String BAMIndexSuffix = ".bai";

    /**
     * Gets the compressed chunks which should be searched for the contents of records contained by the span
     * referenceIndex:startPos-endPos, inclusive.  See the BAM spec for more information on how a chunk is
     * represented.
     * 
     * @param referenceIndex The contig.
     * @param startPos Genomic start of query.
     * @param endPos Genomic end of query.
     * @return A file span listing the chunks in the BAM file.
     */
    BAMFileSpan getSpanOverlapping(final int referenceIndex, final int startPos, final int endPos);

    /**
     * Gets the start of the last linear bin in the index.
     * @return The chunk indicating the start of the last bin in the linear index.
     */
    long getStartOfLastLinearBin();

    /**
     * Gets meta data for the given reference including information about number of aligned, unaligned, and noCoordinate records
     * @param reference the reference of interest
     * @return meta data for the reference
     */
    public BAMIndexMetaData getMetaData(int reference);

    /**
     * Close the index and release any associated resources.
     */
    void close();
}
