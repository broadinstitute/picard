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

import net.sf.samtools.util.CloseableIterator;

/**
 * A general interface that adds functionality to a CloseableIterator of
 * SAMRecords.  Currently, this interface is implemented by iterators that
 * want to validate as they are iterating that that the records in the
 * underlying SAM/BAM file are in a particular order.
 */
public interface SAMRecordIterator extends CloseableIterator<SAMRecord> {

    /**
     * Establishes that records returned by this iterator are expected to
     * be in the specified sort order.  If this method has been called,
     * then implementers must throw an IllegalStateException from next()
     * when a record is read that violates the sort order.  This method
     * may be called multiple times over the course of an iteration,
     * changing the expected sort, if desired -- from the time it is called,
     * it validates whatever sort is set, or stops validating if it
     * is set to null or SAMFileHeader.SortOrder.unsorted.  If this method
     * is not called, then no validation of the iterated records is done.
     *
     * @param sortOrder The order in which records are expected to be returned
     * @return  This SAMRecordIterator
     */
    public SAMRecordIterator assertSorted(SAMFileHeader.SortOrder sortOrder);
    
}
