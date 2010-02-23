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
package net.sf.picard.util;

import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.IntervalFilter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.util.List;

/**
 * Create an iterator over a SAMFileReader that only returns reads that overlap one of the intervals
 * in an interval list.
 *
 * @author alecw@broadinstitute.org
 */
public class SamRecordIntervalIteratorFactory {

    /**
     * @param samReader
     * @param uniqueIntervals list of intervals of interest, with overlaps merged, in coordinate order
     * @param useIndex if false, do not use a BAM index even if it is present.
     * @return an iterator that will be filtered so that only SAMRecords overlapping the intervals
     * in uniqueIntervals will be returned.  If a BAM index is available, it will be used to improve performance.
     * Note however that if there are many intervals that cover a great deal of the genome, using the BAM
     * index may actually make performance worse.
     */
    public CloseableIterator<SAMRecord> makeSamRecordIntervalIterator(final SAMFileReader samReader,
                                                               final List<Interval> uniqueIntervals,
                                                               final boolean useIndex) {
        final IntervalFilter intervalFilter = new IntervalFilter(uniqueIntervals, samReader.getFileHeader());
        if (!samReader.hasIndex() || !useIndex) {
            return new FilteringIterator(samReader.iterator(), intervalFilter);
        } else {
            // Note that SamRecordIntervalIterator may return some records that do not overlap the intervals,
            // because it merges intervals that are close to one another in order to reduce I/O.  Thus
            // the IntervalFilter is necessary.
            return new FilteringIterator(new SamRecordIntervalIterator(samReader, uniqueIntervals), intervalFilter);
        }
    }
}
