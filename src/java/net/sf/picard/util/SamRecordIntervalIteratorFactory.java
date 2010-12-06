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
import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

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
            final int stopAfterSequence;
            final int stopAfterPosition;
            if (uniqueIntervals.isEmpty()) {
                stopAfterSequence = -1;
                stopAfterPosition = -1;
            } else {
                final Interval lastInterval = uniqueIntervals.get(uniqueIntervals.size() - 1);
                stopAfterSequence = samReader.getFileHeader().getSequenceIndex(lastInterval.getSequence());
                stopAfterPosition = lastInterval.getEnd();
            }
            return new StopAfterFilteringIterator(samReader.iterator(), intervalFilter, stopAfterSequence, stopAfterPosition);
        } else {
            // Note that SamRecordIntervalIterator may return some records that do not overlap the intervals,
            // because it merges intervals that are close to one another in order to reduce I/O.  Thus
            // the IntervalFilter is necessary.
            return new FilteringIterator(new SamRecordIntervalIterator(samReader, uniqueIntervals), intervalFilter);
        }
    }

    /**
     * Halt iteration after a read is encountered that starts after the given sequence and position.
     * Note that most of this code is copied from FilteringIterator.  It would be nice just to override getNextRecord,
     * but that method is called FilteringIterator ctor, so the stopAfter members can't be initialized before
     * it is called.
     * FilteringIterator ctor could take a boolean "advance" that would tell it whether or not to call getNextRecord
     * in the ctor, so that it could be delayed in the subclass.  If this pattern happens again, we should do that. 
     */
    private class StopAfterFilteringIterator implements CloseableIterator<SAMRecord> {
        private final int stopAfterSequence;
        private final int stopAfterPosition;
        private final Iterator<SAMRecord> iterator;
        private final SamRecordFilter filter;
        private SAMRecord next = null;

        private StopAfterFilteringIterator(Iterator<SAMRecord> iterator, SamRecordFilter filter,
                                           int stopAfterSequence, int stopAfterPosition) {
            this.stopAfterSequence = stopAfterSequence;
            this.stopAfterPosition = stopAfterPosition;
            this.iterator = iterator;
            this.filter = filter;
            next = getNextRecord();
        }


        /**
         * Returns true if the iteration has more elements.
         *
         * @return  true if the iteration has more elements.  Otherwise returns false.
         */
        public boolean hasNext() {
            return next != null;
        }

        /**
         * Returns the next element in the iteration.
         *
         * @return  the next element in the iteration
         * @throws java.util.NoSuchElementException
         */
        public SAMRecord next() {
            if (next == null) {
                throw new NoSuchElementException("Iterator has no more elements.");
            }
            final SAMRecord result = next;
            next = getNextRecord();
            return result;
        }

        /**
         * Required method for Iterator API.
         *
         * @throws UnsupportedOperationException
         */
        public void remove() {
            throw new UnsupportedOperationException("Remove() not supported by FilteringIterator");
        }

        public void close() {
            CloserUtil.close(iterator);
        }

        protected SAMRecord getNextRecord() {
            while (iterator.hasNext()) {
                SAMRecord record = iterator.next();
                if (record.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) return null;
                else if (record.getReferenceIndex() > stopAfterSequence) return null;
                else if (record.getReferenceIndex() == stopAfterSequence && record.getAlignmentStart() > stopAfterPosition) {
                    return null;
                }
                if (!filter.filterOut(record)) {
                    return record;
                }
            }
            return null;
        }
    }
}
