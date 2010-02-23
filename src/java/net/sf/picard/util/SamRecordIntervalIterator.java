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

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Iterate over SAMRecords in an indexed BAM file that overlap the given list of intervals.
 * Note that it is not guaranteed that only reads that overlap one of the given intervals will be
 * returned.  Rather, the interval list is a hint that is used to optimize access to the BAM file.
 * An IntervalFilter can be used to ensure that only SAMRecords overlapping a set of intervals are iterated over.
 * c.f. SamRecordIntervalIteratorFactory.
 *
 * Note that if there are too many intervals or they are too close together, using this class may be slower
 * than iterating through the entire BAM file.
 *
 * @author alecw@broadinstitute.org
 */
class SamRecordIntervalIterator implements CloseableIterator<SAMRecord> {
    private static final int DEFAULT_MAX_READ_LENGTH_GUESS = 16384;

    private final SAMFileReader samReader;
    private final Iterator<Interval> mergedIntervalsIterator;
    /**
     * null implies that there are no more records to return.  Otherwise samIterator.next()
     * will return a record that is appropriate to return to the caller.
     */
    private PeekableIterator<SAMRecord> samIterator = null;

    /**
     * The start of the most recent read returned is remembered.  When advancing to the next interval,
     * it is possible that queryOverlapping() could return a read that was previous returned.  Therefore, when
     * getting a new iterator from queryOverlapping(), it must be advanced so that the next SAMRecord starts
     * after this locus.
     */
    private int lastSequenceIndex = -1;
    private int lastPosition = -1;


    /**
     * Create an iterator that returns reads overlapping one (or more) of the intervals in uniqueIntervals.
     * Note that some of the reads may not be overlap any of the intervals.  The only guarantee is that all
     * reads that overlap will be returned.
     *
     * @param samReader supportsQuery() must return true
     * @param uniqueIntervals must be locus-ordered and non-overlapping.
     */
    public SamRecordIntervalIterator(final SAMFileReader samReader, final List<Interval> uniqueIntervals) {
        this(samReader, uniqueIntervals, DEFAULT_MAX_READ_LENGTH_GUESS);
    }

    /**
     * Create an iterator that returns reads overlapping one (or more) of the intervals in uniqueIntervals.
     * Note that some of the reads may not be overlap any of the intervals.  The only guarantee is that all
     * reads that overlap will be returned.
     *
     * @param samReader supportsQuery() must return true
     * @param uniqueIntervals must be locus-ordered and non-overlapping.
     * @param maxReadLengthGuess Guess for the max read length in the SAM file.  intervals closer together
     * than this are merged when querying in order to avoid reading the same SAMRecord more than once.
     */
    public SamRecordIntervalIterator(final SAMFileReader samReader, final List<Interval> uniqueIntervals, final int maxReadLengthGuess) {
        IntervalUtil.assertOrderedNonOverlapping(uniqueIntervals.iterator(), samReader.getFileHeader().getSequenceDictionary());
        if (!samReader.hasIndex()) {
            throw new IllegalArgumentException("SAMFileReader does not support query");
        }
        this.samReader = samReader;
        this.mergedIntervalsIterator = mergeCloseIntervals(uniqueIntervals, maxReadLengthGuess).iterator();
        advanceInterval();
    }

    private List<Interval> mergeCloseIntervals(final List<Interval> uniqueIntervals, final int maxReadLengthGuess) {
        final List<Interval> ret = new ArrayList<Interval>();
        if (uniqueIntervals.isEmpty()) {
            return ret;
        }
        Interval accumulatingInterval = uniqueIntervals.get(0);
        for (int i = 1; i < uniqueIntervals.size(); ++i) {
            final Interval thisInterval = uniqueIntervals.get(i);
            if (!accumulatingInterval.getSequence().equals(thisInterval.getSequence()) ||
                    thisInterval.getStart() - accumulatingInterval.getEnd() > maxReadLengthGuess) {
                ret.add(accumulatingInterval);
                accumulatingInterval = thisInterval;
            } else {
               accumulatingInterval = new Interval(accumulatingInterval.getSequence(),
                                                   accumulatingInterval.getStart(), thisInterval.getEnd());
            }
        }
        ret.add(accumulatingInterval);
        return ret;
    }


    /**
     * Called when iterator for the current interval has been exhausted.  Get an iterator for the next interval
     * for which there are SAMRecords, and advance it past any already seen.
     */
    private void advanceInterval() {
        if (samIterator != null) {
            samIterator.close();
        }
        samIterator = null;
        while (mergedIntervalsIterator.hasNext()) {
            final Interval nextInterval = mergedIntervalsIterator.next();
            samIterator = new PeekableIterator<SAMRecord>(samReader.queryOverlapping(nextInterval.getSequence(),
                    nextInterval.getStart(), nextInterval.getEnd()));
            // Skip over any SAMRecords already seen.
            advanceSamIterator();
            if (samIterator.hasNext()) {
                // This iterator has some SAMRecords to return.
                break;
            } else {
                // Nothing valid for this interval.  Try the next interval.
                samIterator.close();
                samIterator = null;
            }
        }
    }

    /**
     * Advance the current samIterator past any SAMRecords previously returned.
     */
    private void advanceSamIterator() {
        for (; samIterator.hasNext(); samIterator.next()) {
            final SAMRecord rec = samIterator.peek();
            if (rec.getReferenceIndex() > lastSequenceIndex ||
                    rec.getAlignmentStart() > lastPosition) {
                break;
            }
        }
    }

    public void close() {
        samIterator.close();
        samIterator = null;
    }

    public boolean hasNext() {
        return samIterator != null && samIterator.hasNext();
    }

    public SAMRecord next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        final SAMRecord rec = samIterator.next();
        lastSequenceIndex = rec.getReferenceIndex();
        lastPosition = rec.getAlignmentStart();
        if (!samIterator.hasNext()) {
            advanceInterval();
        }
        return rec;
    }

    public void remove() {
        throw new UnsupportedOperationException("Not supported: remove");
    }
}
