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

import net.sf.picard.PicardException;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import java.util.Iterator;

/**
 * @author alecw@broadinstitute.org
 */
public class IntervalUtil {

    /** Return true if the sequence/position lie in the provided interval. */
    public static boolean contains(final Interval interval, final String sequenceName, final long position) {
        return interval.getSequence().equals(sequenceName) && (position >= interval.getStart() && position <= interval.getEnd());
    }

    /** Return true if the sequence/position lie in the provided interval list. */
    public static boolean contains(final IntervalList intervalList, final String sequenceName, final long position) {
        for (final Interval interval : intervalList.getUniqueIntervals()) {
           if (contains(interval, sequenceName, position))
               return true;
        }
        return false;
    }
    
    /**
     * Throws RuntimeException if the given intervals are not locus ordered and non-overlapping
     * @param intervals
     * @param sequenceDictionary used to determine order of sequences
     */
    public static void assertOrderedNonOverlapping(final Iterator<Interval> intervals, final SAMSequenceDictionary sequenceDictionary) {
        if (!intervals.hasNext()) {
            return;
        }
        Interval prevInterval = intervals.next();
        int prevSequenceIndex = sequenceDictionary.getSequenceIndex(prevInterval.getSequence());
        while (intervals.hasNext()) {
            final Interval interval = intervals.next();
            if (prevInterval.intersects(interval)) {
                throw new PicardException("Intervals should not overlap: " + prevInterval + "; " + interval);
            }
            final int thisSequenceIndex = sequenceDictionary.getSequenceIndex(interval.getSequence());
            if (prevSequenceIndex > thisSequenceIndex ||
                (prevSequenceIndex == thisSequenceIndex && prevInterval.compareTo(interval) >= 0)) {
                throw new PicardException("Intervals not in order: " + prevInterval + "; " + interval);
            }
            prevInterval = interval;
            prevSequenceIndex = thisSequenceIndex;
        }
    }
}
