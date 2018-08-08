/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

package picard.util.IntervalList;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.util.Collections;
import java.util.List;

/**
 * Created by farjoun on 6/14/18.
 */
public class IntervalListScattererWithoutSubdivision extends IntervalListScattererByBaseCount {

    @Override
    public List<Interval> takeSome(final Interval interval, final long idealSplitWeight, final long currentSize, final double projectedSizeOfRemaining) {
        final long projectedSize = currentSize + intervalWeight(interval);
        if (shouldIncludeInterval(idealSplitWeight, projectedSizeOfRemaining, projectedSize)) {
            return CollectionUtil.makeList(interval, null);
        } else {
            return CollectionUtil.makeList(null, interval);
        }
    }

    protected boolean shouldIncludeInterval(long idealSplitWeight, double projectedSizeOfRemaining, long projectedSize) {
        return projectedSize <= idealSplitWeight;
    }

    @Override
    public int deduceIdealSplitWeight(final IntervalList intervalList, final int nCount) {
        final int splitWidth = super.deduceIdealSplitWeight(intervalList, nCount);
        final int widestIntervalLength = Collections.max(intervalList.uniqued().getIntervals(), (o1, o2) -> Integer.valueOf(o1.length()).compareTo(o2.length())).length();

        // There is no purpose to splitting with more granularity than the widest interval, so do not.
        return Math.max(widestIntervalLength, splitWidth);
    }
}

