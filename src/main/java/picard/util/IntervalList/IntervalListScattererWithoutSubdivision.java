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
import java.util.Comparator;
import java.util.List;

/**
 * A BaseCount Scatterer that avoid breaking-up intervals. This is done by by only adding intervals to the current list if
 * the resulting size no larger than the "ideal" size. In addition, the ideal length will not be small than the largest sub-interval
 * in the input list.
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
        final int widestIntervalLength = intervalList.getIntervals().stream()
                .map(Interval::length)
                .max(Comparator.comparing(Integer::valueOf))
                .orElse(1);

        // There is no purpose to splitting with more granularity than the widest interval
        return Math.max(widestIntervalLength, splitWidth);
    }
}

