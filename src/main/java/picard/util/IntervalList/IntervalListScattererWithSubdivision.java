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

import java.util.List;

/**
 * Created by farjoun on 6/14/18.
 */
public class IntervalListScattererWithSubdivision extends IntervalListScattererByBaseCount {

    @Override
    public List<Interval> takeSome(final Interval interval, final long idealSplitWeight, final long currentSize, final double projectSizeOfRemaining) {
        final long amount = idealSplitWeight - currentSize;

        if (amount >= interval.length()) {
            return CollectionUtil.makeList(interval, null);
        }

        if (amount == 0) {
            return CollectionUtil.makeList(null, interval);
        }

        final Interval left = new Interval(
                interval.getContig(),
                interval.getStart(),
                interval.getStart() + (int) amount - 1,
                interval.isNegativeStrand(),
                interval.getName()
        );
        final Interval right = new Interval(
                interval.getContig(),
                interval.getStart() + (int) amount,
                interval.getEnd(),
                interval.isNegativeStrand(),
                interval.getName()
        );
        return CollectionUtil.makeList(left, right);
    }
}
