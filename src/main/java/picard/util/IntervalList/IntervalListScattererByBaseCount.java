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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.util.List;

/**
 * a Baseclass for scatterers that scatter by uniqued base count
 */
abstract public class IntervalListScattererByBaseCount implements IntervalListScatterer {

    @Override
    public long intervalWeight(final Interval interval) {
        return interval.length();
    }

    @Override
    public long listWeight(final IntervalList intervalList) {
        return intervalList.getBaseCount();
    }

    @Override
    public int deduceIdealSplitWeight(final IntervalList intervalList, final int nCount) {
        return Math.max(1, (int) Math.floorDiv(intervalList.getUniqueBaseCount(), nCount));
    }
}
