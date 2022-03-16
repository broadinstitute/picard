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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * An interface for a class that scatters IntervalLists.
 */
public interface IntervalListScatterer {

    /**
     * Scatter an {@link IntervalList} into several IntervalLists. The default implementation
     * makes use of the other interfaced methods, and aims to provide a universal way to
     * scatter an IntervalList.
     *
     * @param inputList    IntervalList to be scattered
     * @param scatterCount ideal number of scatters generated.
     * @return Scattered {@link List} of IntervalLists,
     */
    default List<IntervalList> scatter(final IntervalList inputList, final int scatterCount) {
        final Iterator<IntervalList> iterator = new IntervalListScatter(this, inputList, scatterCount).iterator();
        final ArrayList<IntervalList> intervalLists = new ArrayList<>();
        iterator.forEachRemaining(intervalLists::add);
        return intervalLists;
    }

    /**
     * A function that will be called on an IntervalList prior to splitting it into sub-lists, and is a point where
     * implementations can chose to impose some conditions on the lists, for example, merging overlapping/abutting intervals,
     * removing duplicates, etc.
     * @param inputList the original {@link IntervalList}
     * @return the  IntervalList that will be split up by the scatterer.
     */
    default IntervalList preprocessIntervalList(final IntervalList inputList) {
        return inputList.sorted();
    }

    /**
     * A method that defines the "weight" of an interval list for the purpose of scattering. The class will attempt to create
     * sublists that all have similar weights.
     */
    long intervalWeight(final Interval interval);

    /**
     *
     * A method that defines the "weight" of an interval for the purpose of scattering. The class will attempt to create
     * sublists that all have similar weights. This method need to estimate the change in any sublists weight due to the possible
     * of the provided interval.
     */
    long listWeight(final IntervalList intervalList);

    /**
     * Figure out how much of the input interval to put into current list and how much to leave for the next interval list.
     *
     * @param interval
     * @return a list of two (possibly null) elements. The first element should be added to the current interval list, the second
     * should be offered to the next interval list.
     */
    List<Interval> takeSome(final Interval interval, final long idealSplitWeight, final long currentSize, final double projectSizeOfRemaining);

    /**
     * A method that determines the ideal target "weight" of the output IntervalList.
     * @param intervalList the {@link IntervalList} that is about to get split
     * @param nCount the scatter count into which to split intervalList
     * @return The ideal "weight" of the output {@link IntervalList}'s
     */
    int deduceIdealSplitWeight(final IntervalList intervalList, final int nCount);
}

