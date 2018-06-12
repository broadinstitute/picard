package picard.util.IntervalList;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.util.Collections;
import java.util.List;

/**
 * Scatters intervalLists by Interval so that resulting interval lists have same number of intervals in them.
 * Final interval can have up to the scatter number extra intervals.
 */
public class IntervalListScattererByInterval implements IntervalListScatterer {

    @Override
    public long intervalWeight(final Interval interval) {
        return 1;
    }

    @Override
    public long listWeight(final IntervalList intervalList) {
        return intervalList.size();
    }

    @Override
    public List<Interval> takeSome(final Interval interval, final long idealSplitWeight, final long currentSize, final double projectSizeOfRemaining) {
        final long amount = idealSplitWeight - currentSize;
        if (amount > 0) {
            return CollectionUtil.makeList(interval, null);
        } else {
            return CollectionUtil.makeList(null, interval);
        }
    }

    @Override
    public int deduceIdealSplitWeight(final IntervalList intervalList, final int nCount) {
        return (int) Math.max(1, Math.floorDiv(listWeight(intervalList), nCount));
    }
}
