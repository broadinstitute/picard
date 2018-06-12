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

