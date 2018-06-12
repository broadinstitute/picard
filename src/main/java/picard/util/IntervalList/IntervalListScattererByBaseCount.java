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
