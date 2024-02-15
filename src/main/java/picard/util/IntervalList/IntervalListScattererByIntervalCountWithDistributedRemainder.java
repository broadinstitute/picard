package picard.util.IntervalList;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.util.Arrays;
import java.util.List;

/**
 * Scatters {@link IntervalList} by into `interval count` shards so that resulting {@link IntervalList}'s have
 * approximately same number of intervals in them. The "remainder" intervals are distributed over the last lists.
 */
public class IntervalListScattererByIntervalCountWithDistributedRemainder extends IntervalListScattererByIntervalCount {

    @Override
    public List<Interval> takeSome(final Interval interval, final long idealSplitWeight, final long currentSize, final double projectSizeOfRemaining) {
        if (projectSizeOfRemaining > currentSize) {
            return Arrays.asList(interval, null);
        } else {
            return Arrays.asList(null, interval);
        }
    }
}
