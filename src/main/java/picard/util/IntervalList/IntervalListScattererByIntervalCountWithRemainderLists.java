package picard.util.IntervalList;

import htsjdk.samtools.util.IntervalList;

/**
 * Scatters {@link IntervalList} by into at least `interval count` shards so that resulting {@link IntervalList}'s have the same number of intervals in them.
 * If the remainder exceeds the count per scatter, then create more lists to divide the remainder evenly
 */
public class IntervalListScattererByIntervalCountWithRemainderLists extends IntervalListScattererByIntervalCount {

    /**
     * For this @{link IntervalListScatterer} we are not limited by the number of lists, so we never reach the limit
     * @param intervalsReturned the number of {@link IntervalList}s already returned
     * @param scatterCount the requested scatter count
     * @return false
     */
    @Override
    public boolean reachedOutputListLimit(final long intervalsReturned, final int scatterCount) {
        return false;
    }
}
