package picard.util;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import picard.PicardException;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author mccowan
 */
public class IntervalListScatterer {

    public enum Mode {
        /**
         * A simple scatter approach in which all output intervals have size equal to the total base count of the source list divide by the
         * scatter count (except for possible variance in the final interval list).
         */
        INTERVAL_SUBDIVISION,
        /**
         * A scatter approach that differs from {@link Mode#INTERVAL_SUBDIVISION} in a few ways.
         * <ol>
         * <li>No interval will be subdivided, and consequently, the requested scatter count is an upper bound of scatter count, not a
         * guarantee as to how many {@link IntervalList}s will be produced (e.g., if scatterCount = 10 but there is only one input interval,
         * only 1 interval list will be emitted).</li>
         * <li>When an interval would otherwise be split, it is instead deferred to the next scatter list.</li>
         * <li>The "target width" of each scatter list may be wider than what is computed for {@link Mode#INTERVAL_SUBDIVISION}.
         * Specifically, if the widest interval in the source interval list is larger than what would otherwise be the target width, that
         * interval's width is used.<br/><br/>The reasoning for this is that this approach produces more consistently-sized interval lists,
         * which is one of the objectives of scattering.</li>
         * </ol>
         */
        BALANCING_WITHOUT_INTERVAL_SUBDIVISION,
        /**
         * A scatter approach that differs from {@link Mode#BALANCING_WITHOUT_INTERVAL_SUBDIVISION}.
         * <ol>
         * <li>We try to balance the number of unique bases in each interval list by estimating the remaining interval lists sizes.  This is 
         * computed from the total number of unique bases and the bases we have consumed.  This means that the interval list with the most
         * number of unique bases is at most the ideal split length larger than the smallest interval list (unique # of bases).</li>
         * </ol>         
         */
        BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
    }

    private final Mode mode;

    public IntervalListScatterer(final Mode mode) {this.mode = mode;}

    private int deduceIdealSplitLength(final IntervalList uniquedList, final int scatterCount) {
        final int splitWidth = Math.max((int) Math.floor(uniquedList.getBaseCount() / (1.0 * scatterCount)), 1);
        switch (mode) {
            case INTERVAL_SUBDIVISION:
                return splitWidth;
            case BALANCING_WITHOUT_INTERVAL_SUBDIVISION:
            case BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW:
                final int widestIntervalLength = Collections.max(uniquedList.getIntervals(), new Comparator<Interval>() {
                    @Override
                    public int compare(final Interval o1, final Interval o2) {
                        return Integer.valueOf(o1.length()).compareTo(o2.length());
                    }
                }).length();

                // There is no purpose to splitting more granularly than the widest interval, so do not.
                return Math.max(widestIntervalLength, splitWidth);
            default:
                throw new IllegalStateException();
        }
    }

    public List<IntervalList> scatter(final IntervalList uniquedIntervalList, final int scatterCount) {
        return scatter(uniquedIntervalList, scatterCount, false);
    }
    
    /** Helper for the scatter method */
    private boolean shouldAddToRunningIntervalList(final long idealSplitLength, final long projectedSize, final double projectedSizeOfRemainingDivisions) {
        switch (mode) {
            case BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW:
                return (projectedSize <= idealSplitLength || idealSplitLength < projectedSizeOfRemainingDivisions);
            default:
                return (projectedSize <= idealSplitLength);
        }
    }


    public List<IntervalList> scatter(final IntervalList sourceIntervalList, final int scatterCount, final boolean isUniqued) {
        if (scatterCount < 1) throw new IllegalArgumentException("scatterCount < 1");

        final IntervalList uniquedList = isUniqued ? sourceIntervalList : sourceIntervalList.uniqued();
        final long idealSplitLength = deduceIdealSplitLength(uniquedList, scatterCount);
        System.err.println("idealSplitLength=" + idealSplitLength);

        final List<IntervalList> accumulatedIntervalLists = new ArrayList<IntervalList>();

        IntervalList runningIntervalList = new IntervalList(uniquedList.getHeader());
        final ArrayDeque<Interval> intervalQueue = new ArrayDeque<Interval>(uniquedList.getIntervals());
        
        long numBasesLeft = uniquedList.getBaseCount();

        while (!intervalQueue.isEmpty() && accumulatedIntervalLists.size() < scatterCount - 1) {
            final Interval interval = intervalQueue.pollFirst();
            final long projectedSize = runningIntervalList.getBaseCount() + interval.length();

            // The mean expected size of the remaining divisions
            // NOTE: that this looks like double counting but isn't, we subtract here the bases that are in the _current_ running intervalList,
            // and when we create a new intervalList (below) we modify numBasesLeft.
            // Another Note: the -1 in the denominator is for "runningIntervalList" that isn't yet counted in  accumulatedIntervalLists.size()
            final double projectedSizeOfRemainingDivisions = (numBasesLeft - runningIntervalList.getBaseCount()) / ((double)(scatterCount - accumulatedIntervalLists.size() - 1));

            // should we add this interval to the list of running intervals?
            if (shouldAddToRunningIntervalList(idealSplitLength, projectedSize, projectedSizeOfRemainingDivisions)) {
                runningIntervalList.add(interval);
            }
            else {
                switch (mode) {
                    case INTERVAL_SUBDIVISION:
                        final int amountToConsume = (int) (idealSplitLength - runningIntervalList.getBaseCount());
                        final Interval left = new Interval(
                                interval.getContig(),
                                interval.getStart(),
                                interval.getStart() + amountToConsume - 1,
                                interval.isNegativeStrand(),
                                interval.getName()
                        );
                        final Interval right = new Interval(
                                interval.getContig(),
                                interval.getStart() + amountToConsume,
                                interval.getEnd(),
                                interval.isNegativeStrand(),
                                interval.getName()
                        );
                        runningIntervalList.add(left);

                        // Push back the excess back onto our queue for reconsideration.
                        intervalQueue.addFirst(right);
                        break;

                    case BALANCING_WITHOUT_INTERVAL_SUBDIVISION:
                    case BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW:
                        if (runningIntervalList.getIntervals().isEmpty()) {
                            runningIntervalList.add(interval);
                        } else {
                            // Push this interval into the next scatter; re-inject it into the queue, then advance the scatter.
                            intervalQueue.addFirst(interval);
                            numBasesLeft -= runningIntervalList.getBaseCount();
                            accumulatedIntervalLists.add(runningIntervalList.uniqued());
                            runningIntervalList = new IntervalList(uniquedList.getHeader());
                        }
                        break;
                }
            }

            if (runningIntervalList.getBaseCount() >= idealSplitLength) {
                numBasesLeft -= runningIntervalList.getBaseCount(); // keep track of the number of *unique* bases left
                accumulatedIntervalLists.add(runningIntervalList.uniqued());
                runningIntervalList = new IntervalList(uniquedList.getHeader());
            }
        }

        // Flush the remaining intervals into the last split.
        while (!intervalQueue.isEmpty()) {
            runningIntervalList.add(intervalQueue.pollFirst());
        }
        if (!runningIntervalList.getIntervals().isEmpty()) {
            accumulatedIntervalLists.add(runningIntervalList.uniqued());
        }
        
        long maximumIntervalSize = -1, minimumIntervalSize = Integer.MAX_VALUE;
        for (final IntervalList intervalList : accumulatedIntervalLists) {
            final long baseCount = intervalList.getBaseCount();
            if (baseCount < minimumIntervalSize) minimumIntervalSize = baseCount;
            if (maximumIntervalSize < baseCount) maximumIntervalSize = baseCount;
        }

        return accumulatedIntervalLists;
    }
}
