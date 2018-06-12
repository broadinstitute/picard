package picard.util.IntervalList;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import picard.util.IntervalListTools;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;

/**
 * An interface for a class that scatters IntervalLists.
 */
public interface IntervalListScatterer {
    /**
     * Scatter an {@link IntervalList} into several {@link IntervalList}s. The default implementation
     * makes use of the other interfaced methods, and aims to provide a universal way to
     * scatter an {@link IntervalList}.
     *
     * @param inputList IntervalList to be scattered
     * @param scatterCount ideal number of scatters generated.
     * @return Scattered {@link List} of {@link IntervalList}s,
     */
    default List<IntervalList> scatter(final IntervalList inputList, final int scatterCount) {
        if (scatterCount < 1) throw new IllegalArgumentException("scatterCount < 1");

        // How much "weight" should go into each sublist
        final long idealSplitWeight = deduceIdealSplitWeight(inputList, scatterCount);

        Log.getInstance(IntervalListScatterer.class).info("idealSplitWeight=" + idealSplitWeight);

        final List<IntervalList> accumulatedIntervalLists = new ArrayList<>();

        // The IntervalList to which interval are currently being added to.
        IntervalList runningIntervalList = new IntervalList(inputList.getHeader());

        // Use a DeQueue since algo will be adding and removing elements from head.
        final ArrayDeque<Interval> intervalQueue = new ArrayDeque<>(inputList.getIntervals());

        long weightRemaining = listWeight(inputList);

        // Continue processing as long as the queue is not empty, and still haven't generated all scattered lists
        while (!intervalQueue.isEmpty() && accumulatedIntervalLists.size() < scatterCount - 1) {
            final Interval interval = intervalQueue.pollFirst();
            final long currentSize = listWeight(runningIntervalList);

            // The mean expected size of the remaining divisions
            // Note 1: While this looks like double counting, it isn't. We subtract here the bases that are in the _current_ running intervalList,
            // and when we create a new intervalList (below) we modify weightRemaining.
            // Note 2: The -1 in the denominator is for "runningIntervalList" that isn't yet counted in accumulatedIntervalLists.size()

            final double projectedSizeOfRemainingDivisions = (weightRemaining - listWeight(runningIntervalList)) / ((double) scatterCount - accumulatedIntervalLists.size() - 1);

            // split current interval into part that will go into current list (first) and
            // other part that will get put back into queue for next list.
            final List<Interval> split = takeSome(interval, idealSplitWeight, currentSize, projectedSizeOfRemainingDivisions);
            assert split.size() == 2;

            // push second element back to queue (if exists).
            if (split.get(1) != null) {
                intervalQueue.addFirst(split.get(1));
            }

            // if first is null, we are done with the current list, so pack it in.
            if (split.get(0) == null) {
                weightRemaining -= listWeight(runningIntervalList);
                //add running list to return value, and create new running list
                accumulatedIntervalLists.add(runningIntervalList);
                runningIntervalList = new IntervalList(inputList.getHeader());
            } else {
                runningIntervalList.add(split.get(0));
            }
        }

        // Flush the remaining intervals into the last list.
        runningIntervalList.addall(intervalQueue);

        // if last list isn't empty, add it to return value.
        if (!runningIntervalList.getIntervals().isEmpty()) {
            accumulatedIntervalLists.add(runningIntervalList);
        }

        return accumulatedIntervalLists;
    }

    long intervalWeight(final Interval interval);

    long listWeight(final IntervalList intervalList);

    /**
     * figure out how much of input interval to put into current list and how how to leave for the next interval list
     *
     * @param interval
     * @return a list of two (possibly null) elements. The first element should be added to the current interval list, the second
     * should be offered to the next interval list.
     */
    List<Interval> takeSome(final Interval interval, final long idealSplitWeight, final long currentSize, final double projectSizeOfRemaining);

    int deduceIdealSplitWeight(final IntervalList intervalList, final int nCount);
}

