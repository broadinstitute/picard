package picard.util.IntervalList;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

public class IntervalListScatter implements Iterable<IntervalList> {

    private final IntervalListScatterer scatterer;
    private final IntervalList intervals;
    private final int scatterCount;

    public IntervalListScatter(final IntervalListScatterer scatterer, final IntervalList intervals, final int scatterCount) {
        if (scatterCount < 1) {
            throw new IllegalArgumentException("scatterCount < 1");
        }
        this.scatterer = scatterer;
        this.intervals = scatterer.preprocessIntervalList(intervals);
        this.scatterCount = scatterCount;
    }

    @Override
    public Iterator<IntervalList> iterator() {
        return new ScatterState(scatterer, intervals, scatterCount);
    }

    class ScatterState implements Iterator<IntervalList> {

        private final int scatterCount;
        private final int idealSplitWeight;
        private final IntervalListScatterer scatterer;
        private final SAMFileHeader header;

        // Use a DeQueue since algo will be adding and removing elements from head.
        final ArrayDeque<Interval> intervalQueue;
        long weightRemaining;
        long intervalsReturned = 0;

        public ScatterState(final IntervalListScatterer scatterer, final IntervalList inputIntervals, final int scatterCount) {
            this.scatterCount = scatterCount;
            final IntervalList processedIntervals = scatterer.preprocessIntervalList(inputIntervals);
            this.idealSplitWeight = scatterer.deduceIdealSplitWeight(processedIntervals, scatterCount);
            this.intervalQueue = new ArrayDeque<>(processedIntervals.getIntervals());
            this.scatterer = scatterer;
            this.header = processedIntervals.getHeader();
            this.weightRemaining = scatterer.listWeight(processedIntervals);
        }

        @Override
        public boolean hasNext() {
            return !intervalQueue.isEmpty();
        }

        @Override
        public IntervalList next() {
            intervalsReturned++;
            final IntervalList runningIntervalList = new IntervalList(header);

            while (!intervalQueue.isEmpty() && intervalsReturned < scatterCount ) {
                final Interval interval = intervalQueue.pollFirst();
                final long currentSize = scatterer.listWeight(runningIntervalList);

                // The mean expected size of the remaining divisions
                // Note 1: While this looks like double counting, it isn't. We subtract here the bases that are in the _current_ running intervalList,
                // and when we create a new intervalList (below) we modify weightRemaining.
                // Note 2: The -1 in the denominator is for "runningIntervalList" that isn't yet counted in accumulatedIntervalLists.size()

                final double projectedSizeOfRemainingDivisions = (weightRemaining - scatterer.listWeight(runningIntervalList)) / ((double) scatterCount - intervalsReturned);

                // split current interval into part that will go into current list (first) and
                // other part that will get put back into queue for next list.
                final List<Interval> split = scatterer.takeSome(interval, idealSplitWeight, currentSize, projectedSizeOfRemainingDivisions);
                if (split.size() != 2) {
                    throw new IllegalStateException("takeSome should always return exactly 2 (possibly null) intervals.");
                }

                // push second element back to queue (if exists).
                if (split.get(1) != null) {
                    intervalQueue.addFirst(split.get(1));
                }

                // if first is null, we are done with the current list, so pack it in.
                if (split.get(0) == null) {
                    weightRemaining -= scatterer.listWeight(runningIntervalList);
                    //add running list to return value, and create new running list
                    return runningIntervalList;
                } else {
                    runningIntervalList.add(split.get(0));
                }
            }
            // Flush the remaining intervals into the last list.
            while(!intervalQueue.isEmpty()){
                runningIntervalList.add(intervalQueue.pollFirst());
            }

            if (!runningIntervalList.getIntervals().isEmpty()){
                return runningIntervalList;
            }
            // you asked for a next() when hasNext() was false
            throw new NoSuchElementException();
        }
    }
}
