package picard.util.IntervalList;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Interval;

import java.util.List;

/**
 * Created by farjoun on 6/14/18.
 */
public class IntervalListScattererWithSubdivision extends IntervalListScattererByBaseCount {

    @Override
    public List<Interval> takeSome(final Interval interval, final long idealSplitWeight, final long currentSize, final double projectSizeOfRemaining) {
        final long amount = idealSplitWeight - currentSize;

        if (amount >= interval.length()) {
            return CollectionUtil.makeList(interval, null);
        }

        if (amount == 0) {
            return CollectionUtil.makeList(null, interval);
        }

        final Interval left = new Interval(
                interval.getContig(),
                interval.getStart(),
                interval.getStart() + (int) amount - 1,
                interval.isNegativeStrand(),
                interval.getName()
        );
        final Interval right = new Interval(
                interval.getContig(),
                interval.getStart() + (int) amount,
                interval.getEnd(),
                interval.isNegativeStrand(),
                interval.getName()
        );
        return CollectionUtil.makeList(left, right);
    }
}
