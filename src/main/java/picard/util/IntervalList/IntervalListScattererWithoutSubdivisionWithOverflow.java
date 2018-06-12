package picard.util.IntervalList;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.util.Collections;
import java.util.List;

/**
 * Created by farjoun on 6/14/18.
 */
public class IntervalListScattererWithoutSubdivisionWithOverflow extends IntervalListScattererWithoutSubdivision {

    @Override
    protected boolean shouldIncludeInterval(long idealSplitWeight, double projectedSizeOfRemaining, long projectedSize) {
        return projectedSize <= idealSplitWeight || idealSplitWeight < projectedSizeOfRemaining;
    }

}
