package picard.util.IntervalList;

import java.util.function.Supplier;

import htsjdk.samtools.util.IntervalList;

/**
 * Created by farjoun on 6/20/18.
 */
public enum IntervalListScatterMode {
    /**
     * A simple scatter approach in which all output intervals have size equal to the total base count of the source list divide by the
     * scatter count (except for possible variance in the final interval list).
     */
    INTERVAL_SUBDIVISION(IntervalListScattererWithSubdivision::new),
    /**
     * A scatter approach that differs from {@link IntervalListScatterMode#INTERVAL_SUBDIVISION} in a few ways.
     * <ol>
     * <li>No interval will be subdivided, and consequently, the requested scatter count is an upper bound of scatter count, not a
     * guarantee as to how many {@link IntervalList}s will be produced (e.g., if scatterCount = 10 but there is only one input interval,
     * only 1 interval list will be emitted).</li>
     * <li>When an interval would otherwise be split, it is instead deferred to the next scatter list.</li>
     * <li>The "target width" of each scatter list may be wider than what is computed for {@link IntervalListScatterMode#INTERVAL_SUBDIVISION}.
     * Specifically, if the widest interval in the source interval list is larger than what would otherwise be the target width, that
     * interval's width is used.<br/><br/>The reasoning for this is that this approach produces more consistently-sized interval lists,
     * which is one of the objectives of scattering.</li>
     * </ol>
     */
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION(IntervalListScattererWithoutSubdivision::new),
    /**
     * A scatter approach that differs from {@link IntervalListScatterMode#BALANCING_WITHOUT_INTERVAL_SUBDIVISION}.
     * <ol>
     * <li>We try to balance the number of unique bases in each interval list by estimating the remaining interval lists sizes.  This is
     * computed from the total number of unique bases and the bases we have consumed.  This means that the interval list with the most
     * number of unique bases is at most the ideal split length larger than the smallest interval list (unique # of bases).</li>
     * </ol>
     */
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW(IntervalListScattererWithoutSubdivisionWithOverflow::new),

    /**
     * fill this!!!
     */
    SCATTER_BY_INTERVAL_COUNT(IntervalListScattererByInterval::new);

    final private Supplier<IntervalListScatterer> scattererSupplier;

    IntervalListScatterMode(final Supplier<IntervalListScatterer> supplier) {
        scattererSupplier = supplier;
    }

    public IntervalListScatterer make() {
        IntervalList intervalList;
        return scattererSupplier.get();
    }
}