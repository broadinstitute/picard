package picard.filter;

import htsjdk.samtools.SAMRecord;

/** Counting filter that discards reads below a configurable mapping quality threshold. */
public class CountingMapQFilter extends CountingFilter {
    private final int minMapq;

    public CountingMapQFilter(final int minMapq) { this.minMapq = minMapq; }

    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return record.getMappingQuality() < minMapq; }
}
