package picard.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.util.Iterator;

/**
 * Takes a VCFFileReader and an IntervalList and provides a single iterator over all variants in all the intervals.
 *
 * //TODO Currently this uses the VCFFileReader.query method - could be useful to make a version of this iterator that uses the .iterator method
 *
 * @author Tim Fennell
 * @author George Grant
 */
public class ByIntervalListVariantContextIterator implements Iterator<VariantContext> {
    private final VCFFileReader reader;
    private final Iterator<Interval> intervals;
    private CloseableIterator<VariantContext> currentCloseableIterator = null;
    private Iterator<VariantContext> currentIterator = null;
    private Interval lastInterval = null;

    /**
     * @param reader the source of variants.
     * @param intervals the intervals to which to restrict variants.
     */
    public ByIntervalListVariantContextIterator(final VCFFileReader reader, final IntervalList intervals) {
        this.reader    = reader;
        this.intervals = intervals.uniqued().iterator();
        this.advance();
    }

    /** Returns true if the variant and interval overlap. */
    private boolean overlapsInterval(final VariantContext ctx, final Interval interval) {
        if (!ctx.getContig().equals(interval.getContig())) return false;
        else if (CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), interval.getStart(), interval.getEnd())) return true;
        else return false;
    }

    /** If the current iterator is null or exhausted, move to the next interval. */
    private void advance() {
        while ((currentIterator == null || !currentIterator.hasNext()) && this.intervals.hasNext()) {
            if (currentIterator != null) currentCloseableIterator.close();
            final Interval interval         = this.intervals.next();
            final Interval previousInterval = this.lastInterval;

            this.currentCloseableIterator = this.reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
            this.currentIterator          = this.currentCloseableIterator.stream().filter ( ctx ->
                null == previousInterval || !overlapsInterval(ctx, previousInterval)
            ).iterator();

            this.lastInterval = interval;
        }
    }

    @Override public boolean hasNext() {
        return (this.currentIterator != null && this.currentIterator.hasNext());
    }

    @Override public VariantContext next() {
        final VariantContext ctx = this.currentIterator.next();
        advance();
        return ctx;
    }

    @Override public void remove() {
        throw new UnsupportedOperationException();
    }
}
