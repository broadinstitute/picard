package picard.filter;

import htsjdk.samtools.SAMRecord;

/** Counting filter that discards reads that have been marked as duplicates. */
public class CountingDuplicateFilter extends CountingFilter {
    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return record.getDuplicateReadFlag(); }
}
