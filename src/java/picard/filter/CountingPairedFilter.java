package picard.filter;

import htsjdk.samtools.SAMRecord;

/** Counting filter that discards reads that are unpaired in sequencing and paired reads who's mates are not mapped. */
public class CountingPairedFilter extends CountingFilter {
    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return !record.getReadPairedFlag() || record.getMateUnmappedFlag(); }
}
