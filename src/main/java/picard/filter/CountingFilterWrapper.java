package picard.filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * A CountingFilter that wraps any SamRecordFilter and provides a count of the reads and bases filtered
 */
public class CountingFilterWrapper extends CountingFilter {
    private final SamRecordFilter wrappedFilter;
    public CountingFilterWrapper(SamRecordFilter wrappedFilter) {
        this.wrappedFilter = wrappedFilter;
    }

    @Override
    public boolean reallyFilterOut(SAMRecord record) {
        return wrappedFilter.filterOut(record);
    }
}
