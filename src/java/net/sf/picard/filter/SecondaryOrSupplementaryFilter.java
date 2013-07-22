package net.sf.picard.filter;

import net.sf.samtools.SAMRecord;

/**
 * Filter out SAMRecords with NotPrimaryAlignment or Supplementary flag set
 * This class should be viewed as a replacement for NotPrimarySkippingIterator,
 * in that we did not want to change the functionality of NPSI to no longer match its name
 * $Id$
 */
public class SecondaryOrSupplementaryFilter  implements SamRecordFilter {
    /**
     * @param record the SAMRecord to evaluate
     * @return true if the SAMRecord matches the filter, otherwise false
     */
    public boolean filterOut(final SAMRecord record) {
        return record.isSecondaryOrSupplementary();
    }

    /**
     * Determines whether a pair of SAMRecord matches this filter
     *
     * @param first  the first SAMRecord to evaluate
     * @param second the second SAMRecord to evaluate
     *
     * @return true if the SAMRecords matches the filter, otherwise false
     */
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        // if either fails, exclude them both
        return first.isSecondaryOrSupplementary() || second.isSecondaryOrSupplementary();
    }
}
