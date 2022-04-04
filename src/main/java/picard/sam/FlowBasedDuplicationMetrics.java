package picard.sam;

import htsjdk.samtools.SAMRecord;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;
import picard.util.MathUtil;

public class FlowBasedDuplicationMetrics extends DuplicationMetrics {

    /*
     * count of single end reads where the exact fragment length is known (i.e. clipped)
     */
    @MergeByAdding
    public long UNPAIRED_WITH_TLEN;

    /*
     * count of single end duplicates where the exact fragment length is
     * unknown (quality trimmed, not clipped)
     */
    @MergeByAdding
    public long UNPAIRED_DUPS_WITHOUT_TLEN;

    /*
     * count of duplicates where both ends are known
     */
    @MergeByAdding
    public long UNPAIRED_DUPS_WITH_TLEN;

    /**
     * The fraction of duplicated reads out of all reads with exact
     * fragment length unknown
     */
    @NoMergingIsDerived
    public Double UNPAIRED_DUP_RATE_WITHOUT_TLEN;

    /**
     * The fraction of duplicated reads out of all reads with exact fragment
     * length known
     */
    @NoMergingIsDerived
    public Double UNPAIRED_DUP_RATE_WITH_TLEN;


    @Override
    public void calculateDerivedFields() {
        super.calculateDerivedFields();

        UNPAIRED_DUP_RATE_WITHOUT_TLEN = MathUtil.divide(UNPAIRED_DUPS_WITHOUT_TLEN, UNPAIRED_READS_EXAMINED - UNPAIRED_WITH_TLEN);
        UNPAIRED_DUP_RATE_WITH_TLEN = MathUtil.divide(UNPAIRED_DUPS_WITH_TLEN, UNPAIRED_WITH_TLEN);
    }

    public void addDuplicateReadToMetrics(final SAMRecord rec) {
        super.addDuplicateReadToMetrics(rec);

        if (!rec.isSecondaryOrSupplementary() && !rec.getReadUnmappedFlag()) {
            if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                if ( AbstractMarkDuplicatesCommandLineProgram.isSingleEndReadKnownFragment(rec) ) {
                    ++UNPAIRED_DUPS_WITH_TLEN;
                } else {
                    ++UNPAIRED_DUPS_WITHOUT_TLEN;
                }
            }
        }
    }

    public void addReadToLibraryMetrics(final SAMRecord rec) {

        super.addReadToLibraryMetrics(rec);

        if (AbstractMarkDuplicatesCommandLineProgram.isSingleEndReadKnownFragment(rec)) {
            ++UNPAIRED_WITH_TLEN;
        }
    }

}
