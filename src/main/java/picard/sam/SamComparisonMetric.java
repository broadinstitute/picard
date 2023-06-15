package picard.sam;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.sam.util.SamComparison;
import picard.util.help.HelpConstants;

/**
 * Metric for results of SamComparison.  Used to store results in CompareSAMs.
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public class SamComparisonMetric extends MetricBase {

    /**
     * Left file used in comparison
     */
    public String LEFT_FILE;

    /**
     * Right file used in comparison
     */
    public String RIGHT_FILE;

    /**
     * The number of primary records for which mappings match in files.  If running in strict alignment mode, this counts only records
     * which are mapped in both files with the same alignment start positions and strands.  If running with LENIENT_LOW_MQ_ALIGNMENT=true, then
     * all records which are mapped in both files with mapping quality at most equal to LOW_MQ_THRESHOLD are counted as matching. If
     * running with LENIENT_255_MQ_ALIGNMENT=true, all records which are mapped in both files with mapping quality 255 are counted as matches.
     */
    public int MAPPINGS_MATCH;

    /**
     * The number of primary records which are mapped in both files but do not meet criteria to be counted in MAPPINGS_MATCH.
     */
    public int MAPPINGS_DIFFER;

    /**
     * The number of primary records which are not mapped in either file.
     */
    public int UNMAPPED_BOTH;

    /**
     * The number of primary records which are mapped in right file and found but not mapped in left file
     */
    public int UNMAPPED_LEFT;

    /**
     * The number of primary records which are mapped in left file and found but not mapped in right file
     */
    public int UNMAPPED_RIGHT;

    /**
     * The number of primary records which are found in right file but not found in left file
     */
    public int MISSING_LEFT;

    /**
     * The number of primary records which are found in left file but not found in right file
     */
    public int MISSING_RIGHT;

    /**
     * The number of primary records for which duplicate markings are different.  If running in strict duplicate
     * marking mode, any primary alignment which is marked as a duplicate in one file but not in the other will
     * be counted.  If running with LENIENT_DUP=true, we allow for swaps between duplicate and non-duplicate fragments
     * in the same duplicate set to reduce the number of fragments with different duplicate marking.  Note that this
     * metric is counted on a per read basis, so a paired end fragment which differs in duplicate marking between the two
     * files will increment this metric by 2.
     */
    public int DUPLICATE_MARKINGS_DIFFER;

    /**
     * Whether or not to consider the two input files equal.  The two input files are considered equal iff
     * MAPPINGS_DIFFER == UNMAPPED_LEFT == UNMAPPED_RIGHT == MISSING_LEFT == MISSING_RIGHT == DUPLICATE_MARKINGS_DIFFER == 0 &&
     * the headers have been compared to be equal.  Note that the header comparison result can be dependent on whether
     * the tool is run with LENIENT_HEADER true or false.
     */
    public boolean ARE_EQUAL;

    public boolean allVisitedAlignmentsEqual() {
        return !(MISSING_LEFT > 0 || MISSING_RIGHT > 0 || MAPPINGS_DIFFER > 0 || UNMAPPED_LEFT > 0 || UNMAPPED_RIGHT > 0);
    }

    public void updateMetric(final SamComparison.AlignmentComparison comp) {
        switch (comp) {
            case UNMAPPED_BOTH:
                ++UNMAPPED_BOTH;
                break;
            case UNMAPPED_LEFT:
                ++UNMAPPED_LEFT;
                break;
            case UNMAPPED_RIGHT:
                ++UNMAPPED_RIGHT;
                break;
            case MAPPINGS_DIFFER:
                ++MAPPINGS_DIFFER;
                break;
            case MAPPINGS_MATCH:
                ++MAPPINGS_MATCH;
                break;
            default:
                // unreachable
                throw new PicardException("Unhandled comparison type: " + comp);
        }
    }
}
