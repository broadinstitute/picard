package picard.sam;

import htsjdk.samtools.metrics.MetricBase;
import picard.PicardException;
import picard.sam.util.SamComparison;

/**
 * Metric for results of SamComparison.  Used to store results in CompareSAMs.
 */
public class SamComparisonMetric extends MetricBase {

    /**
     * Left file used in comparison
     */
    public String leftFile;

    /**
     * Right file used in comparison
     */
    public String rightFile;

    /**
     * The number of primary records for which mappings match in files.  If running in strict alignment mode, this counts only records
     * which are mapped in both files with the same alignment start positions and strands.  If running with LENIENT_LOW_MQ_ALIGNMENT=true, then
     * all records which are mapped in both files with mapping quality at most equal to LOW_MQ_THRESHOLD are counted as matching. If
     * running with LENIENT_255_MQ_ALIGNMENT=true, all records which are mapped in both files with mapping quality 255 are counted as matches.
     */
    public int mappingsMatch;

    /**
     * The number of primary records which are mapped in both files but do not meet criteria to be counted in mappingsMatch.
     */
    public int mappingsDiffer;

    /**
     * The number of primary records which are not mapped in either file.
     */
    public int unmappedBoth;

    /**
     * The number of primary records which are mapped in right file and found but not mapped in left file
     */
    public int unmappedLeft;

    /**
     * The number of primary records which are mapped in left file and found but not mapped in right file
     */
    public int unmappedRight;

    /**
     * The number of primary records which are found in right file but not found in left file
     */
    public int missingLeft;

    /**
     * The number of primary records which are found in left file but not found in right file
     */
    public int missingRight;

    /**
     * The number of primary records for which duplicate markings are different.  If running in strict duplicate
     * marking mode, any primary alignment which is marked as a duplicate in one file but not in the other will
     * be counted.  If running with LENIENT_DUP=true, we allow for swaps between duplicate and non-duplicate fragments
     * in the same duplicate set to reduce the number of fragments with different duplicate marking.  Note that this
     * metric is counted on a per read basis, so a paired end fragment which differs in duplicate marking between the two
     * files will increment this metric by 2.
     */
    public int duplicateMarkingsDiffer;

    /**
     * Whether or not to consider the two input files equal.  The two input files are considered equal iff
     * mappingsDiffer == unmappedLeft == unmappedRight == missingLeft == missingRight == duplicateMarkingsDiffer ==0 &&
     * the headers have been compared to be equal.  Note that the header comparison result can be dependent on whether
     * the tool is run with LENIENT_HEADER true or false.
     */
    public boolean areEqual;

    public boolean allVisitedAlignmentsEqual() {
        return !(missingLeft > 0 || missingRight > 0 || mappingsDiffer > 0 || unmappedLeft > 0 || unmappedRight > 0);
    }

    public void updateMetric(final SamComparison.AlignmentComparison comp) {
        switch (comp) {
            case UNMAPPED_BOTH:
                ++unmappedBoth;
                break;
            case UNMAPPED_LEFT:
                ++unmappedLeft;
                break;
            case UNMAPPED_RIGHT:
                ++unmappedRight;
                break;
            case MAPPINGS_DIFFER:
                ++mappingsDiffer;
                break;
            case MAPPINGS_MATCH:
                ++mappingsMatch;
                break;
            default:
                // unreachable
                throw new PicardException("Unhandled comparison type: " + comp);
        }
    }
}
