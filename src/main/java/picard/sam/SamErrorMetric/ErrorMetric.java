package picard.sam.SamErrorMetric;

import htsjdk.samtools.util.QualityUtil;
import picard.analysis.MergeableMetricBase;

/**
 * Created by farjoun on 6/26/18.
 */
public class ErrorMetric extends MergeableMetricBase {


    protected static double PRIOR_ERROR = 0.001;
    /**
     * The value of the covariate define the bases included in this metric
     */
    @MergingIsManual
    public String COVARIATE;
    /**
     * The total number of bases included in the calculation of this metric
     */
    @MergeByAdding
    public long TOTAL_BASES;

    public ErrorMetric(final String covariate, final long totalBases) {
        this.TOTAL_BASES = totalBases;
        this.COVARIATE = covariate;
    }

    // required to enable reading metric from a file.
    public ErrorMetric() {
    }

    public static void setPriorError(double priorError) {
        PRIOR_ERROR = priorError;
    }
    /**
     * compute a qscore given the number of errors and the total number of bases.
     * Uses a false count of 1 int the numerator and 1/PRIOR_ERROR in the denominator.
     */
    protected int computeQScore(final long numberOfErrors) {
        return computeQScore(numberOfErrors, TOTAL_BASES);
    }

    /**
     * compute a qscore given the number of errors and the total number of bases.
     * Uses a false count of 1 int the denominator and 1 in the numerator.
     */
    protected int computeQScore(final long numberOfErrors, final long nTotalBases) {
        return QualityUtil.getPhredScoreFromErrorProbability(
                (numberOfErrors + PRIOR_ERROR) / (nTotalBases + 1.0D));
    }
}
