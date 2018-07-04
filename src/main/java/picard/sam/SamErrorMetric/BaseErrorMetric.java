package picard.sam.SamErrorMetric;

/**
 * An error metric for the errors in bases.
 */
public class BaseErrorMetric extends ErrorMetric {
    /** The number of bases that disagree with the reference */
    @MergeByAdding
    public long ERROR_BASES;

    /** The (phred) rate of bases that disagree with the reference */
    @NoMergingIsDerived
    public int Q_SCORE;

    @Override
    public void calculateDerivedFields() {
        this.Q_SCORE = computeQScore(ERROR_BASES);
    }

    public BaseErrorMetric(final String covariate, final long totalBases, final long errorBases) {
        super(covariate, totalBases);
        this.ERROR_BASES = errorBases;
    }

    // needed for reading in a metric from a file
    public BaseErrorMetric() {
    }
}
