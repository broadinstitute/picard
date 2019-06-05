package picard.sam.SamErrorMetric;

/**
 * Metric to be used for InDel errors
 */
public class IndelErrorMetric extends BaseErrorMetric {

    /**
     * The number of insertions. Note: This is not the number of bases that have been inserted.
     */
    @MergeByAdding
    public long NUM_INSERTS = 0;

    /** The (phred) rate of insertions. TODO mgatzen Does this make sense? */
    @NoMergingIsDerived
    public int INSERTS_Q = 0;

    /**
     * The number of deletions. Note: This is not the number of bases that have been deleted.
     */
    @MergeByAdding
    public long NUM_DELETIONS = 0L;

    /** The (phred) rate of deletions. TODO mgatzen Does this make sense? */
    @NoMergingIsDerived
    public int DELETIONS_Q = 0;

    @Override
    public void calculateDerivedFields() {
        this.INSERTS_Q = computeQScore(NUM_INSERTS, TOTAL_BASES);
        this.DELETIONS_Q = computeQScore(NUM_DELETIONS, TOTAL_BASES);

    }

    public IndelErrorMetric(final String covariate,
                                  final long nTotalBases,
                                  final long nInserts,
                                  final long nDeletions) {
        super(covariate, nTotalBases, nInserts + nDeletions);

        this.NUM_INSERTS = nInserts;
        this.NUM_DELETIONS = nDeletions;
    }

    // needed for reading in a metric from a file
    public IndelErrorMetric() {}
}
