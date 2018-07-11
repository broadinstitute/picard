package picard.sam.SamErrorMetric;

/** An error metric for the errors invovling bases in the overlapping region of a read-pair.
 * The resulting metric includes error rate information which can be assigned to the reading
 * of the molecular insert {@link #DISAGREES_WITH_REF_AND_MATE_ONLY_Q}, error rate which can be
 * assigned to events that occured to to the molecular insert before it was loaded onto the
 * flowcell/sequencer {@link #DISAGREES_WITH_REFERENCE_ONLY_Q}, and an error rate which
 * cannot be explained nicely {@link #THREE_WAYS_DISAGREEMENT_ONLY_Q}.
 *
 */
public class OverlappingErrorMetric extends ErrorMetric {

    /** The number of bases for which an overlapping base from the mate read was found*/
    @MergeByAdding
    public long NUM_BASES_WITH_OVERLAPPING_READS = 0;

    /** The number of bases that disagree with the reference, but agree with their mate */
    @MergeByAdding
    public long NUM_DISAGREES_WITH_REFERENCE_ONLY = 0L;

    /** The (phred) rate of bases that disagree with the reference, but agree with their mate */
    @NoMergingIsDerived
    public int DISAGREES_WITH_REFERENCE_ONLY_Q = 0;

    /** The number of bases that disagree with both the reference and their mate (which agree with each other) */
    @MergeByAdding
    public long NUM_DISAGREES_WITH_REF_AND_MATE = 0L;

    /** The (phred) rate of bases that disagree with both the reference and their mate (which agree with each other)*/
    @NoMergingIsDerived
    public int DISAGREES_WITH_REF_AND_MATE_ONLY_Q = 0;

    /** The number of bases that disagree with both the reference and their mate (which also disagree) */
    @MergeByAdding
    public long NUM_THREE_WAYS_DISAGREEMENT = 0L;

    /** The (phred) rate of bases that disagree with both the reference and their mate (which also disagree) */
    @NoMergingIsDerived
    public int THREE_WAYS_DISAGREEMENT_ONLY_Q = 0;

    @Override
    public void calculateDerivedFields() {
        this.DISAGREES_WITH_REFERENCE_ONLY_Q = computeQScore(NUM_DISAGREES_WITH_REFERENCE_ONLY, NUM_BASES_WITH_OVERLAPPING_READS);
        this.DISAGREES_WITH_REF_AND_MATE_ONLY_Q = computeQScore(NUM_DISAGREES_WITH_REF_AND_MATE, NUM_BASES_WITH_OVERLAPPING_READS);
        this.THREE_WAYS_DISAGREEMENT_ONLY_Q = computeQScore(NUM_THREE_WAYS_DISAGREEMENT, NUM_BASES_WITH_OVERLAPPING_READS);

    }

    public OverlappingErrorMetric(final String covariate,
                                  final long nTotalBases,
                                  final long nTotalBasesWithOverlappingReads,
                                  final long nDisagreeWithRefAndMate,
                                  final long nDisagreeWithReferenceOnly,
                                  final long nThreeWaysDisagreement) {
        super(covariate, nTotalBases);

        this.NUM_BASES_WITH_OVERLAPPING_READS = nTotalBasesWithOverlappingReads;
        this.NUM_DISAGREES_WITH_REFERENCE_ONLY = nDisagreeWithReferenceOnly;
        this.NUM_DISAGREES_WITH_REF_AND_MATE = nDisagreeWithRefAndMate;
        this.NUM_THREE_WAYS_DISAGREEMENT = nThreeWaysDisagreement;
    }

    // needed for reading in a metric from a file
    public OverlappingErrorMetric() {}

}
