/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


package picard.sam.SamErrorMetric;

import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.util.help.HelpConstants;

/** An error metric for the errors invovling bases in the overlapping region of a read-pair.
 * The resulting metric includes error rate information which can be assigned to the reading
 * of the molecular insert {@link #DISAGREES_WITH_REF_AND_MATE_ONLY_Q}, error rate which can be
 * assigned to events that occured to to the molecular insert before it was loaded onto the
 * flowcell/sequencer {@link #DISAGREES_WITH_REFERENCE_ONLY_Q}, and an error rate which
 * cannot be explained nicely {@link #THREE_WAYS_DISAGREEMENT_ONLY_Q}.
 *
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
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
