/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

package picard.sam;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Metrics that are calculated during the process of marking duplicates
 * within a stream of SAMRecords using the UmiAwareDuplicateSetIterator.
 */
public class UmiMetrics extends MetricBase {
    // Number of bases in each UMI
    public int UMI_LENGTH;

    // Number of different UMI sequences observed
    public long OBSERVED_UNIQUE_UMIS = 0;

    // Number of different inferred UMI sequences derived
    public long INFERRED_UNIQUE_UMIS = 0;

    // Number of errors inferred by comparing the observed and inferred UMIs
    public long OBSERVED_BASE_ERRORS = 0;

    // Number of duplicate sets found before taking UMIs into account
    public long DUPLICATE_SETS_WITHOUT_UMI = 0;

    // Number of duplicate sets found after taking UMIs into account
    public long DUPLICATE_SETS_WITH_UMI = 0;

    // This is a measure of entropy (in base 4) to indicate the amount of
    // information (given as effective bases) provided by each observed UMI
    public double EFFECTIVE_LENGTH_OF_OBSERVED_UMIS = 0;

    public double EFFECTIVE_LENGTH_OF_INFERRED_UMIS = 0;

    // Estimation of Phred scaled quality scores for UMIs
    public double ESTIMATED_BASE_QUALITY_OF_UMIS;

    // MLE estimation of reads that will be falsely labeled as being part of a duplicate set due to UMI collisions
    public double EXPECTED_READS_WITH_UMI_COLLISION;

    // Phred scale of MLE estimate of collision rate
    public double UMI_COLLISION_Q;

    public void estimateBaseQualities(final int observedUmiBases) {
        ESTIMATED_BASE_QUALITY_OF_UMIS = -10.0*Math.log10((double) OBSERVED_BASE_ERRORS / (double) observedUmiBases);
    }

    public UmiMetrics() {}

    public UmiMetrics(final int length, final int observedUniqueUmis, final int inferredUniqueUmis,
                      final int observedBaseErrors, final int duplicateSetsWithoutUmi,
                      final int duplicateSetsWithUmi, final double effectiveLengthOfInferredUmis,
                      final double effectiveLengthOfObservedUmis, final double estimatedBaseQualityOfUmis,
                      final double expectedUmiCollisions, final double umiCollisionQ) {
        UMI_LENGTH = length;
        OBSERVED_UNIQUE_UMIS = observedUniqueUmis;
        INFERRED_UNIQUE_UMIS = inferredUniqueUmis;
        OBSERVED_BASE_ERRORS = observedBaseErrors;
        DUPLICATE_SETS_WITHOUT_UMI = duplicateSetsWithoutUmi;
        DUPLICATE_SETS_WITH_UMI = duplicateSetsWithUmi;
        EFFECTIVE_LENGTH_OF_INFERRED_UMIS = effectiveLengthOfInferredUmis;
        EFFECTIVE_LENGTH_OF_OBSERVED_UMIS = effectiveLengthOfObservedUmis;
        ESTIMATED_BASE_QUALITY_OF_UMIS = estimatedBaseQualityOfUmis;
        EXPECTED_READS_WITH_UMI_COLLISION = expectedUmiCollisions;
        UMI_COLLISION_Q = umiCollisionQ;
    }
}
