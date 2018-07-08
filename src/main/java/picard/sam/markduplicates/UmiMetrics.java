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

package picard.sam.markduplicates;

import java.util.stream.Collectors;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.QualityUtil;
import picard.PicardException;
import picard.util.MathUtil;

/**
 * Metrics that are calculated during the process of marking duplicates
 * within a stream of SAMRecords using the UmiAwareDuplicateSetIterator.
 */
public class UmiMetrics extends MetricBase {
    private final Histogram<String> observedUmis = new Histogram<>();
    private final Histogram<String> inferredUmis = new Histogram<>();
    private long observedUmiBases = 0;
    private long observedUmiWithNs = 0;
    private long totalObservedUmisWithoutNs = 0;

    /** Library that was used to generate UMI data. */
    public String LIBRARY;

    /** Number of bases in each UMI */
    public double MEAN_UMI_LENGTH = 0.0;

    /** Number of different UMI sequences observed */
    public long OBSERVED_UNIQUE_UMIS = 0;

    /** Number of different inferred UMI sequences derived */
    public long INFERRED_UNIQUE_UMIS = 0;

    /** Number of errors inferred by comparing the observed and inferred UMIs */
    public long OBSERVED_BASE_ERRORS = 0;

    /** Number of duplicate sets found before taking UMIs into account */
    public long DUPLICATE_SETS_IGNORING_UMI = 0;

    /** Number of duplicate sets found after taking UMIs into account */
    public long DUPLICATE_SETS_WITH_UMI = 0;

    /**
     * Entropy (in base 4) of the observed UMI sequences, indicating the
     * effective number of bases in the UMIs.  If this is significantly
     * smaller than UMI_LENGTH, it indicates that the UMIs are not
     * distributed uniformly.
     */
    public double OBSERVED_UMI_ENTROPY = 0;

    /** Entropy (in base 4) of the inferred UMI sequences, indicating the
     * effective number of bases in the inferred UMIs.  If this is significantly
     * smaller than UMI_LENGTH, it indicates that the UMIs are not
     * distributed uniformly.
     */
    public double INFERRED_UMI_ENTROPY = 0;

    /** Estimation of Phred scaled quality scores for UMIs */
    public double UMI_BASE_QUALITIES = 0.0;

    /** The percentage of reads that contain an UMI that contains at least one N */
    public double PCT_UMI_WITH_N = 0.0;

    public UmiMetrics() { }

    public UmiMetrics(final String library) {
        LIBRARY = library;
    }

    public UmiMetrics(final String library, final double length, final int observedUniqueUmis, final int inferredUniqueUmis,
                      final int observedBaseErrors, final int duplicateSetsWithoutUmi,
                      final int duplicateSetsWithUmi, final double effectiveLengthOfInferredUmis,
                      final double effectiveLengthOfObservedUmis, final double estimatedBaseQualityOfUmis,
                      final double percentUmiWithN) {

        LIBRARY = library;
        MEAN_UMI_LENGTH = length;
        OBSERVED_UNIQUE_UMIS = observedUniqueUmis;
        INFERRED_UNIQUE_UMIS = inferredUniqueUmis;
        OBSERVED_BASE_ERRORS = observedBaseErrors;
        DUPLICATE_SETS_IGNORING_UMI = duplicateSetsWithoutUmi;
        DUPLICATE_SETS_WITH_UMI = duplicateSetsWithUmi;
        INFERRED_UMI_ENTROPY = effectiveLengthOfInferredUmis;
        OBSERVED_UMI_ENTROPY = effectiveLengthOfObservedUmis;
        UMI_BASE_QUALITIES = estimatedBaseQualityOfUmis;
        PCT_UMI_WITH_N = percentUmiWithN;
    }

    public void calculateDerivedFields() {
        OBSERVED_UNIQUE_UMIS = observedUmis.size();
        INFERRED_UNIQUE_UMIS = inferredUmis.size();

        // Percent of UMIs that contain at least one N of total number of UMIs observed (observed UMIs is a count of UMIs without N)
        PCT_UMI_WITH_N = (double) observedUmiWithNs/ ((double) observedUmiWithNs + (double) totalObservedUmisWithoutNs);
        OBSERVED_UMI_ENTROPY = effectiveNumberOfBases(observedUmis);
        INFERRED_UMI_ENTROPY = effectiveNumberOfBases(inferredUmis);
        UMI_BASE_QUALITIES = QualityUtil.getPhredScoreFromErrorProbability((double) OBSERVED_BASE_ERRORS / (double) observedUmiBases);
    }

    /**
     * Add an observation of a UMI to the metrics
     * @param observedUmi String containing the observed UMI
     * @param inferredUmi String containing the UMI inferred after error correcting the observed UMI
     */
    public void addUmiObservation(String observedUmi, String inferredUmi) {
        observedUmis.increment(observedUmi);
        inferredUmis.increment(inferredUmi);
        observedUmiBases += observedUmi.length();
        totalObservedUmisWithoutNs++;
    }

    /** Add an observation of a UMI containing at least one N to the metrics */
    public void addUmiObservationN() {
        observedUmiWithNs++;
    }

    private double effectiveNumberOfBases(Histogram<?> observations) {
        double totalObservations = observations.getSumOfValues();

        // Convert to log base 4 so that the entropy is now a measure
        // of the effective number of DNA bases.  If we used log(2.0)
        // our result would be in bits.
        double entropyBaseE = observations.values().stream().collect(Collectors.summingDouble(
                v -> {double p = v.getValue() / totalObservations;
                        return -p * Math.log(p);}));

        return entropyBaseE / MathUtil.LOG_4_BASE_E;
    }
}

