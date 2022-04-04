/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Histogram;
import picard.analysis.MergeableMetricBase;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;
import picard.util.MathUtil;

import java.util.List;

/**
 * Metrics that are calculated during the process of marking duplicates
 * within a stream of SAMRecords.
 */
public class DuplicationMetrics extends MergeableMetricBase {
    /**
     * The library on which the duplicate marking was performed.
     */
    @MergeByAssertEquals
    public String LIBRARY;

    /**
     * The number of mapped reads examined which did not have a mapped mate pair,
     * either because the read is unpaired, or the read is paired to an unmapped mate.
     */
    @MergeByAdding
    public long UNPAIRED_READS_EXAMINED;

    /**
     * The number of mapped read pairs examined. (Primary, non-supplemental)
     */
    @MergeByAdding
    public long READ_PAIRS_EXAMINED;

    /**
     * The number of reads that were either secondary or supplementary
     */
    @MergeByAdding
    public long SECONDARY_OR_SUPPLEMENTARY_RDS;

    /**
     * The total number of unmapped reads examined. (Primary, non-supplemental)
     */
    @MergeByAdding
    public long UNMAPPED_READS;

    /**
     * The number of fragments that were marked as duplicates.
     */
    @MergeByAdding
    public long UNPAIRED_READ_DUPLICATES;

    /**
     * The number of read pairs that were marked as duplicates.
     */
    @MergeByAdding
    public long READ_PAIR_DUPLICATES;

    /**
     * The number of read pairs duplicates that were caused by optical duplication.
     * Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.
     */
    @MergeByAdding
    public long READ_PAIR_OPTICAL_DUPLICATES;

    /**
     * The fraction of mapped sequence that is marked as duplicate.
     */
    @NoMergingIsDerived
    public Double PERCENT_DUPLICATION;

    /**
     * The estimated number of unique molecules in the library based on PE duplication.
     */
    @NoMergingIsDerived
    public Long ESTIMATED_LIBRARY_SIZE;

    /**
     * Fills in the ESTIMATED_LIBRARY_SIZE based on the paired read data examined where
     * possible and the PERCENT_DUPLICATION.
     */
    @Override
    public void calculateDerivedFields() {
        this.ESTIMATED_LIBRARY_SIZE = estimateLibrarySize(this.READ_PAIRS_EXAMINED - this.READ_PAIR_OPTICAL_DUPLICATES,
                this.READ_PAIRS_EXAMINED - this.READ_PAIR_DUPLICATES);

        if (UNPAIRED_READS_EXAMINED + READ_PAIRS_EXAMINED != 0) {
            PERCENT_DUPLICATION = (UNPAIRED_READ_DUPLICATES + READ_PAIR_DUPLICATES * 2) / (double) (UNPAIRED_READS_EXAMINED + READ_PAIRS_EXAMINED * 2);
        } else {
            PERCENT_DUPLICATION = (double) 0;
        }
    }

    /**
     * Fills in the ESTIMATED_LIBRARY_SIZE based on the paired read data examined where
     * possible and the PERCENT_DUPLICATION.
     * <p>
     * Deprecated, use {@link #calculateDerivedFields()} instead.
     */
    @Deprecated
    public void calculateDerivedMetrics() {
        this.calculateDerivedFields();
    }

    /**
     * Estimates the size of a library based on the number of paired end molecules observed
     * and the number of unique pairs observed.
     * <p>
     * Based on the Lander-Waterman equation that states:
     * C/X = 1 - exp( -N/X )
     * where
     * X = number of distinct molecules in library
     * N = number of read pairs
     * C = number of distinct fragments observed in read pairs
     */
    public static Long estimateLibrarySize(final long readPairs, final long uniqueReadPairs) {
        final long readPairDuplicates = readPairs - uniqueReadPairs;

        if (readPairs > 0 && readPairDuplicates > 0) {

            double m = 1.0;
            double M = 100.0;

            if (uniqueReadPairs >= readPairs || f(m * uniqueReadPairs, uniqueReadPairs, readPairs) < 0) {
                throw new IllegalStateException("Invalid values for pairs and unique pairs: "
                        + readPairs + ", " + uniqueReadPairs);
            }

            // find value of M, large enough to act as other side for bisection method
            while (f(M * uniqueReadPairs, uniqueReadPairs, readPairs) > 0) {
                M *= 10.0;
            }

            // use bisection method (no more than 40 times) to find solution
            for (int i = 0; i < 40; i++) {
                double r = (m + M) / 2.0;
                double u = f(r * uniqueReadPairs, uniqueReadPairs, readPairs);
                if (u == 0) {
                    break;
                } else if (u > 0) {
                    m = r;
                } else if (u < 0) {
                    M = r;
                }
            }

            return (long) (uniqueReadPairs * (m + M) / 2.0);
        } else {
            return null;
        }
    }

    /**
     * Method that is used in the computation of estimated library size.
     */
    private static double f(double x, double c, double n) {
        return c / x - 1 + Math.exp(-n / x);
    }

    /**
     * Estimates the ROI (return on investment) that one would see if a library was sequenced to
     * x higher coverage than the observed coverage.
     *
     * @param estimatedLibrarySize the estimated number of molecules in the library
     * @param x                    the multiple of sequencing to be simulated (i.e. how many X sequencing)
     * @param pairs                the number of pairs observed in the actual sequencing
     * @param uniquePairs          the number of unique pairs observed in the actual sequencing
     * @return a number z <= x that estimates if you had pairs*x as your sequencing then you
     * would observe uniquePairs*z unique pairs.
     */
    public static double estimateRoi(long estimatedLibrarySize, double x, long pairs, long uniquePairs) {
        return estimatedLibrarySize * (1 - Math.exp(-(x * pairs) / estimatedLibrarySize)) / uniquePairs;
    }

    /**
     * Calculates a histogram using the estimateRoi method to estimate the effective yield
     * doing x sequencing for x=1..10.
     */
    public Histogram<Double> calculateRoiHistogram() {
        if (ESTIMATED_LIBRARY_SIZE == null) {
            try {
                calculateDerivedFields();
                if (ESTIMATED_LIBRARY_SIZE == null) {
                    return null;
                }
            } catch (IllegalStateException ise) {
                return null;
            }
        }

        long uniquePairs = READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES;
        Histogram<Double> histo = new Histogram<>();

        for (double x = 1; x <= 100; x += 1) {
            histo.increment(x, estimateRoi(ESTIMATED_LIBRARY_SIZE, x, READ_PAIRS_EXAMINED, uniquePairs));
        }
        histo.setValueLabel("CoverageMult");
        return histo;
    }

    // Main method used for debugging the derived metrics
    // Usage = DuplicationMetrics READ_PAIRS READ_PAIR_DUPLICATES
    public static void main(String[] args) {
        DuplicationMetrics m = new DuplicationMetrics();
        m.READ_PAIRS_EXAMINED = Integer.parseInt(args[0]);
        m.READ_PAIR_DUPLICATES = Integer.parseInt(args[1]);
        m.calculateDerivedFields();
        System.out.println("Percent Duplication: " + m.PERCENT_DUPLICATION);
        System.out.println("Est. Library Size  : " + m.ESTIMATED_LIBRARY_SIZE);
        System.out.println();

        System.out.println("X Seq\tX Unique");
        for (Histogram.Bin<Double> bin : m.calculateRoiHistogram().values()) {
            System.out.println(bin.getId() + "\t" + bin.getValue());
        }
    }

    public void addDuplicateReadToMetrics(final SAMRecord rec) {
        // only update duplicate counts for "decider" reads, not tag-a-long reads
        if (!rec.isSecondaryOrSupplementary() && !rec.getReadUnmappedFlag()) {
            // Update the duplication metrics
            if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                ++UNPAIRED_READ_DUPLICATES;

            } else {
                ++READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
            }
        }
    }

    public void addReadToLibraryMetrics(final SAMRecord rec) {

        // First bring the simple metrics up to date
        if (rec.getReadUnmappedFlag()) {
            ++UNMAPPED_READS;
        } else if (rec.isSecondaryOrSupplementary()) {
            ++SECONDARY_OR_SUPPLEMENTARY_RDS;
        } else if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
            ++UNPAIRED_READS_EXAMINED;
        } else {
            ++READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
        }
    }
}
