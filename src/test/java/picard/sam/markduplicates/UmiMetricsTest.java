/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Tests of {@link UmiMetrics}: its constructors and the fields {@link UmiMetrics#calculateDerivedFields()} computes
 * from the observed and inferred UMIs, the UMIs with Ns, and the base errors.
 */
public class UmiMetricsTest {
    private static final double TOLERANCE = 1e-9;

    private static double log4(final double x) {
        return Math.log(x) / Math.log(4);
    }

    @Test
    public void testConstructors() {
        final UmiMetrics empty = new UmiMetrics();
        Assert.assertNull(empty.LIBRARY);
        Assert.assertEquals(empty.MEAN_UMI_LENGTH, 0.0);
        Assert.assertEquals(empty.OBSERVED_UNIQUE_UMIS, 0);
        Assert.assertEquals(empty.DUPLICATE_SETS_WITH_UMI, 0);

        Assert.assertEquals(new UmiMetrics("lib").LIBRARY, "lib");

        final UmiMetrics full = new UmiMetrics("lib", 6.5, 11, 7, 3, 5, 9, 1.25, 1.5, 17.0, 0.125);
        Assert.assertEquals(full.LIBRARY, "lib");
        Assert.assertEquals(full.MEAN_UMI_LENGTH, 6.5);
        Assert.assertEquals(full.OBSERVED_UNIQUE_UMIS, 11);
        Assert.assertEquals(full.INFERRED_UNIQUE_UMIS, 7);
        Assert.assertEquals(full.OBSERVED_BASE_ERRORS, 3);
        Assert.assertEquals(full.DUPLICATE_SETS_IGNORING_UMI, 5);
        Assert.assertEquals(full.DUPLICATE_SETS_WITH_UMI, 9);
        Assert.assertEquals(full.INFERRED_UMI_ENTROPY, 1.25);
        Assert.assertEquals(full.OBSERVED_UMI_ENTROPY, 1.5);
        Assert.assertEquals(full.UMI_BASE_QUALITIES, 17.0);
        Assert.assertEquals(full.PCT_UMI_WITH_N, 0.125);
    }

    @DataProvider(name = "observations")
    public Object[][] observations() {
        return new Object[][]{
                // observed and inferred UMIs, number of UMIs with N, base errors;
                // expected unique observed, unique inferred, observed entropy, inferred entropy, base quality, fraction with N
                {new String[][]{{"AAAA", "AAAA"}, {"AAAA", "AAAA"}, {"AAAA", "AAAA"}, {"AAAT", "AAAA"}}, 2, 1,
                        2, 1, -(0.75 * log4(0.75) + 0.25 * log4(0.25)), 0.0, 12.0, 2.0 / 6.0},
                // Sixteen equally frequent UMIs carry two bases of entropy; 1 error in 32 bases is Q15.
                {new String[][]{
                        {"AA", "AA"}, {"AC", "AC"}, {"AG", "AG"}, {"AT", "AT"}, {"CA", "CA"}, {"CC", "CC"}, {"CG", "CG"}, {"CT", "CT"},
                        {"GA", "GA"}, {"GC", "GC"}, {"GG", "GG"}, {"GT", "GT"}, {"TA", "TA"}, {"TC", "TC"}, {"TG", "TG"}, {"TT", "TT"}},
                        0, 1, 16, 16, 2.0, 2.0, 15.0, 0.0},
                // 3 errors in 4 bases is Q1.
                {new String[][]{{"A", "A"}, {"C", "A"}, {"G", "A"}, {"T", "A"}}, 0, 3, 4, 1, 1.0, 0.0, 1.0, 0.0},
                {new String[][]{{"ACGT", "ACGT"}, {"ACGA", "ACGT"}}, 1, 1, 2, 1, 0.5, 0.0, 9.0, 1.0 / 3.0},
                {new String[][]{}, 3, 0, 0, 0, 0.0, 0.0, null, 1.0},
        };
    }

    @Test(dataProvider = "observations")
    public void testCalculateDerivedFields(final String[][] observedAndInferred, final int umisWithN, final int baseErrors,
                                           final long expectedObservedUnique, final long expectedInferredUnique,
                                           final double expectedObservedEntropy, final double expectedInferredEntropy,
                                           final Double expectedBaseQuality, final double expectedFractionWithN) {
        final UmiMetrics metrics = new UmiMetrics("lib");
        for (final String[] umis : observedAndInferred) {
            metrics.addUmiObservation(umis[0], umis[1]);
        }
        for (int i = 0; i < umisWithN; i++) {
            metrics.addUmiObservationN();
        }
        metrics.OBSERVED_BASE_ERRORS = baseErrors;
        metrics.MEAN_UMI_LENGTH = 4;
        metrics.DUPLICATE_SETS_IGNORING_UMI = 3;
        metrics.DUPLICATE_SETS_WITH_UMI = 5;

        // Calculating twice gives the same values.
        for (int i = 0; i < 2; i++) {
            metrics.calculateDerivedFields();
            Assert.assertEquals(metrics.OBSERVED_UNIQUE_UMIS, expectedObservedUnique);
            Assert.assertEquals(metrics.INFERRED_UNIQUE_UMIS, expectedInferredUnique);
            Assert.assertEquals(metrics.OBSERVED_UMI_ENTROPY, expectedObservedEntropy, TOLERANCE);
            Assert.assertEquals(metrics.INFERRED_UMI_ENTROPY, expectedInferredEntropy, TOLERANCE);
            if (expectedBaseQuality != null) {
                Assert.assertEquals(metrics.UMI_BASE_QUALITIES, expectedBaseQuality, TOLERANCE);
            }
            Assert.assertEquals(metrics.PCT_UMI_WITH_N, expectedFractionWithN, TOLERANCE);

            // Fields set by the caller are left alone.
            Assert.assertEquals(metrics.LIBRARY, "lib");
            Assert.assertEquals(metrics.OBSERVED_BASE_ERRORS, baseErrors);
            Assert.assertEquals(metrics.MEAN_UMI_LENGTH, 4.0);
            Assert.assertEquals(metrics.DUPLICATE_SETS_IGNORING_UMI, 3);
            Assert.assertEquals(metrics.DUPLICATE_SETS_WITH_UMI, 5);
        }
    }
}
