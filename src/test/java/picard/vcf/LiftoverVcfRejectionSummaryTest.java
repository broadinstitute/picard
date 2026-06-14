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
package picard.vcf;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Tests for {@link LiftoverVcf#formatRejectionSummary(long, long)}, the helper that produces
 * the per-milestone "Rejected so far" log line emitted during LiftoverVcf's read phase.
 */
public class LiftoverVcfRejectionSummaryTest {

    @DataProvider(name = "rejectionSummaries")
    public Object[][] rejectionSummaries() {
        return new Object[][]{
                // {rejected, processed, expected}
                {0L, 0L, "Rejected so far: 0 of 0 (0.00% reject rate)"},
                {0L, 1_000_000L, "Rejected so far: 0 of 1,000,000 (0.00% reject rate)"},
                {12_345L, 1_000_000L, "Rejected so far: 12,345 of 1,000,000 (1.23% reject rate)"},
                {1L, 100L, "Rejected so far: 1 of 100 (1.00% reject rate)"},
                {1_000_000L, 1_000_000L, "Rejected so far: 1,000,000 of 1,000,000 (100.00% reject rate)"},
        };
    }

    @Test(dataProvider = "rejectionSummaries")
    public void testFormatRejectionSummary(final long rejected, final long processed, final String expected) {
        Assert.assertEquals(LiftoverVcf.formatRejectionSummary(rejected, processed), expected);
    }
}
