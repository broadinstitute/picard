/*
 * The MIT License
 *
 * Copyright (c) 2016 Nils Homer
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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Tests for DuplicationMetrics.
 */
public class DuplicationMetricsTest {

    private DuplicationMetrics emptyMetrics() {
        final DuplicationMetrics metric       = new DuplicationMetrics();
        metric.LIBRARY                        = "LIBRARY";
        metric.UNPAIRED_READS_EXAMINED        = 0;
        metric.READ_PAIRS_EXAMINED            = 0;
        metric.SECONDARY_OR_SUPPLEMENTARY_RDS = 0;
        metric.UNMAPPED_READS                 = 0;
        metric.UNPAIRED_READ_DUPLICATES       = 0;
        metric.READ_PAIR_DUPLICATES           = 0;
        metric.READ_PAIR_OPTICAL_DUPLICATES   = 0;
        metric.calculateDerivedFields();
        return metric;
    }
    
    private DuplicationMetrics nonEmptyMetrics(final int scale) {
        final DuplicationMetrics metric       = new DuplicationMetrics();
        metric.LIBRARY                        = "LIBRARY";
        metric.UNPAIRED_READS_EXAMINED        = 1000 * scale;
        metric.READ_PAIRS_EXAMINED            = 1000 * scale;
        metric.SECONDARY_OR_SUPPLEMENTARY_RDS = scale;
        metric.UNMAPPED_READS                 = 10 * scale;
        metric.UNPAIRED_READ_DUPLICATES       = 100 * scale;
        metric.READ_PAIR_DUPLICATES           = 110 * scale;
        metric.READ_PAIR_OPTICAL_DUPLICATES   = 10 * scale;
        metric.calculateDerivedFields();
        return metric;
    }

    @Test(dataProvider="testMergeDataProvider")
    public void testMerge(final DuplicationMetrics left, final DuplicationMetrics right, final DuplicationMetrics expected) {
        left.merge(right);
        left.calculateDerivedFields();

        Assert.assertEquals(left.LIBRARY,                        expected.LIBRARY);
        Assert.assertEquals(left.UNPAIRED_READS_EXAMINED,        expected.UNPAIRED_READS_EXAMINED);
        Assert.assertEquals(left.READ_PAIRS_EXAMINED,            expected.READ_PAIRS_EXAMINED);
        Assert.assertEquals(left.SECONDARY_OR_SUPPLEMENTARY_RDS, expected.SECONDARY_OR_SUPPLEMENTARY_RDS);
        Assert.assertEquals(left.UNMAPPED_READS,                 expected.UNMAPPED_READS);
        Assert.assertEquals(left.UNPAIRED_READ_DUPLICATES,       expected.UNPAIRED_READ_DUPLICATES);
        Assert.assertEquals(left.READ_PAIR_DUPLICATES,           expected.READ_PAIR_DUPLICATES);
        Assert.assertEquals(left.READ_PAIR_OPTICAL_DUPLICATES,   expected.READ_PAIR_OPTICAL_DUPLICATES);
        Assert.assertEquals(left.PERCENT_DUPLICATION,            expected.PERCENT_DUPLICATION);
        Assert.assertEquals(left.ESTIMATED_LIBRARY_SIZE,         expected.ESTIMATED_LIBRARY_SIZE);
    }
    
    @DataProvider(name="testMergeDataProvider")
    public Object[][] testMergeDataProvider() {
        return new Object[][] {
                {emptyMetrics(),     emptyMetrics(),     emptyMetrics()},
                {emptyMetrics(),     nonEmptyMetrics(1), nonEmptyMetrics(1)},
                {nonEmptyMetrics(1), emptyMetrics(),     nonEmptyMetrics(1)},
                {nonEmptyMetrics(1), nonEmptyMetrics(1), nonEmptyMetrics(2)},
                {nonEmptyMetrics(1), nonEmptyMetrics(2), nonEmptyMetrics(3)}
        };
    }
}
