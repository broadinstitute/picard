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

package picard.sam;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Tests {@link DuplicationMetricsFactory}: which metrics class each factory method returns, and that every call
 * returns a new, empty instance.
 */
public class DuplicationMetricsFactoryTest {

    private static void assertEmpty(final DuplicationMetrics metrics) {
        Assert.assertNull(metrics.LIBRARY);
        Assert.assertEquals(DuplicationMetricsBehaviorTest.counts(metrics), new long[7]);
        Assert.assertNull(metrics.PERCENT_DUPLICATION);
        Assert.assertNull(metrics.ESTIMATED_LIBRARY_SIZE);
    }

    @Test
    public void testRegularMetrics() {
        final DuplicationMetrics fromFlag = DuplicationMetricsFactory.createMetrics(false);
        Assert.assertEquals(fromFlag.getClass(), DuplicationMetrics.class);
        assertEmpty(fromFlag);

        final DuplicationMetrics fromDefault = DuplicationMetricsFactory.createMetrics();
        Assert.assertEquals(fromDefault.getClass(), DuplicationMetrics.class);
        assertEmpty(fromDefault);

        Assert.assertNotSame(fromFlag, fromDefault);
        Assert.assertNotSame(DuplicationMetricsFactory.createMetrics(), DuplicationMetricsFactory.createMetrics());
        Assert.assertNotSame(DuplicationMetricsFactory.createMetrics(false), DuplicationMetricsFactory.createMetrics(false));
    }

    @Test
    public void testFlowMetrics() {
        final DuplicationMetrics metrics = DuplicationMetricsFactory.createMetrics(true);
        Assert.assertEquals(metrics.getClass(), FlowBasedDuplicationMetrics.class);
        assertEmpty(metrics);

        final FlowBasedDuplicationMetrics flow = (FlowBasedDuplicationMetrics) metrics;
        Assert.assertEquals(flow.UNPAIRED_WITH_TLEN, 0);
        Assert.assertEquals(flow.UNPAIRED_DUPS_WITH_TLEN, 0);
        Assert.assertEquals(flow.UNPAIRED_DUPS_WITHOUT_TLEN, 0);
        Assert.assertNull(flow.UNPAIRED_DUP_RATE_WITH_TLEN);
        Assert.assertNull(flow.UNPAIRED_DUP_RATE_WITHOUT_TLEN);

        Assert.assertNotSame(DuplicationMetricsFactory.createMetrics(true), DuplicationMetricsFactory.createMetrics(true));
    }

    /** Metrics from separate calls do not share state. */
    @Test
    public void testInstancesAreIndependent() {
        final DuplicationMetrics first = DuplicationMetricsFactory.createMetrics(true);
        final DuplicationMetrics second = DuplicationMetricsFactory.createMetrics(true);
        first.addReadToLibraryMetrics(DuplicationMetricsBehaviorTest.record(0, "8M2S"));
        first.LIBRARY = "first";
        Assert.assertEquals(first.UNPAIRED_READS_EXAMINED, 1);
        Assert.assertEquals(((FlowBasedDuplicationMetrics) first).UNPAIRED_WITH_TLEN, 1);
        assertEmpty(second);
        Assert.assertEquals(((FlowBasedDuplicationMetrics) second).UNPAIRED_WITH_TLEN, 0);
    }
}
