/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.AbstractLocusInfo;
import org.testng.annotations.Test;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;
import static org.testng.Assert.assertFalse;
import static picard.analysis.CollectWgsMetricsTestUtils.createIntervalList;

public class AbstractWgsMetricsCollectorTest {

    @Test
    public void testForCollectorWithoutData(){
        long[] templateQualHistogram =  new long[127];
        long[] templateHistogramArray = new long[11];
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        AbstractWgsMetricsCollector collector = new AbstractWgsMetricsCollector(collectWgsMetrics,
                10, createIntervalList()) {
            @Override
            public void addInfo(AbstractLocusInfo info, ReferenceSequence ref, boolean referenceBaseN) {
            }
        };
        assertEquals(templateHistogramArray, collector.highQualityDepthHistogramArray);
        assertEquals(templateQualHistogram, collector.unfilteredBaseQHistogramArray);
        assertEquals(0, collector.basesExcludedByCapping);
        assertEquals(0, collector.basesExcludedByOverlap);
        assertEquals(0, collector.basesExcludedByBaseq);
    }

    @Test (expectedExceptions = IllegalArgumentException.class)
    public void testForExceptionWithNegativeCoverage(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        new AbstractWgsMetricsCollector(collectWgsMetrics, -10, createIntervalList()) {
            @Override
            public void addInfo(AbstractLocusInfo info, ReferenceSequence ref, boolean referenceBaseN) {
            }
        };
    }

    @Test
    public void testForSetCounter(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        AbstractWgsMetricsCollector collector = new AbstractWgsMetricsCollector(collectWgsMetrics,
                10, createIntervalList()) {
            @Override
            public void addInfo(AbstractLocusInfo info, ReferenceSequence ref, boolean referenceBaseN) {
            }
        };
        long counter = 20;
        collector.setCounter(counter);
        assertEquals(20, collector.counter);
    }

    @Test
    public void testForStop(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        collectWgsMetrics.STOP_AFTER = 10;
        AbstractWgsMetricsCollector collector = new AbstractWgsMetricsCollector(collectWgsMetrics,
                10, createIntervalList()) {
            @Override
            public void addInfo(AbstractLocusInfo info, ReferenceSequence ref, boolean referenceBaseN) {
            }
        };
        assertTrue(collector.isTimeToStop(10));
        assertFalse(collector.isTimeToStop(2));
    }

    @Test
    public void testForRefBaseN(){
        byte[] refBasis = {'A', 'C', 'C', 'T', 'A', 'N', 'G', 'T', 'N', 'N'};
        ReferenceSequence ref = new ReferenceSequence("test", 0, refBasis);
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        AbstractWgsMetricsCollector collector = new AbstractWgsMetricsCollector(collectWgsMetrics,
                10, createIntervalList()) {
            @Override
            public void addInfo(AbstractLocusInfo info, ReferenceSequence ref, boolean referenceBaseN) {
            }
        };
        assertTrue(collector.isReferenceBaseN(6, ref));
        assertFalse(collector.isReferenceBaseN(1, ref));
    }
}
