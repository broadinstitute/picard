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

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.AbstractLocusIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import static org.testng.Assert.assertEquals;
import static picard.analysis.CollectWgsMetricsTestUtils.createIntervalList;
import static picard.analysis.CollectWgsMetricsTestUtils.createReadEndsIterator;
import static picard.analysis.CollectWgsMetricsTestUtils.exampleSamOneRead;
import static picard.analysis.CollectWgsMetricsTestUtils.getReferenceSequenceFileWalker;
import static picard.analysis.CollectWgsMetricsTestUtils.createReferenceSequenceFile;

public class WgsMetricsProcessorImplTest {
    private ProgressLogger progress;

    @BeforeTest
    public void setUp(){
        progress = new ProgressLogger(Log.getInstance(WgsMetricsProcessorImpl.class));
    }

    @Test
    public void testForProcessFile(){
        AbstractLocusIterator iterator = createReadEndsIterator(exampleSamOneRead);
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector =  new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        WgsMetricsProcessorImpl wgsMetricsProcessor =
                new WgsMetricsProcessorImpl(iterator, getReferenceSequenceFileWalker(), collector, progress);
        wgsMetricsProcessor.processFile();
        assertEquals(20, collector.counter);
    }

    @Test
    public void testForFilteredBases(){
        AbstractLocusIterator iterator = createReadEndsIterator(exampleSamOneRead);
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector =  new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        String secondReferenceString = ">ref\nNNNNNNNNNNAATATTCTTC";
        ReferenceSequenceFile referenceSequenceFile = createReferenceSequenceFile(secondReferenceString);
        ReferenceSequenceFileWalker refWalker =  new ReferenceSequenceFileWalker(referenceSequenceFile);
        WgsMetricsProcessorImpl wgsMetricsProcessor =
                new WgsMetricsProcessorImpl(iterator, refWalker, collector, progress);
        wgsMetricsProcessor.processFile();
        assertEquals(10, collector.counter);
    }

    @Test
    public void testForExitAfter(){
        AbstractLocusIterator iterator = createReadEndsIterator(exampleSamOneRead);
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        collectWgsMetrics.STOP_AFTER = 16;
        AbstractWgsMetricsCollector collector =  new FastWgsMetricsCollector(collectWgsMetrics, 100, null);
        WgsMetricsProcessorImpl wgsMetricsProcessor =
                new WgsMetricsProcessorImpl(iterator, getReferenceSequenceFileWalker(), collector, progress);
        wgsMetricsProcessor.processFile();
        assertEquals(15, collector.counter);
    }
}
