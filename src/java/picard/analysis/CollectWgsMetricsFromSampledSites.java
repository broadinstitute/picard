/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.filter.CountingFilter;

import java.io.File;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments,
 * but only at a set of sampled positions.
 * It is important that the sampled positions be chosen so that they are spread out at least further than a read's length apart;
 * otherwise, you run the risk of double-counting reads in the metrics.
 *
 * @author ebanks
 */
@CommandLineProgramProperties(
        usage = "Computes a number of metrics that are useful for evaluating coverage and performance of " +
                "whole genome sequencing experiments, but only at a set of sampled positions.  " +
                "It is important that the sampled positions be chosen so that they are spread out " +
                "at least further than a read's length apart; otherwise, you run the risk of double-counting " +
                "reads in the metrics.",
        usageShort = "Writes whole genome sequencing-related metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectWgsMetricsFromSampledSites extends CollectWgsMetrics {

    @Option(shortName = "INTERVALS", doc = "An interval list file that contains the locations of the positions to assess.", optional = false)
    public File INTERVALS = null;

    public static void main(final String[] args) {
        new CollectWgsMetricsFromSampledSites().instanceMainWithExit(args);
    }

    @Override
    protected SamLocusIterator getLocusIterator(final SamReader in) {
        IOUtil.assertFileIsReadable(INTERVALS);
        return new SamLocusIterator(in, IntervalList.fromFile(INTERVALS));
    }

    /**
     * By design we want to count just those bases at the positions we care about, not across the entire read.
     * Therefore, we call filter.getFilteredRecords() so that only the bases in the pileup at a given position
     * are included in the calculations (with filter.getFilteredBases() we would be including other bases in
     * the read too).
     */
    @Override
    protected long getBasesExcludedBy(final CountingFilter filter) {
        return filter.getFilteredRecords();
    }

    // rename the class so that in the metric file it is annotated differently.
    public static class SampledWgsMetrics extends WgsMetrics {}

    @Override
    protected WgsMetrics generateWgsMetrics() {
        return new SampledWgsMetrics();
    }
}

