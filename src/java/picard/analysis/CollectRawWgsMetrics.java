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

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments, same implementation as CollectWgsMetrics, with different defaults: lacks baseQ and mappingQ filters and has much higher coverage cap.
 *
 * @author farjoun
 */
@CommandLineProgramProperties(
        usage = "Computes a number of metrics that are useful for evaluating coverage and performance of " +
                "whole genome sequencing experiments. Defaults are different than for CollectWgsMetrics.",
        usageShort = "Writes whole genome sequencing-related metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectRawWgsMetrics extends CollectWgsMetrics{

    @Option(shortName="MQ", doc="Minimum mapping quality for a read to contribute coverage.")
    public int MINIMUM_MAPPING_QUALITY = 0;

    @Option(shortName="Q", doc="Minimum base quality for a base to contribute coverage.")
    public int MINIMUM_BASE_QUALITY = 3;

    @Option(shortName="CAP", doc="Treat bases with coverage exceeding this value as if they had coverage at this value.")
    public int COVERAGE_CAP = 100000;

    // rename the class so that in the metric file it is annotated differently.
    public static class RawWgsMetrics extends WgsMetrics {}

    @Override
    protected WgsMetrics generateWgsMetrics() {
        return new RawWgsMetrics();
    }

}
