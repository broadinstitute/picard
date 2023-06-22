/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

package picard.analysis.directed;

import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.util.help.HelpConstants;

/**
 * TargetMetrics, are metrics to measure how well we hit specific targets (or baits) when using a targeted sequencing process like hybrid selection
 * or Targeted PCR Techniques (TSCA).  TargetMetrics at the moment are the metrics that are shared by both HybridSelection and TargetedPcrMetrics.
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public class TargetMetrics extends TargetMetricsBase {
    /**  The name of the PROBE_SET (BAIT_SET, AMPLICON_SET, ...) used in this metrics collection run */
    public String PROBE_SET;

    /** The number of unique bases covered by the intervals of all probes in the probe set */
    public long PROBE_TERRITORY;

    /** The number of PF aligned probed bases that mapped to a baited region of the genome. */
    public long ON_PROBE_BASES;

    /** The number of PF aligned bases that mapped to within a fixed interval of a probed region, but not on a
     *  probed region. */
    public long NEAR_PROBE_BASES;

    /** The number of PF aligned bases that mapped to neither on or near a probe. */
    public long OFF_PROBE_BASES;

    /** The fraction of bases that map on or near a probe (ON_PROBE_BASES + NEAR_PROBE_BASES)/(ON_PROBE_BASES +
     * NEAR_PROBE_BASES + OFF_PROBE_BASES). */
    public double PCT_SELECTED_BASES;

    /** The fraction of aligned PF bases that mapped neither on or near a probe, OFF_PROBE_BASES/(ON_PROBE_BASES +
     *  NEAR_PROBE_BASES + OFF_PROBE_BASES). */
    public double PCT_OFF_PROBE;

    /** The fraction of on+near probe bases that are on as opposed to near, ON_PROBE_BASES/(ON_PROBE_BASES +
     * NEAR_PROBE_BASES). */
    public double ON_PROBE_VS_SELECTED;

    /** The mean coverage of all probes in the experiment, ON_PROBE_BASES/PROBE_TERRITORY. */
    public double MEAN_PROBE_COVERAGE;

    /** The fold by which the probed region has been amplified above genomic background,
     * (ON_PROBE_BASES/(ON_PROBE_BASES + NEAR_PROBE_BASES + OFF_PROBE_BASES))/(PROBE_TERRITORY/GENOME_SIZE) */
    public double FOLD_ENRICHMENT;
}


