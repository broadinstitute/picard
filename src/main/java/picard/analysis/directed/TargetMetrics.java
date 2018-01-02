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

import picard.metrics.MultilevelMetrics;

/**
 * TargetMetrics, are metrics to measure how well we hit specific targets (or baits) when using a targeted sequencing process like hybrid selection
 * or Targeted PCR Techniques (TSCA).  TargetMetrics at the moment are the metrics that are shared by both HybridSelection and TargetedPcrMetrics.
 */
public class TargetMetrics extends MultilevelMetrics {
    /**  The name of the PROBE_SET (BAIT_SET, AMPLICON_SET, ...) used in this metrics collection run */
    public String PROBE_SET;

    /** The number of unique bases covered by the intervals of all probes in the probe set */
    public long PROBE_TERRITORY;

    /** The number of unique bases covered by the intervals of all targets that should be covered */
    public long TARGET_TERRITORY;

    /** The number of bases in the reference genome used for alignment. */
    public long GENOME_SIZE;

    /** The total number of reads in the SAM or BAM file examined. */
    public long TOTAL_READS;

    /** The number of passing filter reads (PF). */
    public long PF_READS;

    /** The number of bases in the PF_READS of a SAM or BAM file */
    public long PF_BASES;

    /** The number of PF_READS that are not marked as duplicates. */
    public long PF_UNIQUE_READS;

    /** Tracks the number of read pairs that we see that are PF (used to calculate library size) */
    public long PF_SELECTED_PAIRS;

    /** Tracks the number of unique PF_SELECTED_PAIRS we see (used to calc library size) */
    public long PF_SELECTED_UNIQUE_PAIRS;

    /** The number of PF_UNIQUE_READS that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_READS_ALIGNED;

    /** The number of PF_BASES that are aligned with mapping score > 0 to the reference genome. */
    public long PF_BASES_ALIGNED;

    /** The number of PF unique bases that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_BASES_ALIGNED;

    /** The number of PF aligned probed bases that mapped to a baited region of the genome. */
    public long ON_PROBE_BASES;

    /** The number of PF aligned bases that mapped to within a fixed interval of a probed region, but not on a
     *  baited region. */
    public long NEAR_PROBE_BASES;

    /** The number of PF aligned bases that mapped to neither on or near a probe. */
    public long OFF_PROBE_BASES;

    /** The number of PF aligned bases that mapped to a targeted region of the genome. */
    public long ON_TARGET_BASES;

    /** The number of PF aligned bases that are mapped in pair to a targeted region of the genome. */
    public long ON_TARGET_FROM_PAIR_BASES;

    //metrics below here are derived after collection

    /** The fraction of reads passing filter, PF_READS/TOTAL_READS.   */
    public double PCT_PF_READS;

    /** The fraction of unique reads passing filter, PF_UNIQUE_READS/TOTAL_READS. */
    public double PCT_PF_UQ_READS;

    /** The fraction of unique reads passing filter that align to the reference,
     * PF_UQ_READS_ALIGNED/PF_UNIQUE_READS. */
    public double PCT_PF_UQ_READS_ALIGNED;

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

    /** The mean coverage of targets. */
    public double MEAN_TARGET_COVERAGE;

    /** The median coverage of targets. */
    public double MEDIAN_TARGET_COVERAGE;

    /** The maximum coverage of reads that mapped to target regions of an experiment. */
    public long MAX_TARGET_COVERAGE;

    /** The fraction of targets that did not reach coverage=1 over any base. */
    public double ZERO_CVG_TARGETS_PCT;

    /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
    public double PCT_EXC_DUPE;

    /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality. */
    public double PCT_EXC_MAPQ;

    /** The fraction of aligned bases that were filtered out because they were of low base quality. */
    public double PCT_EXC_BASEQ;

    /** The fraction of aligned bases that were filtered out because they were the second observation from
     *  an insert with overlapping reads. */
    public double PCT_EXC_OVERLAP;

    /** The fraction of aligned bases that were filtered out because they did not align over a target base. */
    public double PCT_EXC_OFF_TARGET;

    /**
     * The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to
     * the mean coverage level in those targets.
     */
    public double FOLD_80_BASE_PENALTY;

    /** The fraction of all target bases achieving 1X or greater coverage. */
    public double PCT_TARGET_BASES_1X;
    /** The fraction of all target bases achieving 2X or greater coverage. */
    public double PCT_TARGET_BASES_2X;
    /** The fraction of all target bases achieving 10X or greater coverage. */
    public double PCT_TARGET_BASES_10X;
    /** The fraction of all target bases achieving 20X or greater coverage. */
    public double PCT_TARGET_BASES_20X;
    /** The fraction of all target bases achieving 30X or greater coverage. */
    public double PCT_TARGET_BASES_30X;
    /** The fraction of all target bases achieving 40X or greater coverage. */
    public double PCT_TARGET_BASES_40X;
    /** The fraction of all target bases achieving 50X or greater coverage. */
    public double PCT_TARGET_BASES_50X;
    /** The fraction of all target bases achieving 100X or greater coverage. */
    public double PCT_TARGET_BASES_100X;

    /**
     * A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50]
     * we calculate a = % of target territory, and b = % of aligned reads aligned to these targets.
     * AT DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total
     * reads that should have mapped to GC<=50% regions mapped elsewhere.
     */
    public double AT_DROPOUT;

    /**
     * A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100]
     * we calculate a = % of target territory, and b = % of aligned reads aligned to these targets.
     * GC DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total
     * reads that should have mapped to GC>=50% regions mapped elsewhere.
     */
    public double GC_DROPOUT;

    /** The theoretical HET SNP sensitivity. */
    public double HET_SNP_SENSITIVITY;

    /** The Phred Scaled Q Score of the theoretical HET SNP sensitivity. */
    public double HET_SNP_Q;
}


