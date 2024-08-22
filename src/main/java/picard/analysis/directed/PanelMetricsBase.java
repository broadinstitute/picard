/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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
 * A base class for Metrics for targeted panels. Metrics for library construction protocols based on
 * PCR amplicons/probes and Hybrid Selection baits are derived from this class.
 */
public class PanelMetricsBase extends MultilevelMetrics {
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

    /** The number of PF_UNIQUE_READS that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_READS_ALIGNED;

    /** The number of PF_BASES that are aligned with mapping score > 0 to the reference genome. */
    public long PF_BASES_ALIGNED;

    /** The number of PF unique bases that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_BASES_ALIGNED;

    /** The number of PF aligned bases that mapped to a targeted region of the genome. */
    public long ON_TARGET_BASES;

    // Metrics below here are derived after collection

    /** The fraction of reads passing filter, PF_READS/TOTAL_READS.   */
    public double PCT_PF_READS;

    /** The fraction of unique reads passing filter, PF_UNIQUE_READS/TOTAL_READS. */
    public double PCT_PF_UQ_READS;

    /** The fraction of unique reads passing filter that align to the reference,
     * PF_UQ_READS_ALIGNED/PF_UNIQUE_READS. */
    public double PCT_PF_UQ_READS_ALIGNED;

    /** The mean coverage of targets. */
    public double MEAN_TARGET_COVERAGE;

    /** The median coverage of targets. */
    public double MEDIAN_TARGET_COVERAGE;

    /** The maximum coverage of reads that mapped to target regions of an experiment. */
    public long MAX_TARGET_COVERAGE;

    /** The minimum coverage of reads that mapped to target regions of an experiment. */
    public long MIN_TARGET_COVERAGE;

    /** The fraction of targets that did not reach coverage=1 over any base. */
    public double ZERO_CVG_TARGETS_PCT;

    /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
    public double PCT_EXC_DUPE;

    /** The fraction of aligned bases in reads that have MQ=0 and whose 5' end consists of adapter sequence. */
    public double PCT_EXC_ADAPTER;

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
    /** The fraction of all target bases achieving 250X or greater coverage. */
    public double PCT_TARGET_BASES_250X;
    /** The fraction of all target bases achieving 500X or greater coverage. */
    public double PCT_TARGET_BASES_500X;
    /** The fraction of all target bases achieving 1000X or greater coverage. */
    public double PCT_TARGET_BASES_1000X;
    /** The fraction of all target bases achieving 2500X or greater coverage. */
    public double PCT_TARGET_BASES_2500X;
    /** The fraction of all target bases achieving 5000X or greater coverage. */
    public double PCT_TARGET_BASES_5000X;
    /** The fraction of all target bases achieving 10000X or greater coverage. */
    public double PCT_TARGET_BASES_10000X;
    /** The fraction of all target bases achieving 25000X or greater coverage. */
    public double PCT_TARGET_BASES_25000X;
    /** The fraction of all target bases achieving 50000X or greater coverage. */
    public double PCT_TARGET_BASES_50000X;
    /** The fraction of all target bases achieving 100000X or greater coverage. */
    public double PCT_TARGET_BASES_100000X;

    /**
     * A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50]
     * we calculate a = % of target territory, and b = % of aligned reads aligned to these targets.
     * AT_DROPOUT is then sum(a-b if a-b > 0 else 0). E.g. if the value is 5%, this implies that 5% of total
     * reads that should have mapped to GC<=50% regions mapped elsewhere.
     */
    public double AT_DROPOUT;

    /**
     * A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100]
     * we calculate a = % of target territory, and b = % of aligned reads aligned to these targets.
     * GC_DROPOUT is then sum(a-b if a-b > 0 else 0). E.g. if the value is 5% this implies that 5% of total
     * reads that should have mapped to GC>=50% regions mapped elsewhere.
     */
    public double GC_DROPOUT;

    /** The theoretical HET SNP sensitivity. */
    public double HET_SNP_SENSITIVITY;

    /** The Phred Scaled Q Score of the theoretical HET SNP sensitivity. */
    public double HET_SNP_Q;
}


