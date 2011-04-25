/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

package net.sf.picard.analysis;

import net.sf.picard.metrics.MetricBase;

/**
 * High level metrics about the alignment of reads within a SAM file, produced by
 * the CollectAlignmentSummaryMetrics program and usually stored in a file with
 * the extension ".alignment_summary_metrics".
 */
public class AlignmentSummaryMetrics extends MetricBase {
    public enum Category { UNPAIRED, FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR }

    /**
     * One of either UNPAIRED (for a fragment run), FIRST_OF_PAIR when metrics are for only the
     * first read in a paired run, SECOND_OF_PAIR when the metrics are for only the second read
     * in a paired run or PAIR when the metrics are aggregeted for both first and second reads
     * in a pair.
     */
    public Category CATEGORY;

    /**
     * The total number of reads including all PF and non-PF reads. When CATEGORY equals PAIR
     * this value will be 2x the number of clusters.
     */
    public long TOTAL_READS;

    /** The number of PF reads where PF is defined as passing Illumina's filter. */
    public long PF_READS;

    /** The percentage of reads that are PF (PF_READS / TOTAL_READS) */
    public double PCT_PF_READS;

    /**
     * The number of PF reads that are marked as noise reads.  A noise read is one which is composed
     * entirey of A bases and/or N bases. These reads are marked as they are usually artifactual and
     * are of no use in downstream analysis.
     */
    public long PF_NOISE_READS;

    /**
     * The number of PF reads that were aligned to the reference sequence. This includes reads that
     * aligned with low quality (i.e. their alignments are ambiguous).
     */
    public long PF_READS_ALIGNED;

    /**
     * The percentage of PF reads that aligned to the reference sequence. PF_READS_ALIGNED / PF_READS
     */
    public double PCT_PF_READS_ALIGNED;
    
    /**
     * The total number of aligned bases, in all mapped PF reads, that are aligned to the reference sequence.
     */
    public long PF_ALIGNED_BASES;

    /**
     * The number of PF reads that were aligned to the reference sequence with a mapping quality of
     * Q20 or higher signifying that the aligner estimates a 1/100 (or smaller) chance that the
     * alignment is wrong.
     */
    public long PF_HQ_ALIGNED_READS;

    /**
     * The number of bases aligned to the reference sequence in reads that were mapped at high
     * quality.  Will usually approximate PF_HQ_ALIGNED_READS * READ_LENGTH but may differ when
     * either mixed read lengths are present or many reads are aligned with gaps.
     */
    public long PF_HQ_ALIGNED_BASES;

    /**
     * The subest of PF_HQ_ALIGNED_BASES where the base call quality was Q20 or higher.
     */
    public long PF_HQ_ALIGNED_Q20_BASES;

    /**
     * The median number of mismatches versus the reference sequence in reads that were aligned
     * to the reference at high quality (i.e. PF_HQ_ALIGNED READS).
     */
    public double PF_HQ_MEDIAN_MISMATCHES;

    /**
     * The rate of bases mismatching the reference for all bases aligned to the reference sequence.
     */
    public double PF_MISMATCH_RATE;

    /**
     * The percentage of bases that mismatch the reference in PF HQ aligned reads.
     */
    public double PF_HQ_ERROR_RATE;

    /**
     * The mean read length of the set of reads examined.  When looking at the data for a single lane with
     * equal length reads this number is just the read length.  When looking at data for merged lanes with
     * differing read lengths this is the mean read length of all reads.
     */
    public double MEAN_READ_LENGTH;

    /**
     * The number of aligned reads who's mate pair was also aligned to the reference.
     */
    public long READS_ALIGNED_IN_PAIRS;

    /**
     * The percentage of reads who's mate pair was also aligned to the reference.
     * READS_ALIGNED_IN_PAIRS / PF_READS_ALIGNED
     */
    public double PCT_READS_ALIGNED_IN_PAIRS;

    /**
     * The number of instrument cycles in which 80% or more of base calls were no-calls.
     */
    public long BAD_CYCLES;

    /**
     * The number of PF reads aligned to the positive strand of the genome divided by the number of
     * PF reads aligned to the genome.
     */
    public double STRAND_BALANCE;

    /**
     * The percentage of reads that map outside of a maximum insert size (usually 100kb) or that have
     * the two ends mapping to different chromosomes.
     */
    public double PCT_CHIMERAS;

    /**
     * The percentage of PF reads that are unaligned and match to a known adapter sequence right from the
     * start of the read.
     */
    public double PCT_ADAPTER;
}
