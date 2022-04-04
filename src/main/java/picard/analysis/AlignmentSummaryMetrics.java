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

package picard.analysis;

import picard.metrics.MultilevelMetrics;

/**
 * High level metrics about the alignment of reads within a SAM file, produced by
 * the CollectAlignmentSummaryMetrics program and usually stored in a file with
 * the extension ".alignment_summary_metrics".
 */
public class AlignmentSummaryMetrics extends MultilevelMetrics {
    public enum Category {UNPAIRED, FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR}

    /**
     * One of either UNPAIRED (for a fragment run), FIRST_OF_PAIR when metrics are for only the
     * first read in a paired run, SECOND_OF_PAIR when the metrics are for only the second read
     * in a paired run or PAIR when the metrics are aggregated for both first and second reads
     * in a pair.
     */
    public Category CATEGORY;

    /**
     * The total number of reads including all PF and non-PF reads. When CATEGORY equals PAIR
     * this value will be 2x the number of clusters.
     */
    public long TOTAL_READS;

    /**
     * The number of PF reads where PF is defined as passing Illumina's filter.
     */
    public long PF_READS;

    /**
     * The fraction of reads that are PF (PF_READS / TOTAL_READS)
     */
    public double PCT_PF_READS;

    /**
     * The number of PF reads that are marked as noise reads.  A noise read is one which is composed
     * entirely of A bases and/or N bases. These reads are marked as they are usually artifactual and
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
     * The subset of PF_HQ_ALIGNED_BASES where the base call quality was Q20 or higher.
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
     * The fraction of bases that mismatch the reference in PF HQ aligned reads.
     */
    public double PF_HQ_ERROR_RATE;

    /**
     * The number of insertion and deletion events per 100 aligned bases.  Uses the number of events
     * as the numerator, not the number of inserted or deleted bases.
     */
    public double PF_INDEL_RATE;

    /**
     * The mean read length of the set of reads examined.  When looking at the data for a single lane with
     * equal length reads this number is just the read length.  When looking at data for merged lanes with
     * differing read lengths this is the mean read length of all reads. Computed using all read lengths
     * including clipped bases.
     */
    public double MEAN_READ_LENGTH;

    /** The standard deviation of the read lengths. Computed using all read lengths including clipped bases. */
    public double SD_READ_LENGTH;

    /**
     * The median read length of the set of reads examined.  When looking at the data for a single lane with
     * equal length reads this number is just the read length.  When looking at data for merged lanes with
     * differing read lengths this is the median read length of all reads. Computed using all bases in reads,
     * including clipped bases.
     */
    public double MEDIAN_READ_LENGTH;

    /**
     * The median absolute deviation of the distribution of all read lengths.  If the distribution is
     * essentially normal then the standard deviation can be estimated as ~1.4826 * MAD. Computed using all
     * read lengths including clipped bases.
     */
    public double MAD_READ_LENGTH;

    /** The minimum read length. Computed using all read lengths including clipped bases. */
    public double MIN_READ_LENGTH;

    /** The maximum read length. Computed using all read lengths including clipped bases. */
    public double MAX_READ_LENGTH;

    /**
     * The mean aligned read length of the set of reads examined.  When looking at the data for a single lane with
     * equal length reads this number is just the read length.  When looking at data for merged lanes with
     * differing read lengths this is the mean read length of all reads. Clipped bases are not counted.
     */
    public double MEAN_ALIGNED_READ_LENGTH;

    /**
     * The number of aligned reads whose mate pair was also aligned to the reference.
     */
    public long READS_ALIGNED_IN_PAIRS;

    /**
     * The fraction of aligned reads whose mate pair was also aligned to the reference.
     * READS_ALIGNED_IN_PAIRS / PF_READS_ALIGNED
     */
    public double PCT_READS_ALIGNED_IN_PAIRS;

    /**
     * The number of (primary) aligned reads that are **not** "properly" aligned in pairs (as per SAM flag 0x2).
     */
    public long PF_READS_IMPROPER_PAIRS;

    /**
     * The fraction of (primary) reads that are *not* "properly" aligned in pairs (as per SAM flag 0x2).
     * PF_READS_IMPROPER_PAIRS / PF_READS_ALIGNED
     */
    public double PCT_PF_READS_IMPROPER_PAIRS;

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
     * The fraction of reads that map outside of a maximum insert size (usually 100kb) or that have
     * the two ends mapping to different chromosomes.
     */
    public double PCT_CHIMERAS;

    /**
     * The fraction of PF reads that are unaligned or aligned with MQ0 and match to a known adapter sequence right from the
     * start of the read (indication of adapter-dimer pairs).
     */
    public double PCT_ADAPTER;

    /**
     * the fraction of PF bases that are on (primary) aligned reads and are soft-clipped, as a fraction of the
     * PF_ALIGNED_BASES (even though these are not aligned!)
     */
    public double PCT_SOFTCLIP;

    /**
     * The fraction of PF bases that are (on primary, aligned reads and) hard-clipped, as a fraction of the
     * PF_ALIGNED_BASES (even though these are not aligned!)
     */
    public double PCT_HARDCLIP;

    /**
     * The average length of the soft-clipped bases at the 3' end of reads. This could be used as an estimate for
     * the amount by which the insert-size must be increased in order to obtain a significant reduction in bases
     * lost due to reading off the end of the insert.
     */
    public double AVG_POS_3PRIME_SOFTCLIP_LENGTH;
}
