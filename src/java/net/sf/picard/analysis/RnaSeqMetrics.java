/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
 * Metrics about the alignment of RNA-seq reads within a SAM file to genes, produced by the CollectRnaSeqMetrics
 * program and usually stored in a file with the extension ".rna_metrics".
 */
public class RnaSeqMetrics extends MetricBase{
    /** The total number of PF bases including non-aligned reads. */
    public long PF_BASES;

    /** The total number of aligned PF bases.  Non-primary alignments are not counted. */
    public long PF_ALIGNED_BASES;

    /** Number of bases in primary aligments that align to ribosomal sequence. */
    public long RIBOSOMAL_BASES;

    /** Number of bases in primary aligments that align to a non-UTR coding base for some gene, and not ribosomal sequence. */
    public long CODING_BASES;

    /** Number of bases in primary aligments that align to a UTR base for some gene, and not a coding base. */
    public long UTR_BASES;

    /** Number of bases in primary aligments that align to an intronic base for some gene, and not a coding or UTR base. */
    public long INTRONIC_BASES;

    /** Number of bases in primary aligments that do not align to any gene. */
    public long INTERGENIC_BASES;

    /**
     * Number of primary alignments that map to a sequence specified on command-line as IGNORED_SEQUENCE.  These are not
     * counted in PF_ALIGNED_BASES, CORRECT_STRAND_READS, INCORRECT_STRAND_READS, or any of the base-counting metrics.
     * These reads are counted in PF_BASES.
     */
    public long IGNORED_READS;

    /** Number of aligned reads that map to the correct strand.  0 if library is not strand-specific. */
    public long CORRECT_STRAND_READS;

    /** Number of aligned reads that map to the incorrect strand.  0 if library is not strand-specific. */
    public long INCORRECT_STRAND_READS;

    /** RIBOSOMAL_BASES / PF_ALIGNED_BASES */
    public double PCT_RIBOSOMAL_BASES;

    /** CODING_BASES / PF_ALIGNED_BASES */
    public double PCT_CODING_BASES;

    /** UTR_BASES / PF_ALIGNED_BASES */
    public double PCT_UTR_BASES;

    /** INTRONIC_BASES / PF_ALIGNED_BASES */
    public double PCT_INTRONIC_BASES;

    /** INTERGENIC_BASES / PF_ALIGNED_BASES */
    public double PCT_INTERGENIC_BASES;

    /** PCT_UTR_BASES + PCT_CODING_BASES */
    public double PCT_MRNA_BASES;

    /** The percentage of bases mapping to mRNA divided by the total number of PF bases. */
    public double PCT_USABLE_BASES;

    /** CORRECT_STRAND_READS/(CORRECT_STRAND_READS + INCORRECT_STRAND_READS).  0 if library is not strand-specific. */
    public double PCT_CORRECT_STRAND_READS;

    /** The median CV of coverage of the 1000 most highly expressed transcripts. */
    public double MEDIAN_CV_COVERAGE;

    /**
     * The median 5 prime bias of the 1000 most highly expressed transcripts, where 5 prime bias is calculated per
     * transcript as: mean coverage of the 5' most 100 bases divided by the mean coverage of the whole transcript.
     */
    public double MEDIAN_5PRIME_BIAS;

    /**
     * The median 3 prime bias of the 1000 most highly expressed transcripts, where 3 prime bias is calculated per
     * transcript as: mean coverage of the 3' most 100 bases divided by the mean coverage of the whole transcript.
     */
    public double MEDIAN_3PRIME_BIAS;

    /** The ratio of coverage at the 5' end of to the 3' end based on the 1000 most highly expressed transcripts. */
    public double MEDIAN_5PRIME_TO_3PRIME_BIAS;
}
