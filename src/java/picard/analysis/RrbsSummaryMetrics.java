/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
 * Holds summary statistics from RRBS processing QC
 *
 * @author jgentry
 */
public final class RrbsSummaryMetrics extends MultilevelMetrics {
	/** Number of mapped reads processed */
	public Integer READS_ALIGNED;
	/** Number of times a non-CpG cytosine was encountered */
	public Integer NON_CPG_BASES;
	/** Number of times a non-CpG cytosine was converted (C->T for +, G->A for -) */
	public Integer NON_CPG_CONVERTED_BASES;
	/** NON_CPG_BASES / NON_CPG_CONVERTED_BASES */
	public Double PCT_NON_CPG_BASES_CONVERTED;
	/** Number of CpG sites encountered */
	public Integer CPG_BASES_SEEN;
	/** Number of CpG sites that were converted (TG for +, CA for -) */
	public Integer CPG_BASES_CONVERTED;
	/** CPG_BASES_SEEN / CPG_BASES_CONVERTED */
	public Double PCT_CPG_BASES_CONVERTED;
	/** Mean coverage of CpG sites */
	public Double MEAN_CPG_COVERAGE;
	/** Median coverage of CpG sites */
	public Integer MEDIAN_CPG_COVERAGE;
	/** Number of reads discarded for having no CpG sites */
	public Integer READS_WITH_NO_CPG;
	/** Number of reads discarded due to being too short */
	public Integer READS_IGNORED_SHORT;
	/** Number of reads discarded for exceeding the mismatch threshold */
	public Integer READS_IGNORED_MISMATCHES;
}
