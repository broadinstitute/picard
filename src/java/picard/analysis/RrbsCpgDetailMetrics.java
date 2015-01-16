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
 * Holds information about CpG sites encountered for RRBS processing QC
 * @author jgentry
 */
public final class RrbsCpgDetailMetrics extends MultilevelMetrics {
	/** Sequence the CpG is seen in */
	public String SEQUENCE_NAME;
	/** Position within the sequence of the CpG site */
	public Integer POSITION;
	/** Number of times this CpG site was encountered */
	public Integer TOTAL_SITES;
	/** Number of times this CpG site was converted (TG for + strand, CA for - strand) */
	public Integer CONVERTED_SITES;
	/** TOTAL_BASES / CONVERTED_BASES */
	public Double PCT_CONVERTED;
}
