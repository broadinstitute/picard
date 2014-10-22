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

import htsjdk.samtools.metrics.MetricBase;

import java.util.List;

/**
 * Holds per-MetricAccumulationLevel metric information for the RRBS metrics. Required as the MultiLevelCollector
 * is designed around having a single metrics object and we have two being calculated so RrbsMetricsCollector builds
 * this object which can be teased apart downstream
 *
 * NB: This is purely for internal use, if used as a proper metric object it likely won't do what you want it to
 *
 * @author jgentry@broadinstitute.org
 */
class RrbsMetrics extends MetricBase {
	private final RrbsSummaryMetrics summaryMetrics;
	private final List<RrbsCpgDetailMetrics> detailMetrics;

	public RrbsMetrics(final RrbsSummaryMetrics summaryMetrics, final List<RrbsCpgDetailMetrics> detailMetrics) {
		this.summaryMetrics = summaryMetrics;
		this.detailMetrics = detailMetrics;
	}

	public List<RrbsCpgDetailMetrics> getDetailMetrics() {
		return detailMetrics;
	}

	public RrbsSummaryMetrics getSummaryMetrics() {
		return summaryMetrics;
	}
}
