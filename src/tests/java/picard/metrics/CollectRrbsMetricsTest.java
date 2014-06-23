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

package picard.metrics;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
import picard.analysis.CollectRrbsMetrics;
import picard.analysis.RrbsSummaryMetrics;

import java.io.File;
import java.io.FileReader;
import java.lang.Exception;import java.lang.Integer;import java.lang.String;import java.util.ArrayList;
import java.util.List;

/**
 * @author jgentry@broadinstitute.org
 */

public class CollectRrbsMetricsTest {
	public static final String CHR_M_SAM = "testdata/picard/metrics/chrMReads.sam";
	public static final String CHR_M_REFERENCE ="testdata/picard/metrics/chrM.reference.fasta";

	private File rootTestDir;

	@BeforeTest
	private void setUp() throws Exception {
		rootTestDir = File.createTempFile("crmt.", ".tmp");
		Assert.assertTrue(rootTestDir.delete());
		Assert.assertTrue(rootTestDir.mkdir());
	}

	@AfterTest
	private void tearDown() {
		IOUtil.deleteDirectoryTree(rootTestDir);
	}

	@Test
	public void chrMReads() throws Exception {
		final MetricsFile<RrbsSummaryMetrics, ?> metricsFile = getSummaryFile(CHR_M_SAM, CHR_M_REFERENCE, rootTestDir + "/READ_TEST", new ArrayList<String>());
		final RrbsSummaryMetrics metrics = metricsFile.getMetrics().get(0);
		Assert.assertEquals(metrics.READS_ALIGNED.intValue(), 5);
		Assert.assertEquals(metrics.NON_CPG_BASES.intValue(), 15);
		Assert.assertEquals(metrics.NON_CPG_CONVERTED_BASES.intValue(), 11);
		Assert.assertEquals(metrics.PCT_NON_CPG_BASES_CONVERTED, 0.733333);
		Assert.assertEquals(metrics.CPG_BASES_SEEN.intValue(), 5);
		Assert.assertEquals(metrics.CPG_BASES_CONVERTED.intValue(), 1);
		Assert.assertEquals(metrics.PCT_CPG_BASES_CONVERTED, 0.2);
		Assert.assertEquals(metrics.MEAN_CPG_COVERAGE, 1.666667);
		Assert.assertEquals(metrics.MEDIAN_CPG_COVERAGE.intValue(), 2);
		Assert.assertEquals(metrics.READS_WITH_NO_CPG.intValue(), 1);
		Assert.assertEquals(metrics.READS_IGNORED_SHORT.intValue(), 1);
		Assert.assertEquals(metrics.READS_IGNORED_MISMATCHES.intValue(), 1);
	}

	private MetricsFile<RrbsSummaryMetrics, ?> getSummaryFile(final String input, final String reference, final String prefix,
															  final List<String> sequences) throws Exception {
		final List<String> argList = new ArrayList<String>();
		argList.add("INPUT=" + input);
		argList.add("METRICS_FILE_PREFIX=" + prefix);
		argList.add("REFERENCE=" + reference);
		for (final String sequence : sequences) {
			argList.add("SEQUENCE_NAMES=" + sequence);
		}

		final String[] args = new String[argList.size()];
		argList.toArray(args);

		Assert.assertEquals(new CollectRrbsMetrics().instanceMain(args), 0);

		final MetricsFile<RrbsSummaryMetrics, ?> retVal = new MetricsFile<RrbsSummaryMetrics, Integer>();
		retVal.read(new FileReader(prefix + ".rrbs_summary_metrics"));
		return retVal;
	}


}
