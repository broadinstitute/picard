package picard.vcf;

/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Test for AccumulateVariantCallingMetrics
 *
 * @author Eric Banks
 */
public class AccumulateVariantCallingMetricsTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/vcf");

    @DataProvider(name = "shardDataProvider")
    public Object[][] shardDataProvider() {
        final File filePrefix1 = new File(TEST_DATA_DIR, "mergeTest.shard1");
        final File filePrefix2 = new File(TEST_DATA_DIR, "mergeTest.shard2");
        final File filePrefix3 = new File(TEST_DATA_DIR, "mergeTest.emptyShard");

        return new Object[][] {
                {Arrays.asList(filePrefix1, filePrefix2)},
                {Arrays.asList(filePrefix1, filePrefix2, filePrefix3)},
        };
    }


    @Test(dataProvider = "shardDataProvider")
    public void testMerge(final List<File> inputs) throws IOException {
        final File mergedFilePrefix = new File(TEST_DATA_DIR + "mergeTest");
        final File mergedSummaryFile = new File(mergedFilePrefix.getAbsolutePath() + ".variant_calling_summary_metrics");
        final File mergedDetailFile = new File(mergedFilePrefix.getAbsolutePath() + ".variant_calling_detail_metrics");
        mergedSummaryFile.deleteOnExit();
        mergedDetailFile.deleteOnExit();

        final AccumulateVariantCallingMetrics program = new AccumulateVariantCallingMetrics();
        program.INPUT = inputs;
        program.OUTPUT = mergedFilePrefix;

        Assert.assertEquals(program.doWork(), 0);

        final MetricsFile<CollectVariantCallingMetrics.VariantCallingDetailMetrics, Comparable<?>> detail = new MetricsFile<>();
        detail.read(new FileReader(mergedDetailFile));

        final MetricsFile<CollectVariantCallingMetrics.VariantCallingSummaryMetrics, Comparable<?>> summary = new MetricsFile<>();
        summary.read(new FileReader(mergedSummaryFile));

        checkResults(detail, summary);
    }

    private void checkResults(final MetricsFile<CollectVariantCallingMetrics.VariantCallingDetailMetrics, Comparable<?>> detail,
                              final MetricsFile<CollectVariantCallingMetrics.VariantCallingSummaryMetrics, Comparable<?>> summary) {

        int parsedDetail = 0;
        for (final CollectVariantCallingMetrics.VariantCallingDetailMetrics metrics : detail.getMetrics()) {
            if (metrics.SAMPLE_ALIAS.equals("FOO1")) {
                Assert.assertEquals(metrics.HET_HOMVAR_RATIO, 2.0);
                Assert.assertEquals(metrics.TOTAL_HET_DEPTH, 30);

                Assert.assertEquals(metrics.TOTAL_SNPS, 15);
                Assert.assertEquals(metrics.NUM_IN_DB_SNP, 10);
                Assert.assertEquals(metrics.NOVEL_SNPS, 5);
                Assert.assertEquals(metrics.FILTERED_SNPS, 7);

                Assert.assertEquals(metrics.PCT_DBSNP, 0.666667, 0.01);
                Assert.assertEquals(metrics.DBSNP_TITV, 2.333333, 0.01);
                Assert.assertEquals(metrics.NOVEL_TITV, 1.5, 0.01);

                Assert.assertEquals(metrics.TOTAL_INDELS, 9);
                Assert.assertEquals(metrics.NOVEL_INDELS, 3);
                Assert.assertEquals(metrics.FILTERED_INDELS, 12);
                Assert.assertEquals(metrics.NUM_IN_DB_SNP_INDELS, 6);

                Assert.assertEquals(metrics.PCT_DBSNP_INDELS, 0.666667, 0.01);
                Assert.assertEquals(metrics.DBSNP_INS_DEL_RATIO, 1.0, 0.01);
                Assert.assertEquals(metrics.NOVEL_INS_DEL_RATIO, 0.0, 0.01);

                Assert.assertEquals(metrics.SNP_REFERENCE_BIAS, 0.466667, 0.01);
                Assert.assertEquals(metrics.NUM_SINGLETONS, 10);
            } else if (metrics.SAMPLE_ALIAS.equals("FOO2")) {
                Assert.assertEquals(metrics.HET_HOMVAR_RATIO, 1.571429);
                Assert.assertEquals(metrics.TOTAL_HET_DEPTH, 33);

                Assert.assertEquals(metrics.TOTAL_SNPS, 18);
                Assert.assertEquals(metrics.NUM_IN_DB_SNP, 13);
                Assert.assertEquals(metrics.NOVEL_SNPS, 5);
                Assert.assertEquals(metrics.FILTERED_SNPS, 5);

                Assert.assertEquals(metrics.PCT_DBSNP, 0.722222, 0.01);
                Assert.assertEquals(metrics.DBSNP_TITV, 2.25, 0.01);
                Assert.assertEquals(metrics.NOVEL_TITV, 0.666667, 0.01);

                Assert.assertEquals(metrics.TOTAL_INDELS, 6);
                Assert.assertEquals(metrics.NOVEL_INDELS, 3);
                Assert.assertEquals(metrics.FILTERED_INDELS, 6);
                Assert.assertEquals(metrics.NUM_IN_DB_SNP_INDELS, 3);

                Assert.assertEquals(metrics.PCT_DBSNP_INDELS, 0.5, 0.01);
                Assert.assertEquals(metrics.DBSNP_INS_DEL_RATIO, 0.5, 0.01);
                Assert.assertEquals(metrics.NOVEL_INS_DEL_RATIO, 0.5, 0.01);

                Assert.assertEquals(metrics.SNP_REFERENCE_BIAS, 0.696969, 0.01);
                Assert.assertEquals(metrics.NUM_SINGLETONS, 9);
            } else {
                Assert.assertTrue(false, "Unexpected sample name in detailed metrics: " + metrics.SAMPLE_ALIAS);
            }
            parsedDetail++;
        }
        Assert.assertEquals(parsedDetail, 2, "Did not parse enough detail metrics.");

        boolean parsedSummary = false;
        for (final CollectVariantCallingMetrics.VariantCallingSummaryMetrics metrics : summary.getMetrics()) {
            Assert.assertEquals(metrics.TOTAL_SNPS, 33);
            Assert.assertEquals(metrics.NOVEL_SNPS, 10);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 23);
            Assert.assertEquals(metrics.FILTERED_SNPS, 12);

            Assert.assertEquals(metrics.PCT_DBSNP, 0.696969, 0.01);
            Assert.assertEquals(metrics.DBSNP_TITV, 2.285714, 0.01);
            Assert.assertEquals(metrics.NOVEL_TITV, 1.0, 0.01);

            Assert.assertEquals(metrics.TOTAL_INDELS, 15);
            Assert.assertEquals(metrics.NOVEL_INDELS, 6);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_INDELS, 9);
            Assert.assertEquals(metrics.FILTERED_INDELS, 18);

            Assert.assertEquals(metrics.PCT_DBSNP_INDELS, 0.6, 0.01);
            Assert.assertEquals(metrics.DBSNP_INS_DEL_RATIO, 0.8, 0.01);
            Assert.assertEquals(metrics.NOVEL_INS_DEL_RATIO, 0.2, 0.01);

            Assert.assertEquals(metrics.SNP_REFERENCE_BIAS, 0.587302, 0.01);
            Assert.assertEquals(metrics.NUM_SINGLETONS, 19);

            parsedSummary = true;
        }

        Assert.assertTrue(parsedSummary, "Did not parse summary metrics.");
    }
}
