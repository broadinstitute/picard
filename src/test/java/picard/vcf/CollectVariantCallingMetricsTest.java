package picard.vcf;

/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

/**
 * Test for CollectVariantCallingMetrics
 *
 * @author Joel Thibault (thibault at broadinstitute dot org)
 */
public class CollectVariantCallingMetricsTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/vcf");

    @Test
    public void testMetricsTiny() throws IOException {
        final File dbSnpFile = new File(TEST_DATA_DIR, "mini.dbsnp.vcf");
        final File vcfFile = new File(TEST_DATA_DIR, "mini.vcf");

        final File outFile = new File(TEST_DATA_DIR, "vcmetrics_tiny");
        final File summaryFile = new File(TEST_DATA_DIR, "vcmetrics_tiny.variant_calling_summary_metrics");
        final File detailFile = new File(TEST_DATA_DIR, "vcmetrics_tiny.variant_calling_detail_metrics");

        summaryFile.deleteOnExit();
        detailFile.deleteOnExit();

        final CollectVariantCallingMetrics program = new CollectVariantCallingMetrics();
        program.INPUT = vcfFile;
        program.DBSNP = dbSnpFile;
        program.OUTPUT = outFile;

        Assert.assertEquals(program.doWork(), 0);

        final MetricsFile<CollectVariantCallingMetrics.VariantCallingSummaryMetrics, Comparable<?>> summary = new MetricsFile<>();
        summary.read(new FileReader(summaryFile));

        boolean parsedSummary = false;
        for (final CollectVariantCallingMetrics.VariantCallingSummaryMetrics metrics : summary.getMetrics()) {
            Assert.assertEquals(metrics.TOTAL_SNPS, 597);
            Assert.assertEquals(metrics.NOVEL_SNPS, 265);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 332);

            Assert.assertEquals(metrics.PCT_DBSNP, 0.5561140179634094, 0.01);
            Assert.assertEquals(metrics.DBSNP_TITV, 3.955224, 0.01);
            Assert.assertEquals(metrics.NOVEL_TITV, 3.206349, 0.01);

            Assert.assertEquals(metrics.TOTAL_INDELS, 29);
            Assert.assertEquals(metrics.NOVEL_INDELS, 11);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_INDELS, 18);

            Assert.assertEquals(metrics.PCT_DBSNP_INDELS, 0.62069, 0.01);
            Assert.assertEquals(metrics.DBSNP_INS_DEL_RATIO, 0.125, 0.01);
            Assert.assertEquals(metrics.NOVEL_INS_DEL_RATIO, 0.375, 0.01);
            Assert.assertEquals(metrics.NUM_SINGLETONS, 245);

            parsedSummary = true;
        }

        Assert.assertTrue(parsedSummary, "Did not parse summary metrics.");

        final MetricsFile<CollectVariantCallingMetrics.VariantCallingDetailMetrics, Comparable<?>> detail = new MetricsFile<>();
        detail.read(new FileReader(detailFile));
        final List<CollectVariantCallingMetrics.VariantCallingDetailMetrics> detailMetrics = detail.getMetrics();
        detail.getMetrics().stream().filter(metrics -> metrics.SAMPLE_ALIAS.equals("HG00160")).forEach(metrics -> {
            Assert.assertEquals(metrics.HET_HOMVAR_RATIO, 0.72549, 0.0001);
            Assert.assertEquals(metrics.TOTAL_SNPS, 81);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 44);
            Assert.assertEquals(metrics.NOVEL_SNPS, 37);
            Assert.assertEquals(metrics.PCT_DBSNP, 0.543210, 0.01);
            Assert.assertEquals(metrics.DBSNP_TITV, 6.333333, 0.01);
            Assert.assertEquals(metrics.NOVEL_TITV, 2.7, 0.01);
            Assert.assertEquals(metrics.TOTAL_INDELS, 6);
            Assert.assertEquals(metrics.NOVEL_INDELS, 3);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_INDELS, 3);
            Assert.assertEquals(metrics.PCT_DBSNP_INDELS, 0.5, 0.01);
            Assert.assertEquals(metrics.DBSNP_INS_DEL_RATIO, 0.0, 0.01);
            Assert.assertEquals(metrics.NOVEL_INS_DEL_RATIO, 0.0, 0.01);
            Assert.assertEquals(metrics.TOTAL_MULTIALLELIC_SNPS, 0.0, 0.01);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_MULTIALLELIC, 0, 0.01);
            Assert.assertEquals(metrics.TOTAL_COMPLEX_INDELS, 1.0, 0.01);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_COMPLEX_INDELS, 0, 0.01);
            Assert.assertEquals(metrics.SNP_REFERENCE_BIAS, 0.510204, 0.01);
            Assert.assertEquals(metrics.NUM_SINGLETONS, 3);
        });

        Assert.assertEquals(detailMetrics.size(), 50, "Did not parse the desired number of detail metrics.");
    }


    @Test
    public void testMetricsTinyGVCF() throws IOException {
        final File dbSnpFile = new File(TEST_DATA_DIR, "mini.dbsnp.vcf");
        final File vcfFile = new File(TEST_DATA_DIR, "mini_gvcf.vcf");

        final File outFile = new File(TEST_DATA_DIR, "vcmetrics_tiny_gvcf");
        final File summaryFile = new File(outFile+".variant_calling_summary_metrics");
        final File detailFile = new File(outFile+".variant_calling_detail_metrics");

        summaryFile.deleteOnExit();
        detailFile.deleteOnExit();

        final CollectVariantCallingMetrics program = new CollectVariantCallingMetrics();
        program.INPUT = vcfFile;
        program.DBSNP = dbSnpFile;
        program.OUTPUT = outFile;
        program.GVCF_INPUT = true;
        Assert.assertEquals(program.doWork(), 0);

        final MetricsFile<CollectVariantCallingMetrics.VariantCallingSummaryMetrics, Comparable<?>> summary = new MetricsFile<>();
        summary.read(new FileReader(summaryFile));

        boolean parsedSummary = false;
        for (final CollectVariantCallingMetrics.VariantCallingSummaryMetrics metrics : summary.getMetrics()) {
            Assert.assertEquals(metrics.TOTAL_SNPS, 20);
            Assert.assertEquals(metrics.NOVEL_SNPS, 19);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 1);
            Assert.assertEquals(metrics.FILTERED_SNPS, 0);

            Assert.assertEquals(metrics.PCT_DBSNP, 0.05, 0.001);
            Assert.assertEquals(metrics.DBSNP_TITV, 0D, 0.01);
            Assert.assertEquals(metrics.NOVEL_TITV, 12D/(19-12), 0.01);

            Assert.assertEquals(metrics.TOTAL_INDELS, 7);
            Assert.assertEquals(metrics.NOVEL_INDELS, 7);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_INDELS, 0);
            Assert.assertEquals(metrics.NOVEL_INS_DEL_RATIO,3/4D,0.01);
            Assert.assertEquals(metrics.PCT_DBSNP_INDELS, 0, 0.01);
            Assert.assertEquals(metrics.DBSNP_INS_DEL_RATIO, 0, 0.01);
            Assert.assertEquals(metrics.NUM_SINGLETONS, 8);

            parsedSummary = true;
        }

        Assert.assertTrue(parsedSummary, "Did not parse summary metrics.");

        final MetricsFile<CollectVariantCallingMetrics.VariantCallingDetailMetrics, Comparable<?>> detail = new MetricsFile<>();
        detail.read(new FileReader(detailFile));
        final List<CollectVariantCallingMetrics.VariantCallingDetailMetrics> detailMetrics = detail.getMetrics();
        detail.getMetrics().stream().filter(metrics -> metrics.SAMPLE_ALIAS.equals("HG00160")).forEach(metrics -> {
            Assert.assertEquals(metrics.HET_HOMVAR_RATIO, .6, 0.0001);
            Assert.assertEquals(metrics.TOTAL_SNPS, 20);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 1);
            Assert.assertEquals(metrics.NOVEL_SNPS, 19);
            Assert.assertEquals(metrics.PCT_DBSNP, 1D/20, 0.01);
            Assert.assertEquals(metrics.DBSNP_TITV, 0D, 0.01);
            Assert.assertEquals(metrics.NOVEL_TITV, 12D/(19-12), 0.01);
            Assert.assertEquals(metrics.TOTAL_INDELS, 7);
            Assert.assertEquals(metrics.NOVEL_INDELS, 7);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_INDELS, 0);
            Assert.assertEquals(metrics.PCT_DBSNP_INDELS, 0, 0.01);
            Assert.assertEquals(metrics.DBSNP_INS_DEL_RATIO, 0.0, 0.01);
            Assert.assertEquals(metrics.NOVEL_INS_DEL_RATIO, 3/4D, 0.01);
            Assert.assertEquals(metrics.TOTAL_MULTIALLELIC_SNPS, 0.0, 0.01);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_MULTIALLELIC, 0, 0.01);
            Assert.assertEquals(metrics.TOTAL_COMPLEX_INDELS, 0, 0.01);
            Assert.assertEquals(metrics.NUM_IN_DB_SNP_COMPLEX_INDELS, 0, 0.01);
            Assert.assertEquals(metrics.NUM_SINGLETONS, 8);
        });

        Assert.assertEquals(detailMetrics.size(), 1, "Did not parse the expected number of detail metrics.");
    }
}
