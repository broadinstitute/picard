package picard.arrays;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Iso8601Date;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.pedigree.Sex;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

public class CollectArraysVariantCallingMetricsTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/");
    private static final File TEST_DBSNP_DATA_DIR = new File("testdata/picard/vcf/");

    @DataProvider(name = "collectArraysVariantCallingMetricsTestProvider")
    public Object[][] testCollectArraysVariantCallingMetricsTestProvider() {
        return new Object[][]{
                {"7991775143_R01C01", 636, 1, true, "F", "F", true, "3.0.0", "1.0.0.0"},
                {"7991775143_R01C01_wo_optional_fields", 631, null, false, "M", "N", false, "2.0.0.137", null}
        };
    }

    @Test(dataProvider = "collectArraysVariantCallingMetricsTestProvider")
    public void testCollectArraysVariantCallingMetricsOverkillTest(final String chipWellBarcode,
                                                                   final int expectedNumCalls,
                                                                   final Integer expectedAnalysisVersion,
                                                                   final Boolean expectedIsZcalled,
                                                                   final String expectedAutocallGender,
                                                                   final String expectedReportedGender,
                                                                   final Boolean expectedGenderConcordancePF,
                                                                   final String expectedAutocallVersion,
                                                                   final String expectedZcallVersion) throws IOException {
        final File dbSnpFile = new File(TEST_DBSNP_DATA_DIR, "mini.dbsnp.vcf");
        final File vcfFile = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, chipWellBarcode + ".vcf"), chipWellBarcode);
        vcfFile.deleteOnExit();

        final File outputBaseFile = File.createTempFile("collectArraysVariantCallingMetrics", "");
        final File summaryMetricsFile = new File(outputBaseFile.getAbsolutePath() + "." + CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics.getFileExtension());
        final File detailMetricsFile  = new File(outputBaseFile.getAbsolutePath() + "." + CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics.getFileExtension());
        final File controlCodeMetricsFile  = new File(outputBaseFile.getAbsolutePath() + "." + CollectArraysVariantCallingMetrics.ArraysControlCodesSummaryMetrics.getFileExtension());
        summaryMetricsFile.deleteOnExit();
        detailMetricsFile.deleteOnExit();
        controlCodeMetricsFile.deleteOnExit();

        final CollectArraysVariantCallingMetrics collectArraysVariantCallingMetrics = new CollectArraysVariantCallingMetrics();
        collectArraysVariantCallingMetrics.INPUT = vcfFile;
        collectArraysVariantCallingMetrics.DBSNP = dbSnpFile;
        collectArraysVariantCallingMetrics.OUTPUT = outputBaseFile;
        collectArraysVariantCallingMetrics.NUM_PROCESSORS = 1;

        Assert.assertEquals(collectArraysVariantCallingMetrics.doWork(), 0);

        final MetricsFile<CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics, Comparable<?>> summary = new MetricsFile<>();
        summary.read(new FileReader(summaryMetricsFile));

        boolean parsedSummary = false;
        for (final CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics metrics : summary.getMetrics()) {
            Assert.assertEquals(metrics.NUM_ASSAYS, 638);
            Assert.assertEquals(metrics.NUM_NON_FILTERED_ASSAYS, 637);
            Assert.assertEquals(metrics.NUM_SNPS, 637);
            Assert.assertEquals(metrics.NUM_INDELS, 0);

            Assert.assertEquals(metrics.NUM_CALLS, expectedNumCalls);
            Assert.assertEquals(metrics.NUM_AUTOCALL_CALLS, 631);
            Assert.assertEquals(metrics.NUM_NO_CALLS, metrics.NUM_NON_FILTERED_ASSAYS - metrics.NUM_CALLS);

            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 65);
            Assert.assertEquals(metrics.NOVEL_SNPS, 572);
            Assert.assertEquals(metrics.NUM_FILTERED_ASSAYS, 1);
            Assert.assertEquals(metrics.NUM_ZEROED_OUT_ASSAYS, 1);

            Assert.assertEquals(metrics.PCT_DBSNP, 0.102041, 0.0001);
            Assert.assertEquals(metrics.CALL_RATE, (float) metrics.NUM_CALLS / metrics.NUM_NON_FILTERED_ASSAYS, 0.0001);
            Assert.assertEquals(metrics.AUTOCALL_CALL_RATE, (float) metrics.NUM_AUTOCALL_CALLS / metrics.NUM_NON_FILTERED_ASSAYS, 0.0001);

            Assert.assertEquals(metrics.NUM_SINGLETONS, 7);

            parsedSummary = true;
        }

        Assert.assertTrue(parsedSummary, "Did not parse summary metrics.");

        final MetricsFile<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics, Comparable<?>> detail = new MetricsFile<>();
        detail.read(new FileReader(detailMetricsFile));
        final List<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics> detailMetrics = detail.getMetrics();
        detail.getMetrics().stream().filter(metrics -> metrics.SAMPLE_ALIAS.equals("NA12878")).forEach(metrics -> {
            Assert.assertEquals(metrics.CHIP_WELL_BARCODE, "7991775143_R01C01");
            Assert.assertEquals(metrics.SAMPLE_ALIAS, "NA12878");
            if (expectedAnalysisVersion != null) {
                Assert.assertEquals((int) metrics.ANALYSIS_VERSION, 1);
            } else {
                Assert.assertNull(metrics.ANALYSIS_VERSION);
            }
            Assert.assertEquals(metrics.CHIP_TYPE, "HumanExome-12v1-1_A");
            Assert.assertTrue(metrics.AUTOCALL_PF);
            Assert.assertEquals(metrics.IMAGING_DATE, new Iso8601Date("2012-12-16T11:10:10-0500"));
            Assert.assertEquals(metrics.IS_ZCALLED, expectedIsZcalled);
            Assert.assertEquals(metrics.AUTOCALL_GENDER, expectedAutocallGender);
            Assert.assertEquals(metrics.FP_GENDER, "U");
            Assert.assertEquals(metrics.REPORTED_GENDER, expectedReportedGender);
            Assert.assertEquals(metrics.GENDER_CONCORDANCE_PF, expectedGenderConcordancePF);
            Assert.assertEquals(metrics.HET_PCT, 0.010989, 0.001);
            Assert.assertEquals(metrics.CLUSTER_FILE_NAME, "HumanExomev1_1_CEPH_A.egt");
            Assert.assertEquals((int) metrics.P95_GREEN, 4251);
            Assert.assertEquals((int) metrics.P95_RED, 5479);
            Assert.assertEquals(metrics.AUTOCALL_VERSION, expectedAutocallVersion);
            if (expectedZcallVersion != null) {
                Assert.assertEquals(metrics.ZCALL_VERSION, "1.0.0.0");
            } else {
                Assert.assertNull(metrics.ZCALL_VERSION);
            }
            Assert.assertEquals(metrics.EXTENDED_MANIFEST_VERSION, "1.3");
            Assert.assertEquals(metrics.HET_HOMVAR_RATIO, 0.777778, 0.001);
            Assert.assertEquals(metrics.SCANNER_NAME, "N370");
            Assert.assertEquals(metrics.NUM_ASSAYS, 638);
            Assert.assertEquals(metrics.NUM_NON_FILTERED_ASSAYS, 637);
            Assert.assertEquals(metrics.NUM_SNPS, 637);
            Assert.assertEquals(metrics.NUM_INDELS, 0);
            Assert.assertEquals(metrics.NUM_CALLS, expectedNumCalls);
            Assert.assertEquals(metrics.NUM_AUTOCALL_CALLS, 631);
            Assert.assertEquals(metrics.NUM_NO_CALLS, metrics.NUM_NON_FILTERED_ASSAYS - metrics.NUM_CALLS);

            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 65);
            Assert.assertEquals(metrics.NOVEL_SNPS, 572);
            Assert.assertEquals(metrics.NUM_FILTERED_ASSAYS, 1);
            Assert.assertEquals(metrics.NUM_ZEROED_OUT_ASSAYS, 1);
            Assert.assertEquals(metrics.PCT_DBSNP, 0.102041, 0.0001);
            Assert.assertEquals(metrics.CALL_RATE, (float) metrics.NUM_CALLS / metrics.NUM_NON_FILTERED_ASSAYS, 0.0001);
            Assert.assertEquals(metrics.AUTOCALL_CALL_RATE, (float) metrics.NUM_AUTOCALL_CALLS / metrics.NUM_NON_FILTERED_ASSAYS, 0.0001);

            Assert.assertEquals(metrics.NUM_SINGLETONS, 7);
            Assert.assertEquals(metrics.PIPELINE_VERSION, "foo");
        });

        Assert.assertEquals(detailMetrics.size(), 1, "Did not parse the desired number of detail metrics.");

        final MetricsFile<CollectArraysVariantCallingMetrics.ArraysControlCodesSummaryMetrics, Comparable<?>> controlCodes = new MetricsFile<>();
        controlCodes.read(new FileReader(controlCodeMetricsFile));
        final List<CollectArraysVariantCallingMetrics.ArraysControlCodesSummaryMetrics> controlCodeMetrics = controlCodes.getMetrics();
        controlCodes.getMetrics().stream().filter(metrics -> metrics.CONTROL.equals("DNP(High)")).forEach(metrics -> {
            Assert.assertEquals(metrics.CATEGORY, "Staining");
            Assert.assertEquals(metrics.RED, 6576);
            Assert.assertEquals(metrics.GREEN, 45);
        });

        Assert.assertEquals(controlCodeMetrics.size(), 23, "Did not parse the number of control code metrics.");
    }

    @Test
    public void testCollectArraysVariantCallingMetricsTwoSamples() throws IOException {
        final File dbSnpFile = new File(TEST_DBSNP_DATA_DIR, "mini.dbsnp.vcf");
        final File vcfFile = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "two_samples.vcf"), "two_samples");
        vcfFile.deleteOnExit();

        final File outputBaseFile = File.createTempFile("collectArraysVariantCallingMetrics", "");
        final File summaryMetricsFile = new File(outputBaseFile.getAbsolutePath() + "." + CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics.getFileExtension());
        final File detailMetricsFile  = new File(outputBaseFile.getAbsolutePath() + "." + CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics.getFileExtension());
        final File controlCodeMetricsFile  = new File(outputBaseFile.getAbsolutePath() + "." + CollectArraysVariantCallingMetrics.ArraysControlCodesSummaryMetrics.getFileExtension());
        summaryMetricsFile.deleteOnExit();
        detailMetricsFile.deleteOnExit();
        controlCodeMetricsFile.deleteOnExit();

        final CollectArraysVariantCallingMetrics collectArraysVariantCallingMetrics = new CollectArraysVariantCallingMetrics();
        collectArraysVariantCallingMetrics.INPUT = vcfFile;
        collectArraysVariantCallingMetrics.DBSNP = dbSnpFile;
        collectArraysVariantCallingMetrics.OUTPUT = outputBaseFile;
        collectArraysVariantCallingMetrics.NUM_PROCESSORS = 1;

        Assert.assertEquals(collectArraysVariantCallingMetrics.doWork(), 0);

        final MetricsFile<CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics, Comparable<?>> summary = new MetricsFile<>();
        summary.read(new FileReader(summaryMetricsFile));

        boolean parsedSummary = false;
        for (final CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics metrics : summary.getMetrics()) {
            Assert.assertEquals(metrics.NUM_ASSAYS, 1276);
            Assert.assertEquals(metrics.NUM_NON_FILTERED_ASSAYS, 1274);
            Assert.assertEquals(metrics.NUM_SNPS, 1274);
            Assert.assertEquals(metrics.NUM_INDELS, 0);

            Assert.assertEquals(metrics.NUM_CALLS, 1272);
            Assert.assertEquals(metrics.NUM_AUTOCALL_CALLS, 1262);
            Assert.assertEquals(metrics.NUM_NO_CALLS, metrics.NUM_NON_FILTERED_ASSAYS - metrics.NUM_CALLS);

            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 130);
            Assert.assertEquals(metrics.NOVEL_SNPS, 1144);
            Assert.assertEquals(metrics.NUM_FILTERED_ASSAYS, 2);
            Assert.assertEquals(metrics.NUM_ZEROED_OUT_ASSAYS, 2);

            Assert.assertEquals(metrics.PCT_DBSNP, 0.102041, 0.0001);
            Assert.assertEquals(metrics.CALL_RATE, (float) metrics.NUM_CALLS / metrics.NUM_NON_FILTERED_ASSAYS, 0.0001);
            Assert.assertEquals(metrics.AUTOCALL_CALL_RATE, (float) metrics.NUM_AUTOCALL_CALLS / metrics.NUM_NON_FILTERED_ASSAYS, 0.0001);

            Assert.assertEquals(metrics.NUM_SINGLETONS, 2);

            parsedSummary = true;
        }

        Assert.assertTrue(parsedSummary, "Did not parse summary metrics.");

        final MetricsFile<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics, Comparable<?>> detail = new MetricsFile<>();
        detail.read(new FileReader(detailMetricsFile));
        final List<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics> detailMetrics = detail.getMetrics();
        detail.getMetrics().stream().filter(metrics -> metrics.CHIP_WELL_BARCODE.equals("7991775143_R01C02")).forEach(metrics -> {
            Assert.assertEquals(metrics.CHIP_WELL_BARCODE, "7991775143_R01C02");
            Assert.assertEquals(metrics.SAMPLE_ALIAS, "NA12878");
            Assert.assertEquals(metrics.CHIP_TYPE, "HumanExome-12v1-1_A");
            Assert.assertTrue(metrics.AUTOCALL_PF);
            Assert.assertEquals(metrics.IMAGING_DATE, new Iso8601Date("2012-12-16T11:10:10-0500"));
            Assert.assertTrue(metrics.IS_ZCALLED);
            Assert.assertEquals(metrics.AUTOCALL_GENDER, "F");
            Assert.assertEquals(metrics.FP_GENDER, "N");
            Assert.assertEquals(metrics.REPORTED_GENDER, "N");
            Assert.assertFalse(metrics.GENDER_CONCORDANCE_PF);
            Assert.assertEquals(metrics.HET_PCT, 0.011006, 0.0001);
            Assert.assertEquals(metrics.CLUSTER_FILE_NAME, "HumanExomev1_1_CEPH_A.egt");
            Assert.assertEquals((int) metrics.P95_GREEN, 4251);
            Assert.assertEquals((int) metrics.P95_RED, 5479);
            Assert.assertEquals(metrics.AUTOCALL_VERSION, "3.0.0");
            Assert.assertEquals(metrics.ZCALL_VERSION, "1.0.0.0");
            Assert.assertEquals(metrics.EXTENDED_MANIFEST_VERSION, "1.3");
            Assert.assertEquals(metrics.HET_HOMVAR_RATIO, 0.875, 0.0001);
            Assert.assertEquals(metrics.SCANNER_NAME, "N370");
            Assert.assertEquals(metrics.NUM_ASSAYS, 638);
            Assert.assertEquals(metrics.NUM_NON_FILTERED_ASSAYS, 637);
            Assert.assertEquals(metrics.NUM_SNPS, 637);
            Assert.assertEquals(metrics.NUM_INDELS, 0);
            Assert.assertEquals(metrics.NUM_CALLS, 636);
            Assert.assertEquals(metrics.NUM_AUTOCALL_CALLS, 631);
            Assert.assertEquals(metrics.NUM_NO_CALLS, metrics.NUM_NON_FILTERED_ASSAYS - metrics.NUM_CALLS);

            Assert.assertEquals(metrics.NUM_IN_DB_SNP, 65);
            Assert.assertEquals(metrics.NOVEL_SNPS, 572);
            Assert.assertEquals(metrics.NUM_FILTERED_ASSAYS, 1);
            Assert.assertEquals(metrics.NUM_ZEROED_OUT_ASSAYS, 1);
            Assert.assertEquals(metrics.PCT_DBSNP, 0.102041, 0.0001);
            Assert.assertEquals(metrics.CALL_RATE, (float) metrics.NUM_CALLS / metrics.NUM_NON_FILTERED_ASSAYS, 0.0001);
            Assert.assertEquals(metrics.AUTOCALL_CALL_RATE, (float) metrics.NUM_AUTOCALL_CALLS / metrics.NUM_NON_FILTERED_ASSAYS, 0.0001);

            Assert.assertEquals(metrics.NUM_SINGLETONS, 1);
            Assert.assertNull(metrics.PIPELINE_VERSION);
        });

        Assert.assertEquals(detailMetrics.size(), 2, "Did not parse the desired number of detail metrics.");

        final MetricsFile<CollectArraysVariantCallingMetrics.ArraysControlCodesSummaryMetrics, Comparable<?>> controlCodes = new MetricsFile<>();
        controlCodes.read(new FileReader(controlCodeMetricsFile));
        final List<CollectArraysVariantCallingMetrics.ArraysControlCodesSummaryMetrics> controlCodeMetrics = controlCodes.getMetrics();
        controlCodes.getMetrics().stream().filter(metrics -> metrics.CONTROL.equals("DNP(High)")).forEach(metrics -> {
            Assert.assertEquals(metrics.CATEGORY, "Staining");
            Assert.assertEquals(metrics.RED, 6576);
            Assert.assertEquals(metrics.GREEN, 45);
        });

        Assert.assertEquals(controlCodeMetrics.size(), 23, "Did not parse the number of control code metrics.");
    }

    @DataProvider(name = "genderConcordanceTypes")
    public Object[][] testSexConcordanceTypesProvider() {
        return new Object[][] {
                {Sex.Male, Sex.Male, Sex.Male, true},
                {Sex.Male, Sex.Male, Sex.Female, false},
                {Sex.Male, Sex.Male, Sex.Unknown, true},
                {Sex.Male, Sex.Male, Sex.NotReported, true},
                {Sex.Male, Sex.Female, Sex.Male, false},
                {Sex.Male, Sex.Female, Sex.Female, false},
                {Sex.Male, Sex.Female, Sex.Unknown, false},
                {Sex.Male, Sex.Female, Sex.NotReported, false},
                {Sex.Male, Sex.Unknown, Sex.Male, true},
                {Sex.Male, Sex.Unknown, Sex.Female, false},
                {Sex.Male, Sex.Unknown, Sex.Unknown, false},
                {Sex.Male, Sex.Unknown, Sex.NotReported, false},
                {Sex.Male, Sex.NotReported, Sex.Male, true},
                {Sex.Male, Sex.NotReported, Sex.Female, false},
                {Sex.Male, Sex.NotReported, Sex.Unknown, false},
                {Sex.Male, Sex.NotReported, Sex.NotReported, false},

                {Sex.Female, Sex.Male, Sex.Male, false},
                {Sex.Female, Sex.Male, Sex.Female, false},
                {Sex.Female, Sex.Male, Sex.Unknown, false},
                {Sex.Female, Sex.Male, Sex.NotReported, false},
                {Sex.Female, Sex.Female, Sex.Male, false},
                {Sex.Female, Sex.Female, Sex.Female, true},
                {Sex.Female, Sex.Female, Sex.Unknown, true},
                {Sex.Female, Sex.Female, Sex.NotReported, true},
                {Sex.Female, Sex.Unknown, Sex.Male, false},
                {Sex.Female, Sex.Unknown, Sex.Female, true},
                {Sex.Female, Sex.Unknown, Sex.Unknown, false},
                {Sex.Female, Sex.Unknown, Sex.NotReported, false},
                {Sex.Female, Sex.NotReported, Sex.Male, false},
                {Sex.Female, Sex.NotReported, Sex.Female, true},
                {Sex.Female, Sex.NotReported, Sex.Unknown, false},
                {Sex.Female, Sex.NotReported, Sex.NotReported, false},

                {Sex.Unknown, Sex.Male, Sex.Male, true},
                {Sex.Unknown, Sex.Male, Sex.Female, false},
                {Sex.Unknown, Sex.Male, Sex.Unknown, false},
                {Sex.Unknown, Sex.Male, Sex.NotReported, false},
                {Sex.Unknown, Sex.Female, Sex.Male, false},
                {Sex.Unknown, Sex.Female, Sex.Female, true},
                {Sex.Unknown, Sex.Female, Sex.Unknown, false},
                {Sex.Unknown, Sex.Female, Sex.NotReported, false},
                {Sex.Unknown, Sex.Unknown, Sex.Male, false},
                {Sex.Unknown, Sex.Unknown, Sex.Female, false},
                {Sex.Unknown, Sex.Unknown, Sex.Unknown, false},
                {Sex.Unknown, Sex.Unknown, Sex.NotReported, false},
                {Sex.Unknown, Sex.NotReported, Sex.Male, false},
                {Sex.Unknown, Sex.NotReported, Sex.Female, false},
                {Sex.Unknown, Sex.NotReported, Sex.Unknown, false},
                {Sex.Unknown, Sex.NotReported, Sex.NotReported, false},

                {Sex.NotReported, Sex.Male, Sex.Male, true},
                {Sex.NotReported, Sex.Male, Sex.Female, false},
                {Sex.NotReported, Sex.Male, Sex.Unknown, false},
                {Sex.NotReported, Sex.Male, Sex.NotReported, false},
                {Sex.NotReported, Sex.Female, Sex.Male, false},
                {Sex.NotReported, Sex.Female, Sex.Female, true},
                {Sex.NotReported, Sex.Female, Sex.Unknown, false},
                {Sex.NotReported, Sex.Female, Sex.NotReported, false},
                {Sex.NotReported, Sex.Unknown, Sex.Male, false},
                {Sex.NotReported, Sex.Unknown, Sex.Female, false},
                {Sex.NotReported, Sex.Unknown, Sex.Unknown, false},
                {Sex.NotReported, Sex.Unknown, Sex.NotReported, false},
                {Sex.NotReported, Sex.NotReported, Sex.Male, false},
                {Sex.NotReported, Sex.NotReported, Sex.Female, false},
                {Sex.NotReported, Sex.NotReported, Sex.Unknown, false},
                {Sex.NotReported, Sex.NotReported, Sex.NotReported, false}
        };
    }

    @Test(dataProvider="genderConcordanceTypes")
    public void testSexConcordance(Sex reportedSex, Sex fingerprintSex, Sex autocallSex, boolean expectedToPass) {
        boolean pass = CollectArraysVariantCallingMetrics.getSexConcordance(reportedSex.toSymbol(), fingerprintSex.toSymbol(), autocallSex.toSymbol());
        Assert.assertEquals(pass, expectedToPass);
    }
}
