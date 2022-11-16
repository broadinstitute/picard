package picard.analysis;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.util.Collections;
import java.util.List;

import static java.util.Arrays.asList;

public class CompareMetricsTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/analysis/metrics");

    @DataProvider(name = "testCompareMetricsDataProvider")
    public Object[][] testCompareMetricsDataProvider() {
        return new Object[][]{
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", null, null, null, null, null, 1 },
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", "GENOME_TERRITORY", Collections.singletonList("SD_COVERAGE:0.00002"), true, null, null, 0 },
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", "GENOME_TERRITORY", Collections.singletonList("SD_COVERAGE:0.00002"), false, null, null, 1 },
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", "GENOME_TERRITORY", Collections.singletonList("SD_COVERAGE:0.00002"), null, null, null, 1 },
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", "GENOME_TERRITORY", Collections.singletonList("SD_COVERAGE:0.00001"), null, null, null, 1 },
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", null, asList("GENOME_TERRITORY:0.0001", "SD_COVERAGE:0.00002"), true, null, null, 0 },
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", null, asList("GENOME_TERRITORY:0.0001", "SD_COVERAGE:0.00002"), false, null, null, 1 },
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", null, asList("GENOME_TERRITORY:0.0001", "SD_COVERAGE:0.00002"), null, null, null, 1 },
                { "wgs_test1.wgs_metrics", "wgs_test2.wgs_metrics", null, asList("GENOME_TERRITORY:0.0001", "SD_COVERAGE:0.00001"), null, null, null, 1 },
                { "wgs_test1.wgs_metrics", "wgs_test1.wgs_metrics", null, null, null, null, null, 0 },
                { "wgs_test1.raw_wgs_metrics", "wgs_test2.raw_wgs_metrics", null, null, null, null, null, 1 },
                { "wgs_test1.2rows.raw_wgs_metrics", "wgs_test1.2rows.raw_wgs_metrics", null, null, null, null, null, 0 },
                { "wgs_test1.2rows.raw_wgs_metrics", "wgs_test2.2rows.raw_wgs_metrics", null, null, null, null, null, 1 },
                { "wgs_test1.raw_wgs_metrics", "wgs_test2.2rows.raw_wgs_metrics", null, null, null, null, null, 1 },
                { "wgs_test1.fingerprinting_summary_metrics", "wgs_test2.fingerprinting_summary_metrics", null, null, null, null, null, 1 },
                { "test1.arrays_variant_calling_detail_metrics", "test2.arrays_variant_calling_detail_metrics", null, null, null, null, null, 1 },
                { "test1.arrays_variant_calling_detail_metrics", "test2.arrays_variant_calling_detail_metrics", "AUTOCALL_DATE", null, null, null, null, 0 },
                { "test1.arrays_variant_calling_summary_metrics", "test2.arrays_variant_calling_summary_metrics", null, null, null, null, null, 0 },
                { "wgs_test1.2rows.fingerprinting_summary_metrics", "wgs_test2.2rows.fingerprinting_summary_metrics", null, null, false, null, null, 1},
                { "wgs_test1.2rows.fingerprinting_summary_metrics", "wgs_test2.2rows.fingerprinting_summary_metrics", null, null, false, Collections.singletonList("SAMPLE"), null, 0},
                { "wgs_test1.missing_metric.wgs_metrics", "wgs_test1.wgs_metrics", null, null, false, null, null, 1},
                { "wgs_test1.missing_metric.wgs_metrics", "wgs_test1.wgs_metrics", null, null, false, null, Collections.singletonList("HET_SNP_SENSITIVITY"), 0}
        };
    }

    @Test(dataProvider = "testCompareMetricsDataProvider")
    public void testCompareMetrics(final String file1, final String file2,
                                   final String metricsToIgnore, final List<String> metricsToAllowableRelativeChange,
                                   final Boolean ignoreHistogramDifferences, final List<String> keys, final List<String> metricsNotRequired, final int expectedReturnValue) {
        final File input1 = new File(TEST_DATA_DIR, file1);
        final File input2 = new File(TEST_DATA_DIR, file2);
        final CompareMetrics compareMetrics = new CompareMetrics();
        compareMetrics.INPUT = asList(input1, input2);
        if (metricsToIgnore != null) {
            compareMetrics.METRICS_TO_IGNORE = Collections.singletonList(metricsToIgnore);
        }
        if (metricsToAllowableRelativeChange != null) {
            compareMetrics.METRIC_ALLOWABLE_RELATIVE_CHANGE = metricsToAllowableRelativeChange;
        }
        if (ignoreHistogramDifferences != null) {
            compareMetrics.IGNORE_HISTOGRAM_DIFFERENCES = ignoreHistogramDifferences;
        }
        if (keys != null) {
            compareMetrics.KEY = keys;
        }
        if (metricsNotRequired != null) {
            compareMetrics.METRICS_NOT_REQUIRED = metricsNotRequired;
        }
        compareMetrics.customCommandLineValidation();
        Assert.assertEquals(compareMetrics.instanceMain(new String[0]), expectedReturnValue);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFailCompareMetricsOnDifferentClasses() {
        final File input1 = new File(TEST_DATA_DIR, "wgs_test1.wgs_metrics");
        final File input2 = new File(TEST_DATA_DIR, "wgs_test1.raw_wgs_metrics");
        final CompareMetrics compareMetrics = new CompareMetrics();
        compareMetrics.INPUT = asList(input1, input2);
        Assert.assertEquals(compareMetrics.instanceMain(new String[0]), 1);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFailCompareMetricsInvalidMetricToIgnore() {
        final File input1 = new File(TEST_DATA_DIR, "test1.arrays_variant_calling_detail_metrics");
        final File input2 = new File(TEST_DATA_DIR, "test2.arrays_variant_calling_detail_metrics");
        final CompareMetrics compareMetrics = new CompareMetrics();
        compareMetrics.INPUT = asList(input1, input2);
        compareMetrics.METRICS_TO_IGNORE = Collections.singletonList("NONEXISTENT_METRIC");
        Assert.assertEquals(compareMetrics.instanceMain(new String[0]), 1);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFailCompareMetricsInvalidMetricToAllowableRelativeChange() {
        final File input1 = new File(TEST_DATA_DIR, "test1.arrays_variant_calling_detail_metrics");
        final File input2 = new File(TEST_DATA_DIR, "test2.arrays_variant_calling_detail_metrics");
        final CompareMetrics compareMetrics = new CompareMetrics();
        compareMetrics.INPUT = asList(input1, input2);
        compareMetrics.METRIC_ALLOWABLE_RELATIVE_CHANGE = Collections.singletonList("NONEXISTENT_METRIC:0.3");
        Assert.assertEquals(compareMetrics.instanceMain(new String[0]), 1);
    }

    @DataProvider(name = "testCompareMetricValuesDataProvider")
    public Object[][] testCompareMetricValuesDataProvider() {
        return new Object[][]{
                { 1.0, 1.0, true, "", null },
                { 1.0, 1.1, false, "Changed by 0.10000000000000009 (relative change of 0.10000000000000009) which is outside of the allowable relative change tolerance of 0.05", 0.05 },
                { 1.0, 1.1, false, "Changed by 0.10000000000000009 (relative change of 0.10000000000000009) which is outside of the allowable relative change tolerance of 0.1", 0.1 },
                { 1.0, 1.1, true, "Changed by 0.10000000000000009 (relative change of 0.10000000000000009) which is within the allowable relative change tolerance of 0.1001", 0.1001 },
                { 1.0, 1.1, true, "Changed by 0.10000000000000009 (relative change of 0.10000000000000009) which is within the allowable relative change tolerance of 0.2", 0.2 },
                { 1.0f, 1.0f, true, "", null },
                { 1.0f, 1.1f, false, "Changed by 0.10000002384185791 (relative change of 0.10000002384185791) which is outside of the allowable relative change tolerance of 0.05", 0.05 },
                { 1.0f, 1.1f, false, "Changed by 0.10000002384185791 (relative change of 0.10000002384185791) which is outside of the allowable relative change tolerance of 0.1", 0.1 },
                { 1.0f, 1.1f, true, "Changed by 0.10000002384185791 (relative change of 0.10000002384185791) which is within the allowable relative change tolerance of 0.1001", 0.1001 },
                { 1.0f, 1.1f, true, "Changed by 0.10000002384185791 (relative change of 0.10000002384185791) which is within the allowable relative change tolerance of 0.2", 0.2 },
                { 1, 1, true, "", null },
                { 1, 2, false, "Changed by 1.0 (relative change of 1.0)", null },
                { 1, 0, false, "Changed by -1.0 (relative change of -1.0)", null },
                { 1f, 1f, true, "", null },
                { Float.NaN, Float.NaN, true, "", null },
                { Double.NaN, Double.NaN, true, "", null }
        };
    }

    @Test(dataProvider = "testCompareMetricValuesDataProvider")
    public void testCompareMetrics(final Object value1, final Object value2,
                                   final boolean expectedEqual, final String expectedDescription, final Double marcValue) {
        final CompareMetrics compareMetrics = new CompareMetrics();
        compareMetrics.METRIC_ALLOWABLE_RELATIVE_CHANGE = Collections.singletonList("M:" + marcValue);
        compareMetrics.customCommandLineValidation();
        final CompareMetrics.SimpleResult result = compareMetrics.compareMetricValues(value1, value2, "M");
        Assert.assertEquals(result.equal, expectedEqual);
        Assert.assertEquals(result.description, expectedDescription);
    }
}
