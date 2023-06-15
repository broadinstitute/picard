/*
 * The MIT License
 *
 * Copyright (c) 2020 The Broad Institute
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
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Writer;
import java.lang.reflect.Field;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Compare two metrics files.
 */
@CommandLineProgramProperties(
        summary = CompareMetrics.USAGE_SUMMARY + CompareMetrics.USAGE_DETAIL,
        oneLineSummary = CompareMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CompareMetrics extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Compare two metrics files.";
    static final String USAGE_DETAIL = "This tool compares the metrics and histograms generated from metric tools to determine " +
            "if the generated results are identical.  Note that if there are differences in metric values, this tool describes those differences " +
            "as the change of the second input metric relative to the first. " +
            "<br /><br />  " +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CompareMetrics \\<br />" +
            "      INPUT=metricfile1.txt \\<br />" +
            "      INPUT=metricfile2.txt \\<br />" +
            "      METRICS_TO_IGNORE=INSERT_LENGTH \\<br />" +
            "      METRIC_ALLOWABLE_RELATIVE_CHANGE=HET_HOM_RATIO:0.0005 \\<br />" +
            "      IGNORE_HISTOGRAM_DIFFERENCES=false" +
            "</pre>" +
            "<hr />";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Metric files to compare.",
            minElements = 2,
            maxElements = 2)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file to write comparison results to.", optional = true)
    public File OUTPUT;

    @Argument(
            doc = "Output file to write table of differences to.",
            optional = true
    )
    public File OUTPUT_TABLE;

    @Argument(shortName = "MI",
            doc = "Metrics to ignore. Any metrics specified here will be excluded from comparison by the tool.  " +
                    "Note that while the values of these metrics are not compared, if they are missing from either file that will be considered a difference.  " +
                    "Use METRICS_NOT_REQUIRED to specify metrics which can be missing from either file without being considered a difference.",
            optional = true)
    public List<String> METRICS_TO_IGNORE;

    @Argument(shortName = "MNR",
            doc = "Metrics which are not required.  Any metrics specified here may be missing from either of the files in the comparison, and this will not affect " +
                    "the result of the comparison.  If metrics specified here are included in both files, their results will not be compared (they will be treated as " +
                    "METRICS_TO_IGNORE.",
            optional = true)
    public List<String> METRICS_NOT_REQUIRED;

    @Argument(shortName = "MARC",
            doc = "Metric Allowable Relative Change. A colon separate pair of metric name and an absolute relative change.  For any metric specified here, " +
                    " when the values are compared between the two files, the program will allow that much relative change between the " +
                    " two values.",
            optional = true)
    public List<String> METRIC_ALLOWABLE_RELATIVE_CHANGE;

    @Argument(shortName = "IHD",
            doc = "Ignore any differences between the two metric file's histograms (useful if using the 'METRIC_ALLOWABLE_RELATIVE_CHANGE')",
            optional = true)
    public boolean IGNORE_HISTOGRAM_DIFFERENCES = false;

    @Argument(
            doc = "Columns to use as keys for matching metrics rows that should agree.  If not specified, it is assumed that rows should be in the same order.",
            optional = true
    )
    public List<String> KEY;

    private final List<String> differences = new ArrayList<>();
    private final List<MetricComparisonDifferences> valueDifferences = new ArrayList<>();

    private String metricClassName = "Unknown";

    private static final Log log = Log.getInstance(CompareMetrics.class);

    protected final Map<String, Double> MetricToAllowableRelativeChange = new HashMap<>();

    @Override
    protected int doWork() {
        IOUtil.assertFilesAreReadable(INPUT);
        try {
            final int retVal = compareMetricsFiles(INPUT.get(0), INPUT.get(1));
            final String status = retVal == 0 ? "equal" : "NOT equal";
            log.info(metricClassName + " Metric files " + INPUT.get(0) + " and " + INPUT.get(1) + " are " + status);
            if (!differences.isEmpty()) {
                for (String difference : differences) {
                    log.error(difference);
                }
            }
            if (OUTPUT != null) {
                final String header = "Comparison of " + metricClassName + " metrics between files "
                        + INPUT.get(0).getAbsolutePath() + " and " + INPUT.get(1).getAbsolutePath() + "\n\nMetrics are " + status;
                writeTextToFile(OUTPUT, header, differences);
            }
            if (OUTPUT_TABLE != null) {
                final MetricsFile<MetricComparisonDifferences, ?> metricDifferencesFile = getMetricsFile();
                metricDifferencesFile.addAllMetrics(valueDifferences);
                metricDifferencesFile.write(OUTPUT_TABLE);
            }
            return retVal;
        } catch (final Exception e) {
            throw new PicardException(e.getMessage(), e);
        }
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errs = new ArrayList<>();

        if (OUTPUT != null) {
            IOUtil.assertFileIsWritable(OUTPUT);
        }

        if (METRIC_ALLOWABLE_RELATIVE_CHANGE != null) {
            for (String diffs : METRIC_ALLOWABLE_RELATIVE_CHANGE) {
                String[] pair = diffs.split(":");
                if (pair.length == 2) {
                    String name = pair[0];
                    try {
                        double value = Double.parseDouble(pair[1]);
                        if (value > 0) {
                            MetricToAllowableRelativeChange.put(name, value);
                        } else {
                            errs.add("Value for numeric component of Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE' must be > 0.0");
                        }
                    } catch (NumberFormatException ne) {
                        errs.add("Invalid value for numeric component of Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE'");
                    }
                } else {
                    errs.add("Invalid value for Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE'");
                }
            }
        }

        if (errs.isEmpty()) {
            return null;
        }
        return errs.toArray(new String[0]);
    }

    private int compareMetricsFiles(final File metricFile1, final File metricFile2) throws IOException, IllegalAccessException, NoSuchFieldException {
        final MetricsFile<?, ?> mf1 = new MetricsFile<>();
        final MetricsFile<?, ?> mf2 = new MetricsFile<>();
        mf1.read(new FileReader(metricFile1));
        mf2.read(new FileReader(metricFile2));
        if (!mf1.getMetrics().isEmpty()) {
            metricClassName = mf1.getMetrics().get(0).getClass().getName();
        }
        else if (!mf2.getMetrics().isEmpty()) {
            metricClassName = mf2.getMetrics().get(0).getClass().getName();
        }
        else {
            metricClassName = "Unknown";
        }
        final boolean histogramsEqual = mf1.areHistogramsEqual(mf2);
        if (mf1.areMetricsEqual(mf2)) {
            if (histogramsEqual) {
                return 0;
            }
            else {
                if (IGNORE_HISTOGRAM_DIFFERENCES) {
                    differences.add("Metrics Histograms differ, but the 'IGNORE_HISTOGRAM_DIFFERENCES' flag is set.");
                    return 0;
                } else {
                    differences.add("Metrics Histograms differ");
                    return 1;
                }
            }
        }
        if (mf1.getMetrics().size() != mf2.getMetrics().size()) {
            differences.add("Number of metric rows differ between " + metricFile1.getAbsolutePath() + " and " + metricFile2.getAbsolutePath());
            return 1;
        }
        if (!mf1.getMetrics().get(0).getClass().equals(mf2.getMetrics().get(0).getClass())) {
            throw new PicardException("Metrics are of differing class between " + metricFile1.getAbsolutePath() + " and " + metricFile2.getAbsolutePath());
        }

        final Set<String> columns1 = new HashSet<>(mf1.getMetricsColumnLabels());
        final Set<String> columns2 = new HashSet<>(mf2.getMetricsColumnLabels());
        columns1.removeAll(METRICS_NOT_REQUIRED);
        columns2.removeAll(METRICS_NOT_REQUIRED);
        if (!columns1.equals(columns2)) {
            final Set<String> missingColumns = new HashSet<>(columns1);
            missingColumns.removeAll(columns2);
            final Set<String> inColumns2NotColumns1 = new HashSet<>(columns2);
            inColumns2NotColumns1.removeAll(columns1);
            missingColumns.addAll(inColumns2NotColumns1);
            differences.add("Metric columns differ between " + metricFile1.getAbsolutePath() + " and " + metricFile2.getAbsolutePath() + " (" +
                    StringUtil.join(",", missingColumns) + ")");

            return 1;
        }

        validateMetricNames(mf1, metricFile1, METRICS_TO_IGNORE);
        validateMetricNames(mf1, metricFile1, MetricToAllowableRelativeChange.keySet());

        Set<String> metricsToIgnore = new HashSet<>(METRICS_TO_IGNORE);
        metricsToIgnore.addAll(METRICS_NOT_REQUIRED);

        final Class<? extends MetricBase> metricClass = mf1.getMetrics().get(0).getClass();

        final Field[] fields = metricClass.getFields();

        int retVal = 0;
        if (KEY.size() == 0) {
            //key by row number
            int rowNumber = -1;
            Iterator<?> mf1Iterator = mf1.getMetrics().iterator();
            Iterator<?> mf2Iterator = mf2.getMetrics().iterator();
            while (mf1Iterator.hasNext()) {
                rowNumber++;
                MetricBase metric1 = (MetricBase) mf1Iterator.next();
                MetricBase metric2 = (MetricBase) mf2Iterator.next();
                if (compareMetricsForEntry(metric1, metric2, fields, metricsToIgnore, String.valueOf(rowNumber)) == 1) {
                    retVal = 1;
                }
            }
        } else {
            //build each map of metrics
            final Map<List<Object>, ? extends MetricBase> metricMap1 = buildMetricsMap(mf1.getMetrics());
            final Map<List<Object>, ? extends MetricBase> metricMap2 = buildMetricsMap(mf2.getMetrics());

            for (final Map.Entry<List<Object>, ? extends MetricBase> entry1 : metricMap1.entrySet()) {
                final List<Object> key = entry1.getKey();
                final MetricBase metric1 = entry1.getValue();

                final MetricBase metric2 = metricMap2.remove(key);
                if (metric2 != null) {
                    if (compareMetricsForEntry(metric1, metric2, fields, metricsToIgnore, StringUtil.join(",", key)) == 1) {
                        retVal = 1;
                    }
                } else {
                    differences.add("KEY " + StringUtil.join(",", key) + " found in " + metricFile1 + " but not in " + metricFile2);
                    retVal = 1;
                }
            }
            //check that all entries in metricMap2 have been matched
            for (final Map.Entry<List<Object>, ? extends MetricBase> entry : metricMap2.entrySet()) {
                differences.add("KEY " + StringUtil.join(",", entry.getKey()) + " found in " + metricFile2 + " but not in " + metricFile1);
                retVal = 1;
            }
        }

        if (!IGNORE_HISTOGRAM_DIFFERENCES) {
            if (!histogramsEqual) {
                final String diffString = "Metric Histograms differ";
                differences.add(diffString);
            }
            if (retVal == 0 && !histogramsEqual) {
                retVal = 1;
            }
        }

        return retVal;
    }

    protected Map<List<Object>, MetricBase> buildMetricsMap(final List<? extends MetricBase> metrics) throws NoSuchFieldException, IllegalAccessException {
        final HashMap<List<Object>, MetricBase> retMap = new LinkedHashMap<>();
        final Class<? extends MetricBase> clazz = metrics.get(0).getClass();
        final List<Field> keyFields = new ArrayList<>();
        for (final String key : KEY) {
            final Field keyField = clazz.getField(key);
            keyFields.add(keyField);
        }

        for (final MetricBase metric : metrics) {
            final List<Object> mapKey = new ArrayList<>();
            for (final Field keyField : keyFields) {
                mapKey.add(keyField.get(metric));
            }
            retMap.put(mapKey, metric);
        }

        return retMap;
    }

    protected int compareMetricsForEntry(final MetricBase metric1, final MetricBase metric2, final Field[] fields, final Set<String> metricsToIgnore, final String key) throws IllegalAccessException {
        int retVal = 0;
        for (Field field : fields) {
            if (!metricsToIgnore.contains(field.getName())) {
                final Object value1 = field.get(metric1);
                final Object value2 = field.get(metric2);
                SimpleResult result = compareMetricValues(value1, value2, field.getName());
                if (!result.equal) {
                    retVal = 1;
                    final String diffString = "Key: " + key + " Metric: " + field.getName() +
                            " values differ. Value1: " + value1 + " Value2: " + value2 + " " + result.description;
                    differences.add(diffString);
                    final MetricComparisonDifferences metricDifferences = new MetricComparisonDifferences();
                    metricDifferences.KEY = key;
                    metricDifferences.METRIC = field.getName();
                    metricDifferences.VALUE1 = value1;
                    metricDifferences.VALUE2 = value2;
                    valueDifferences.add(metricDifferences);
                }
            }
        }
        return retVal;
    }

    protected SimpleResult compareMetricValues(final Object value1, final Object value2, final String metricName) {
        boolean equal = true;
        String description = "";
        if (value1 == null || value2 == null) {
            if (value1 != null || value2 != null) {
                equal = false;
                description = "One of the values is null";
            }
        } else {
            if (value1 instanceof Number) {
                double numValue1 = ((Number) value1).doubleValue();
                double numValue2 = ((Number) value2).doubleValue();
                double absoluteChange = 0;
                if (!Double.isNaN(numValue1) || !Double.isNaN(numValue2)) {
                    absoluteChange = numValue2 - numValue1;
                }
                if (absoluteChange != 0) {
                    double relativeChange = numValue1 == 0 ? Double.MAX_VALUE : absoluteChange / numValue1;
                    if (MetricToAllowableRelativeChange.containsKey(metricName)) {
                        double allowableRelativeChange = MetricToAllowableRelativeChange.get(metricName);
                        if (Math.abs(relativeChange) >= allowableRelativeChange) {
                            equal = false;
                            description = "Changed by " + absoluteChange + " (relative change of " + relativeChange +
                                    ") which is outside of the allowable relative change tolerance of " + allowableRelativeChange;
                        } else {
                            equal = true;
                            description = "Changed by " + absoluteChange + " (relative change of " + relativeChange +
                                    ") which is within the allowable relative change tolerance of " + allowableRelativeChange;
                        }
                    } else {
                        equal = false;
                        description = "Changed by " + absoluteChange + " (relative change of " + relativeChange + ")";
                    }
                }
            } else {
                if (!value1.equals(value2)) {
                    equal = false;
                    description = "";
                }
            }
        }
        return new SimpleResult(equal, description);
    }

    private static void writeTextToFile(final File output, final String header, final List<String> textLines) throws IOException {
        try (Writer writer = Files.newBufferedWriter(output.toPath())) {
            writer.write(header + "\n\n");
            writer.write(String.join("\n", textLines));
        }
    }

    private static void validateMetricNames(final MetricsFile<?, ?> metrics, final File metricsFile, final Collection<String> metricNames) {
        Set<String> metricsToIgnore = new HashSet<>(metricNames);
        metricsToIgnore.removeAll(metrics.getMetricsColumnLabels());
        if (!metricsToIgnore.isEmpty()) {
            throw new PicardException("Metric(s) of the name: " + String.join(", ", metricsToIgnore) + " were not found in " + metricsFile.getAbsolutePath());
        }
    }

    static class SimpleResult {
        final boolean equal;
        final String description;

        public SimpleResult(boolean equal, String description) {
            this.equal = equal;
            this.description = description;
        }
    }

    static public class MetricComparisonDifferences extends MetricBase {
        public String KEY;
        public String METRIC;
        public Object VALUE1;
        public Object VALUE2;
    }
}
