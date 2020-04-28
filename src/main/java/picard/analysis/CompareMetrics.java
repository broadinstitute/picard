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

    @Argument(shortName = "MI",
            doc = "Metrics to ignore. Any metrics specified here will be excluded from comparison by the tool.",
            optional = true)
    public List<String> METRICS_TO_IGNORE;

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

    private final List<String> differences = new ArrayList<>();

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
            return retVal;
        } catch (final Exception e) {
            throw new PicardException(e.getMessage());
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

    private int compareMetricsFiles(final File metricFile1, final File metricFile2) throws IOException, IllegalAccessException {
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
        else if (!mf1.getMetrics().get(0).getClass().equals(mf2.getMetrics().get(0).getClass())) {
            throw new PicardException("Metrics are of differing class between " + metricFile1.getAbsolutePath() + " and " + metricFile2.getAbsolutePath());
        }
        else if (!mf1.getMetricsColumnLabels().equals(mf2.getMetricsColumnLabels())) {
            differences.add("Metric columns differ between " + metricFile1.getAbsolutePath() + " and " + metricFile2.getAbsolutePath());
            return 1;
        }

        validateMetricNames(mf1, metricFile1, METRICS_TO_IGNORE);
        validateMetricNames(mf1, metricFile1, MetricToAllowableRelativeChange.keySet());

        Set<String> metricsToIgnore = new HashSet<>(METRICS_TO_IGNORE);

        final Class<? extends MetricBase> metricClass = mf1.getMetrics().get(0).getClass();

        int retVal = 0;
        int rowNumber = -1;
        final Field[] fields = metricClass.getFields();
        Iterator<?> mf1Iterator = mf1.getMetrics().iterator();
        Iterator<?> mf2Iterator = mf2.getMetrics().iterator();
        while (mf1Iterator.hasNext()) {
            rowNumber++;
            MetricBase metric1 = (MetricBase) mf1Iterator.next();
            MetricBase metric2 = (MetricBase) mf2Iterator.next();
            for (Field field : fields) {
                if (!metricsToIgnore.contains(field.getName())) {
                    final Object value1 = field.get(metric1);
                    final Object value2 = field.get(metric2);
                    SimpleResult result = compareMetricValues(value1, value2, field.getName());
                    if (!result.equal) {
                        retVal = 1;
                        final String diffString = "Row: " + rowNumber + " Metric: " + field.getName() +
                                " values differ. Value1: " + value1 + " Value2: " + value2 + " " + result.description;
                        differences.add(diffString);
                    }
                }
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
}
