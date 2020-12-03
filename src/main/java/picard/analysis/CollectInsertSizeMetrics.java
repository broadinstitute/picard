/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.analysis.directed.InsertSizeMetricsCollector;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.util.RExecutor;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Command line program to read non-duplicate insert sizes, create a Histogram
 * and report distribution statistics.
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
@CommandLineProgramProperties(
        summary = CollectInsertSizeMetrics.USAGE_SUMMARY + CollectInsertSizeMetrics.USAGE_DETAILED,
        oneLineSummary = CollectInsertSizeMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CollectInsertSizeMetrics extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "Collect metrics about the insert size distribution of a paired-end library. ";
    static final String USAGE_DETAILED = "This tool provides useful metrics for validating library construction including " +
            "the insert size distribution and read orientation of paired-end libraries.</p>" +
            "" +
            "The expected proportions of these metrics vary depending on the type of library preparation used, resulting from " +
            "technical differences between pair-end libraries and mate-pair libraries. For a brief primer on paired-end sequencing " +
            "and mate-pair reads, see the " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=6327'>GATK Dictionary</a>." +
            "" +
            "<p>The CollectInsertSizeMetrics tool outputs the percentages of read pairs in each of the three orientations " +
            "(FR, RF, and TANDEM) as a histogram. In addition, the insert size distribution is output as both a histogram " +
            "(.insert_size_Histogram.pdf) and as a data table (.insert_size_metrics.txt).</p>" +
            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectInsertSizeMetrics \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=insert_size_metrics.txt \\<br />" +
            "      H=insert_size_histogram.pdf \\<br />" +
            "      M=0.5" +
            "</pre>"    +
            "Note: If processing a small file, set the minimum percentage option (M) to 0.5, otherwise an error may occur. "+
            "<br /><br />" +
            "Please see <a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics'>" +
            "InsertSizeMetrics</a> for detailed explanations of each metric." +
            "<hr />";

    private static final Log log = Log.getInstance(CollectInsertSizeMetrics.class);
    protected static final String Histogram_R_SCRIPT = "picard/analysis/insertSizeHistogram.R";

    @Argument(shortName="H", doc="File to write insert size Histogram chart to.")
    public File Histogram_FILE;

    @Argument(doc="Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. " +
            "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
            "artifacts to make the mean and sd grossly misleading regarding the real distribution.")
    public double DEVIATIONS = 10;

    @Argument(shortName="W", doc="Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. " +
            "Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.", optional=true)
    public Integer HISTOGRAM_WIDTH = null;

    @Argument(shortName="MW", doc="Minimum width of histogram plots. In the case when the histogram would otherwise be" +
            "truncated to a shorter range of sizes, the MIN_HISTOGRAM_WIDTH will enforce a minimum range.", optional=true)
    public Integer MIN_HISTOGRAM_WIDTH = null;

    @Argument(shortName="M", doc="When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this " +
            "percentage of overall reads. (Range: 0 to 1).")
    public float MINIMUM_PCT = 0.05f;

    @Argument(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Argument(doc="If true, also include reads marked as duplicates in the insert size histogram.")
    public boolean INCLUDE_DUPLICATES = false;

    // Calculates InsertSizeMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private InsertSizeMetricsCollector multiCollector;

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     *         to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errorMsgs = new ArrayList<String>();
        if (MINIMUM_PCT < 0 || MINIMUM_PCT > 0.5) {
            errorMsgs.add("MINIMUM_PCT was set to " + MINIMUM_PCT + ". It must be between 0 and 0.5 so all data categories don't get discarded.");
        }

        if (!checkRInstallation(Histogram_FILE != null)) {
            errorMsgs.add("R is not installed on this machine. It is required for creating the chart.");
        }

        return errorMsgs.isEmpty() ? null : errorMsgs.toArray(new String[errorMsgs.size()]);
    }

    @Override protected boolean usesNoRefReads() { return false; }

    @Override protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(Histogram_FILE);

        //Delegate actual collection to InsertSizeMetricCollector
        multiCollector = new InsertSizeMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), MINIMUM_PCT,
                HISTOGRAM_WIDTH, MIN_HISTOGRAM_WIDTH, DEVIATIONS, INCLUDE_DUPLICATES);
    }

    @Override protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
        multiCollector.acceptRecord(record, ref);
    }

    @Override protected void finish() {
        multiCollector.finish();

        final MetricsFile<InsertSizeMetrics, Integer> file = getMetricsFile();
        multiCollector.addAllLevelsToFile(file);

        if(file.getNumHistograms() == 0) {
            //can happen if user sets MINIMUM_PCT = 0.5, etc.
            log.warn("All data categories were discarded because they contained < " + MINIMUM_PCT +
                     " of the total aligned paired data.");
            final InsertSizeMetricsCollector.PerUnitInsertSizeMetricsCollector allReadsCollector = (InsertSizeMetricsCollector.PerUnitInsertSizeMetricsCollector) multiCollector.getAllReadsCollector();
            log.warn("Total mapped pairs in all categories: " + (allReadsCollector == null ? allReadsCollector : allReadsCollector.getTotalInserts()));
        }
        else  {
            file.write(OUTPUT);

            final List<String> plotArgs = new ArrayList<>();
            Collections.addAll(plotArgs, OUTPUT.getAbsolutePath(), Histogram_FILE.getAbsolutePath(), INPUT.getName());

            if (HISTOGRAM_WIDTH != null) {
                plotArgs.add(String.valueOf(HISTOGRAM_WIDTH));
            }
            else if (MIN_HISTOGRAM_WIDTH != null) {
                final int max = (int) file.getAllHistograms().stream().mapToDouble(Histogram::getMax).max().getAsDouble();
                plotArgs.add(String.valueOf(Math.max(max, MIN_HISTOGRAM_WIDTH)));
            }

            final int rResult = RExecutor.executeFromClasspath(Histogram_R_SCRIPT, plotArgs.toArray(new String[0]));
            if (rResult != 0) {
                throw new PicardException("R script " + Histogram_R_SCRIPT + " failed with return code " + rResult);
            }
        }
    }
}
