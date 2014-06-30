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
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.analysis.directed.InsertSizeMetricsCollector;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

import java.io.File;
import java.util.Set;

/**
 * Command line program to read non-duplicate insert sizes, create a Histogram
 * and report distribution statistics.
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
@CommandLineProgramProperties(
        usage = "Reads a SAM or BAM file and writes a file containing metrics about " +
                "the statistical distribution of insert size (excluding duplicates) " +
                "and generates a Histogram plot.",
        usageShort = "Writes insert size distribution metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectInsertSizeMetrics extends SinglePassSamProgram {
    private static final Log log = Log.getInstance(CollectInsertSizeMetrics.class);
    private static final String Histogram_R_SCRIPT = "picard/analysis/insertSizeHistogram.R";

    @Option(shortName="H", doc="File to write insert size Histogram chart to.")
    public File Histogram_FILE;

    @Option(doc="Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. " +
            "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
            "artifacts to make the mean and sd grossly misleading regarding the real distribution.")
    public double DEVIATIONS = 10;

    @Option(shortName="W", doc="Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. " +
            "Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.", optional=true)
    public Integer Histogram_WIDTH = null;

    @Option(shortName="M", doc="When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this " +
            "percentage of overall reads. (Range: 0 to 1).")
    public float MINIMUM_PCT = 0.05f;

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    // Calculates InsertSizeMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private InsertSizeMetricsCollector multiCollector;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectInsertSizeMetrics().instanceMainWithExit(argv);
    }

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
         if (MINIMUM_PCT < 0 || MINIMUM_PCT > 0.5) {
             return new String[]{"MINIMUM_PCT was set to " + MINIMUM_PCT + ". It must be between 0 and 0.5 so all data categories don't get discarded."};
         }

         return super.customCommandLineValidation();
    }

    @Override protected boolean usesNoRefReads() { return false; }

    @Override protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(Histogram_FILE);

        //Delegate actual collection to InsertSizeMetricCollector
        multiCollector = new InsertSizeMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), MINIMUM_PCT, Histogram_WIDTH, DEVIATIONS);
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

            final int rResult;
            if(Histogram_WIDTH == null) {
                rResult = RExecutor.executeFromClasspath(
                    Histogram_R_SCRIPT,
                    OUTPUT.getAbsolutePath(),
                    Histogram_FILE.getAbsolutePath(),
                    INPUT.getName());
            } else {
                rResult = RExecutor.executeFromClasspath(
                    Histogram_R_SCRIPT,
                    OUTPUT.getAbsolutePath(),
                    Histogram_FILE.getAbsolutePath(),
                    INPUT.getName(),
                    String.valueOf( Histogram_WIDTH ) ); //Histogram_WIDTH is passed because R automatically sets Histogram width to the last
                                                         //bin that has data, which may be less than Histogram_WIDTH and confuse the user.
            }

            if (rResult != 0) {
                throw new PicardException("R script " + Histogram_R_SCRIPT + " failed with return code " + rResult);
            }
        }
    }
}
