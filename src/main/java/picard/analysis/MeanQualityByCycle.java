/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

import java.io.File;
import java.util.Arrays;
import java.util.List;


/**
 * Program to generate a data table and chart of mean quality by cycle from a
 * BAM file.  Works best on a single lane/run of data, but can be applied to
 * merged BAMs - the output may just be a little confusing.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = MeanQualityByCycle.USAGE_SUMMARY + MeanQualityByCycle.USAGE_DETAILS,
        oneLineSummary = MeanQualityByCycle.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class MeanQualityByCycle extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "Collect mean quality by cycle.";
    static final String USAGE_DETAILS = "This tool generates a data table and chart of mean quality by cycle from a BAM file. It is " +
            "intended to be used on a single lane or a read group's worth of data, but can be applied to merged BAMs if needed. " +
            "<br /><br />" +
            "This metric gives an overall snapshot of sequencing machine performance. For most types of sequencing data, the output " +
            "is expected to show a slight reduction in overall base quality scores towards the end of each read. Spikes in quality within " +
            "reads are not expected and may indicate that technical problems occurred during sequencing." +
            "<br /><br />" +
            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar MeanQualityByCycle \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=mean_qual_by_cycle.txt \\<br />" +
            "      CHART=mean_qual_by_cycle.pdf" +
            "</pre>" +
            "<hr />";
    @Argument(shortName="CHART", doc="A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Argument(doc="If set to true, calculate mean quality over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Argument(doc="If set to true calculate mean quality over PF reads only.")
    public boolean PF_READS_ONLY = false;

    private final HistogramGenerator q  = new HistogramGenerator(false);
    private final HistogramGenerator oq = new HistogramGenerator(true);

    /**
     * A subtitle for the plot, usually corresponding to a library.
     */
    private String plotSubtitle = "";

    private final Log log = Log.getInstance(MeanQualityByCycle.class);

    /** Required main method. */
    public static void main(String[] args) {
        System.exit(new MeanQualityByCycle().instanceMain(args));
    }

    private static final class HistogramGenerator {
        final boolean useOriginalQualities;
        int maxLengthSoFar = 0;
        double[] firstReadTotalsByCycle  = new double[maxLengthSoFar];
        long[]   firstReadCountsByCycle  = new long[maxLengthSoFar];
        double[] secondReadTotalsByCycle = new double[maxLengthSoFar];
        long[]   secondReadCountsByCycle = new long[maxLengthSoFar];

        private HistogramGenerator(final boolean useOriginalQualities) {
            this.useOriginalQualities = useOriginalQualities;
        }

        void addRecord(final SAMRecord rec) {
            final byte[] quals = (useOriginalQualities ? rec.getOriginalBaseQualities() : rec.getBaseQualities());
            if (quals == null) return;

            final int length = quals.length;
            final boolean rc = rec.getReadNegativeStrandFlag();
            ensureArraysBigEnough(length+1);

            for (int i=0; i<length; ++i) {
                final int cycle = rc ? length-i : i+1;

                if (rec.getReadPairedFlag() && rec.getSecondOfPairFlag()) {
                    secondReadTotalsByCycle[cycle] += quals[i];
                    secondReadCountsByCycle[cycle] += 1;
                }
                else {
                    firstReadTotalsByCycle[cycle] += quals[i];
                    firstReadCountsByCycle[cycle] += 1;
                }
            }
        }

        private void ensureArraysBigEnough(final int length) {
            if (length > maxLengthSoFar) {
                firstReadTotalsByCycle  = Arrays.copyOf(firstReadTotalsByCycle, length);
                firstReadCountsByCycle  = Arrays.copyOf(firstReadCountsByCycle, length);
                secondReadTotalsByCycle = Arrays.copyOf(secondReadTotalsByCycle , length);
                secondReadCountsByCycle = Arrays.copyOf(secondReadCountsByCycle, length);
                maxLengthSoFar = length;
            }
        }

        Histogram<Integer> getMeanQualityHistogram() {
            final String label = useOriginalQualities ? "MEAN_ORIGINAL_QUALITY" : "MEAN_QUALITY";
            final Histogram<Integer> meanQualities = new Histogram<Integer>("CYCLE", label);

            int firstReadLength = 0;

            for (int cycle=0; cycle < firstReadTotalsByCycle.length; ++cycle) {
                if (firstReadTotalsByCycle[cycle] > 0) {
                    meanQualities.increment(cycle, firstReadTotalsByCycle[cycle] / firstReadCountsByCycle[cycle]);
                    firstReadLength = cycle;
                }
            }

            for (int i=0; i< secondReadTotalsByCycle.length; ++i) {
                if (secondReadCountsByCycle[i] > 0) {
                    final int cycle = firstReadLength + i;
                    meanQualities.increment(cycle, secondReadTotalsByCycle[i] / secondReadCountsByCycle[i]);
                }
            }

            return meanQualities;
        }

        boolean isEmpty() {
            return maxLengthSoFar == 0;
        }
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        // If we're working with a single library, assign that library's name
        // as a suffix to the plot title
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // Skip unwanted records
        if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) return;
        if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) return;
        if (rec.isSecondaryOrSupplementary()) return;

        q.addRecord(rec);
        oq.addRecord(rec);
    }

    @Override
    protected void finish() {
        // Generate a "Histogram" of mean quality and write it to the file
        final MetricsFile<?,Integer> metrics = getMetricsFile();
        metrics.addHistogram(q.getMeanQualityHistogram());
        if (!oq.isEmpty()) metrics.addHistogram(oq.getMeanQualityHistogram());
        metrics.write(OUTPUT);

        if (q.isEmpty() && oq.isEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        }
        else {
            // Now run R to generate a chart
            final int rResult = RExecutor.executeFromClasspath(
                    "picard/analysis/meanQualityByCycle.R",
                    OUTPUT.getAbsolutePath(),
                    CHART_OUTPUT.getAbsolutePath(),
                    INPUT.getName(),
                    plotSubtitle);

            if (rResult != 0) {
                throw new PicardException("R script meanQualityByCycle.R failed with return code " + rResult);
            }
        }
    }
}

