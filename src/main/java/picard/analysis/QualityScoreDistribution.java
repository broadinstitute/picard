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
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

import java.io.File;
import java.util.List;

/**
 * Charts quality score distribution within a BAM file.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = QualityScoreDistribution.USAGE_SUMMARY + QualityScoreDistribution.USAGE_DETAILS,
        oneLineSummary = QualityScoreDistribution.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class QualityScoreDistribution extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "Chart the distribution of quality scores.  ";
    static final String USAGE_DETAILS = "<p>This tool is used for determining the overall 'quality' for a library in a given run. To " +
            "that effect, it outputs a chart and tables indicating the range of quality scores and the total numbers of bases " +
            "corresponding to those scores. Options include plotting the distribution of all of the reads, only the aligned reads, " +
            "or reads that have passed the Illumina Chastity filter thresholds as described " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=6329'>here</a>.</p>" +
            "" +
            "<h4>Note on base quality score options</h4>" +
            "If the quality score of read bases has been modified in a previous data processing step such as GATK  " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=44'>Base Recalibration</a> " +
            "and an OQ tag is available, this tool can be set to plot the OQ value as well as the primary quality value for the " +
            "evaluation. <br />" +
            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +
            "<h4>Usage Example:</h4>" +
            "<pre>" +
            "java -jar picard.jar QualityScoreDistribution \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=qual_score_dist.txt \\<br />" +
            "      CHART=qual_score_dist.pdf" +
            "</pre>" +
            "<hr />";
    @Argument(shortName="CHART", doc="A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Argument(doc="If set to true calculate mean quality over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Argument(shortName="PF", doc="If set to true calculate mean quality over PF reads only.")
    public boolean PF_READS_ONLY = false;

    @Argument(doc="If set to true, include quality for no-call bases in the distribution.")
    public boolean INCLUDE_NO_CALLS = false;

    private final long[] qCounts  = new long[128];
    private final long[] oqCounts = new long[128];

    /**
     * A subtitle for the plot, usually corresponding to a library.
     */
    private String plotSubtitle = "";

    private final Log log = Log.getInstance(QualityScoreDistribution.class);

    /** Required main method. */
    public static void main(final String[] args) {
        System.exit(new QualityScoreDistribution().instanceMain(args));
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(CHART_OUTPUT);

        // If we're working with a single library, assign that library's name
        // as a suffix to the plot title
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            this.plotSubtitle = readGroups.get(0).getLibrary();
            if (null == this.plotSubtitle) this.plotSubtitle = "";
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // Skip unwanted records
        if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) return;
        if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) return;
        if (rec.isSecondaryOrSupplementary()) return;

        final byte[] bases = rec.getReadBases();
        final byte[] quals = rec.getBaseQualities();
        final byte[] oq    = rec.getOriginalBaseQualities();

        final int length = quals.length;

        for (int i=0; i<length; ++i) {
            if (INCLUDE_NO_CALLS || !SequenceUtil.isNoCall(bases[i])) {
                qCounts[quals[i]]++;
                if (oq != null) oqCounts[oq[i]]++;
            }
        }
    }

    @Override
    protected void finish() {
        // Built the Histograms out of the long[]s
        final Histogram<Byte> qHisto  = new Histogram<Byte>("QUALITY", "COUNT_OF_Q");
        final Histogram<Byte> oqHisto = new Histogram<Byte>("QUALITY", "COUNT_OF_OQ");

        for (int i=0; i< qCounts.length; ++i) {
            if (qCounts[i]  > 0) qHisto.increment( (byte) i, (double) qCounts[i]);
            if (oqCounts[i] > 0) oqHisto.increment((byte) i, (double) oqCounts[i]);
        }

        final MetricsFile<?,Byte> metrics = getMetricsFile();
        metrics.addHistogram(qHisto);
        if (!oqHisto.isEmpty()) metrics.addHistogram(oqHisto);
        metrics.write(OUTPUT);
        if (qHisto.isEmpty() && oqHisto.isEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        }
        else {
            // Now run R to generate a chart
            final int rResult = RExecutor.executeFromClasspath(
                    "picard/analysis/qualityScoreDistribution.R",
                    OUTPUT.getAbsolutePath(),
                    CHART_OUTPUT.getAbsolutePath(),
                    INPUT.getName(),
                    this.plotSubtitle);

            if (rResult != 0) {
                throw new PicardException("R script qualityScoreDistribution.R failed with return code " + rResult);
            }
        }
    }
}
