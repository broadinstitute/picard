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
import picard.PicardException;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

import java.io.File;
import java.util.Arrays;
import java.util.List;


/**
 * tabulate insert size counts for both duplate and non-duplate reads
 */
@CommandLineProgramProperties(
        usage = " " +
                " ",
        usageShort = "Writes insert size counts for both duplate and non-duplate reads",
        programGroup = Metrics.class
)
public class DuplicationByInsertLength extends SinglePassSamProgram {

    @Option(shortName="CHART", doc="A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Option(doc="If set to true, calculate mean quality over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Option(doc="If set to true calculate mean quality over PF reads only.")
    public boolean PF_READS_ONLY = false;

    private final HistogramGenerator q  = new HistogramGenerator(false);

    /**
     * A subtitle for the plot, usually corresponding to a library.
     */
    private String plotSubtitle = "";

    /** Required main method. */
    public static void main(String[] args) {
        System.exit(new MeanQualityByCycle().instanceMain(args));
    }

    private static class HistogramGenerator {
        final boolean useOriginalQualities;
        Histogram<Integer> insertHistDup = new Histogram<Integer>("Insert_length", "duplicate_reads");
        Histogram<Integer> insertHistNotDup = new Histogram<Integer>("Insert_length", "unique_reads");

        private HistogramGenerator(final boolean useOriginalQualities) {
            this.useOriginalQualities = useOriginalQualities;
        }

        void addRecord(final SAMRecord rec) {
            int iSz = rec.getInferredInsertSize();
            boolean isDup = rec.getDuplicateReadFlag();
            if (isDup){
                insertHistDup.increment(Math.abs(iSz));
            }
            else {
                insertHistNotDup.increment(Math.abs(iSz));
            }
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
    }

    @Override
    protected void finish() {
        // Generate a "Histogram" of insert size length
        final MetricsFile<?,Integer> metrics = getMetricsFile();
        metrics.addHistogram(q.insertHistDup);
        metrics.addHistogram(q.insertHistNotDup);
        metrics.write(OUTPUT);
        /*else {
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
        }*/
    }
}

