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

package net.sf.picard.analysis;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Histogram;
import net.sf.picard.util.RExecutor;
import net.sf.samtools.util.SequenceUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.PicardException;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * Charts quality score distribution within a BAM file.
 *
 * @author Tim Fennell
 */
public class QualityScoreDistribution extends CommandLineProgram {
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input BAM file to process")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="A file to write the table of qualities to")
    public File OUTPUT;

    @Option(shortName="CHART", doc="A file (with .pdf extension) to write the chart to")
    public File CHART_OUTPUT;

    @Option(doc="If set to true calculate mean quality over aligned reads only")
    public boolean ALIGNED_READS_ONLY = false;

    @Option(shortName="PF", doc="If set to true calculate mean quality over PF reads only")
    public boolean PF_READS_ONLY = false;

    @Option(doc="If set to true, include quality for no-call bases in the distribution")
    public boolean INCLUDE_NO_CALLS = false;

    @Option(doc="Stop after processing N reads. Mostly for debugging purposes.")
    public int STOP_AFTER = 0;

    /** Required main method. */
    public static void main(final String[] args) {
        System.exit(new QualityScoreDistribution().instanceMain(args));
    }

    /**
     * Does all the work of calculating the mean qualities and writing out the files.
     */
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        IoUtil.assertFileIsWritable(CHART_OUTPUT);

        final SAMFileReader in = new SAMFileReader(INPUT);
        final Histogram<Byte> qHisto  = new Histogram<Byte>("QUALITY", "COUNT_OF_Q");
        final Histogram<Byte> oqHisto = new Histogram<Byte>("QUALITY", "COUNT_OF_OQ");

        // Read through the SAM file and aggregate total quality and observations by cycle
        int processed = 0;
        for (final SAMRecord rec : in) {
            // Skip unwanted records
            if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) continue;
            if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) continue;

            final byte[] bases = rec.getReadBases();
            final byte[] quals = rec.getBaseQualities();
            final byte[] oq    = rec.getOriginalBaseQualities();

            final int length = quals.length;

            for (int i=0; i<length; ++i) {
                if (INCLUDE_NO_CALLS || !SequenceUtil.isNoCall(bases[i])) {
                    qHisto.increment(quals[i]);
                    if (oq != null) oqHisto.increment(oq[i]);
                }
            }

            if (STOP_AFTER > 0 && ++processed > STOP_AFTER) break;
        }

        final MetricsFile<?,Byte> metrics = getMetricsFile();
        metrics.addHistogram(qHisto);
        if (!oqHisto.isEmpty()) metrics.addHistogram(oqHisto);
        metrics.write(OUTPUT);

        // Now run R to generate a chart
        final int rResult = RExecutor.executeFromClasspath(
                "edu/mit/broad/picard/quality/qualityScoreDistribution.R",
                OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                INPUT.getName());

        if (rResult != 0) {
            throw new PicardException("R script qualityScoreDistribution.R failed with return code " + rResult);
        }

        return 0;
    }
}
