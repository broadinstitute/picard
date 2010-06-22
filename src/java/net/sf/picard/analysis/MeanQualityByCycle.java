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

import net.sf.picard.util.RExecutor;
import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Histogram;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.io.File;


/**
 * Program to generate a data table and chart of mean quality by cycle from a
 * BAM file.  Works best on a single lane/run of data, but can be applied to
 * merged BAMs - the output may just be a little confusing.
 *
 * @author Tim Fennell
 */
public class MeanQualityByCycle extends CommandLineProgram {
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input BAM file to process")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="A file to write the table of qualities to")
    public File OUTPUT;

    @Option(shortName="CHART", doc="A file (with .pdf extension) to write the chart to")
    public File CHART_OUTPUT;

    @Option(doc="If set to true calculate mean quality over aligned reads only")
    public boolean ALIGNED_READS_ONLY = false;

    @Option(doc="If set to true calculate mean quality over PF reads only")
    public boolean PF_READS_ONLY = false;

    @Option(doc="Stop after processing N reads, mainly for debugging.") public int STOP_AFTER = 0;

    /** Required main method. */
    public static void main(String[] args) {
        System.exit(new MeanQualityByCycle().instanceMain(args));
    }

    private static class HistogramGenerator {
        final boolean useOriginalQualities;
        Histogram<Integer> totalQ = new Histogram<Integer>();
        Histogram<Integer> countOfBases = new Histogram<Integer>();

        private HistogramGenerator(boolean useOriginalQualities) {
            this.useOriginalQualities = useOriginalQualities;
        }

        void addRecord(SAMRecord rec) {
            final byte[] quals = (useOriginalQualities ? rec.getOriginalBaseQualities() : rec.getBaseQualities());
            if (quals == null) return;

            final int length = quals.length;
            final boolean rc = !rec.getReadUnmappedFlag() && rec.getReadNegativeStrandFlag();

            for (int i=0; i<length; ++i) {
                int cycle = rc ? length-i : i+1;

                if (rec.getReadPairedFlag() && rec.getSecondOfPairFlag()) {
                    cycle = rec.getReadLength() + cycle;
                }

                totalQ.increment(cycle, quals[i]);
                countOfBases.increment(cycle);
            }
        }

        Histogram<Integer> getMeanQualityHistogram() {
            final String label = useOriginalQualities ? "MEAN_ORIGINAL_QUALITY" : "MEAN_QUALITY";
            final Histogram<Integer> meanQualities = new Histogram<Integer>("CYCLE", label);

            for (Integer cycle : countOfBases.keySet()) {
                meanQualities.increment(cycle, (totalQ.get(cycle).getValue() / countOfBases.get(cycle).getValue()));
            }

            return meanQualities;
        }

        boolean isEmpty() {
            return totalQ.isEmpty();
        }
    }

    /**
     * Does all the work of calcuating the mean qualities and writing out the files.
     */
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        IoUtil.assertFileIsWritable(CHART_OUTPUT);

        SAMFileReader in = new SAMFileReader(INPUT);
        HistogramGenerator q  = new HistogramGenerator(false);
        HistogramGenerator oq = new HistogramGenerator(true);

        // Read through the SAM file and aggregate total quality and observations by cycle
        int i = 0;
        for (SAMRecord rec : in) {
            // Skip unwanted records
            if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) continue;
            if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) continue;

            q.addRecord(rec);
            oq.addRecord(rec);

            if (STOP_AFTER > 0 && ++i >= STOP_AFTER) break;
        }

        // Generate a "histogram" of mean quality and write it to the file
        MetricsFile<?,Integer> metrics = getMetricsFile();
        metrics.addHistogram(q.getMeanQualityHistogram());
        if (!oq.isEmpty()) metrics.addHistogram(oq.getMeanQualityHistogram());
        metrics.write(OUTPUT);

        // Now run R to generate a chart
        final int rResult = RExecutor.executeFromClasspath(
                "net/sf/picard/analysis/meanQualityByCycle.R",
                OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                INPUT.getName());

        if (rResult != 0) {
            throw new PicardException("R script meanQualityByCycle.R failed with return code " + rResult);
        }

        return 0;
    }
}

