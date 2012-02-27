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

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Log;
import net.sf.picard.util.RExecutor;
import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Histogram;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.io.File;
import java.util.Arrays;


/**
 * Program to generate a data table and chart of mean quality by cycle from a
 * BAM file.  Works best on a single lane/run of data, but can be applied to
 * merged BAMs - the output may just be a little confusing.
 *
 * @author Tim Fennell
 */
public class MeanQualityByCycle extends SinglePassSamProgram {
    @Option(shortName="CHART", doc="A file (with .pdf extension) to write the chart to")
    public File CHART_OUTPUT;

    @Option(doc="If set to true calculate mean quality over aligned reads only")
    public boolean ALIGNED_READS_ONLY = false;

    @Option(doc="If set to true calculate mean quality over PF reads only")
    public boolean PF_READS_ONLY = false;

    private HistogramGenerator q  = new HistogramGenerator(false);
    private HistogramGenerator oq = new HistogramGenerator(true);

    private final Log log = Log.getInstance(MeanQualityByCycle.class);

    /** Required main method. */
    public static void main(String[] args) {
        System.exit(new MeanQualityByCycle().instanceMain(args));
    }

    private static class HistogramGenerator {
        final boolean useOriginalQualities;
        int maxLengthSoFar = 0;
        double[] firstReadTotalsByCycle  = new double[maxLengthSoFar];
        long[]   firstReadCountsByCycle  = new long[maxLengthSoFar];
        double[] secondReadTotalsByCycle = new double[maxLengthSoFar];
        long[]   secondReadCountsByCycle = new long[maxLengthSoFar];

        private HistogramGenerator(boolean useOriginalQualities) {
            this.useOriginalQualities = useOriginalQualities;
        }

        void addRecord(SAMRecord rec) {
            final byte[] quals = (useOriginalQualities ? rec.getOriginalBaseQualities() : rec.getBaseQualities());
            if (quals == null) return;

            final int length = quals.length;
            final boolean rc = rec.getReadNegativeStrandFlag();
            ensureArraysBigEnough(length+1);

            for (int i=0; i<length; ++i) {
                int cycle = rc ? length-i : i+1;

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

        private void ensureArraysBigEnough(int length) {
            if (length > maxLengthSoFar) {
                this.firstReadTotalsByCycle  = Arrays.copyOf(this.firstReadTotalsByCycle, length);
                this.firstReadCountsByCycle  = Arrays.copyOf(this.firstReadCountsByCycle, length);
                this.secondReadTotalsByCycle = Arrays.copyOf(this.secondReadTotalsByCycle , length);
                this.secondReadCountsByCycle = Arrays.copyOf(secondReadCountsByCycle, length);
                this.maxLengthSoFar = length;
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
                final int cycle = firstReadLength + i;

                if (secondReadCountsByCycle[i] > 0) {
                    meanQualities.increment(cycle, secondReadTotalsByCycle[i] / firstReadCountsByCycle[i]);
                }
            }

            return meanQualities;
        }

        boolean isEmpty() {
            return this.maxLengthSoFar == 0;
        }
    }


    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IoUtil.assertFileIsWritable(CHART_OUTPUT);
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // Skip unwanted records
        if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) return;
        if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) return;
        if (rec.getNotPrimaryAlignmentFlag()) return;

        q.addRecord(rec);
        oq.addRecord(rec);
    }

    @Override
    protected void finish() {
        // Generate a "histogram" of mean quality and write it to the file
        MetricsFile<?,Integer> metrics = getMetricsFile();
        metrics.addHistogram(q.getMeanQualityHistogram());
        if (!oq.isEmpty()) metrics.addHistogram(oq.getMeanQualityHistogram());
        metrics.write(OUTPUT);

        if (q.isEmpty() && oq.isEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        }
        else {
            // Now run R to generate a chart
            final int rResult = RExecutor.executeFromClasspath(
                    "net/sf/picard/analysis/meanQualityByCycle.R",
                    OUTPUT.getAbsolutePath(),
                    CHART_OUTPUT.getAbsolutePath(),
                    INPUT.getName());

            if (rResult != 0) {
                throw new PicardException("R script meanQualityByCycle.R failed with return code " + rResult);
            }
        }
    }
}

