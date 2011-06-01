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

import java.io.File;
import java.util.EnumMap;
import java.util.Map.Entry;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.sam.SamPairUtil;
import net.sf.picard.util.Histogram;
import net.sf.picard.util.Log;
import net.sf.picard.util.RExecutor;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.picard.sam.SamPairUtil.PairOrientation;

/**
 * Command line program to read non-duplicate insert sizes, create a histogram
 * and report distribution statistics.
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class CollectInsertSizeMetrics extends SinglePassSamProgram {
    private static final Log log = Log.getInstance(CollectInsertSizeMetrics.class);
    private static final String HISTOGRAM_R_SCRIPT = "net/sf/picard/analysis/insertSizeHistogram.R";
    // Usage and parameters
    @Usage
    public String USAGE = "Reads a SAM or BAM file and writes a file containing metrics about " +
            "the statistical distribution of insert size (excluding duplicates) " +
            "and generates a histogram plot.\n";

    @Option(shortName="H", doc="File to write insert size histogram chart to")
    public File HISTOGRAM_FILE;

    @Option(doc="Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. " +
            "This is done because insert size data typically includes enough anomolous values from chimeras and other " +
            "artifacts to make the mean and sd grossly misleading regardling the real distribution.")
    public double DEVIATIONS = 10;

    @Option(shortName="W", doc="Explicitly sets the histogram width, overriding the TAIL_LIMIT option. Also, when calculating " +
                               "mean and stdev, only bins <= HISTOGRAM_WIDTH will be included.", optional=true)
    public Integer HISTOGRAM_WIDTH = null;

    @Option(shortName="M", doc="When generating the histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this " +
            "percentage of overall reads. (Range: 0 to 1)")
    public float MINIMUM_PCT = 0.05f;


    final EnumMap<PairOrientation, Histogram<Integer>> histograms = new EnumMap<PairOrientation, Histogram<Integer>>(PairOrientation.class);

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
        IoUtil.assertFileIsWritable(OUTPUT);
        IoUtil.assertFileIsWritable(HISTOGRAM_FILE);

        histograms.put(PairOrientation.FR,     new Histogram<Integer>("insert_size", "fr_count"));
        histograms.put(PairOrientation.TANDEM, new Histogram<Integer>("insert_size", "tandem_count"));
        histograms.put(PairOrientation.RF,     new Histogram<Integer>("insert_size", "rf_count"));
    }

    @Override protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
        if (!record.getReadPairedFlag() ||
                record.getReadUnmappedFlag() ||
                record.getMateUnmappedFlag() ||
                record.getFirstOfPairFlag() ||
                record.getNotPrimaryAlignmentFlag() ||
                record.getDuplicateReadFlag() ||
                record.getInferredInsertSize() == 0) {
            return;
        }

        //add record to 1 of the 3 data categories
        final int insertSize = Math.abs(record.getInferredInsertSize());
        final PairOrientation orientation = SamPairUtil.getPairOrientation(record);
        histograms.get(orientation).increment(insertSize);
    }

    @Override protected void finish() {
        final MetricsFile<InsertSizeMetrics, Integer> file = getMetricsFile();

        double totalInserts = 0;
        for (final Histogram<Integer> h : this.histograms.values()) totalInserts += h.getCount();

        for(Entry<PairOrientation, Histogram<Integer>> entry : histograms.entrySet()) {
            final PairOrientation pairOrientation = entry.getKey();
            final Histogram<Integer> histogram = entry.getValue();
            final double total = histogram.getCount();

            // Only include a category if it has a sufficient percentage of the data in it
            if( total > totalInserts * MINIMUM_PCT ) {
                final InsertSizeMetrics metrics = new InsertSizeMetrics();
                metrics.PAIR_ORIENTATION   = pairOrientation;
                metrics.READ_PAIRS         = (long) total;
                metrics.MAX_INSERT_SIZE    = (int) histogram.getMax();
                metrics.MIN_INSERT_SIZE    = (int) histogram.getMin();
                metrics.MEDIAN_INSERT_SIZE = histogram.getMedian();
                metrics.MEDIAN_ABSOLUTE_DEVIATION = histogram.getMedianAbsoluteDeviation();

                final double median  = histogram.getMedian();
                double covered = 0;
                double low  = median;
                double high = median;

                while (low >= histogram.getMin() || high <= histogram.getMax()) {
                    final Histogram<Integer>.Bin lowBin = histogram.get((int) low);
                    if (lowBin != null) covered += lowBin.getValue();

                    if (low != high) {
                        final Histogram<Integer>.Bin highBin = histogram.get((int) high);
                        if (highBin != null) covered += highBin.getValue();
                    }

                    final double percentCovered = covered / total;
                    final int distance = (int) (high - low) + 1;
                    if (percentCovered >= 0.1  && metrics.WIDTH_OF_10_PERCENT == 0) metrics.WIDTH_OF_10_PERCENT = distance;
                    if (percentCovered >= 0.2  && metrics.WIDTH_OF_20_PERCENT == 0) metrics.WIDTH_OF_20_PERCENT = distance;
                    if (percentCovered >= 0.3  && metrics.WIDTH_OF_30_PERCENT == 0) metrics.WIDTH_OF_30_PERCENT = distance;
                    if (percentCovered >= 0.4  && metrics.WIDTH_OF_40_PERCENT == 0) metrics.WIDTH_OF_40_PERCENT = distance;
                    if (percentCovered >= 0.5  && metrics.WIDTH_OF_50_PERCENT == 0) metrics.WIDTH_OF_50_PERCENT = distance;
                    if (percentCovered >= 0.6  && metrics.WIDTH_OF_60_PERCENT == 0) metrics.WIDTH_OF_60_PERCENT = distance;
                    if (percentCovered >= 0.7  && metrics.WIDTH_OF_70_PERCENT == 0) metrics.WIDTH_OF_70_PERCENT = distance;
                    if (percentCovered >= 0.8  && metrics.WIDTH_OF_80_PERCENT == 0) metrics.WIDTH_OF_80_PERCENT = distance;
                    if (percentCovered >= 0.9  && metrics.WIDTH_OF_90_PERCENT == 0) metrics.WIDTH_OF_90_PERCENT = distance;
                    if (percentCovered >= 0.99 && metrics.WIDTH_OF_99_PERCENT == 0) metrics.WIDTH_OF_99_PERCENT = distance;

                    --low;
                    ++high;
                }

                // Trim the histogram down to get rid of outliers that would make the chart useless.
                final Histogram<Integer> trimmedHisto = histogram; //alias it
                if (HISTOGRAM_WIDTH == null) {
                    HISTOGRAM_WIDTH = (int) (metrics.MEDIAN_INSERT_SIZE + (DEVIATIONS * metrics.MEDIAN_ABSOLUTE_DEVIATION));
                }

                trimmedHisto.trimByWidth(HISTOGRAM_WIDTH);

                metrics.MEAN_INSERT_SIZE = trimmedHisto.getMean();
                metrics.STANDARD_DEVIATION = trimmedHisto.getStandardDeviation();

                file.addHistogram(trimmedHisto);
                file.addMetric(metrics);
            }
        }

        if(file.getNumHistograms() == 0) {
            //can happen if user sets MINIMUM_PCT = 0.5, etc.
            log.warn("All data categories were discarded because they contained < " + MINIMUM_PCT +
                     " of the total aligned paired data.");
            log.warn("Total mapped pairs in all categories: " + totalInserts);
        }
        else  {
            file.write(OUTPUT);

            final int rResult;
            if(HISTOGRAM_WIDTH == null) {
                rResult = RExecutor.executeFromClasspath(
                    HISTOGRAM_R_SCRIPT,
                    OUTPUT.getAbsolutePath(),
                    HISTOGRAM_FILE.getAbsolutePath(),
                    INPUT.getName());
            } else {
                rResult = RExecutor.executeFromClasspath(
                    HISTOGRAM_R_SCRIPT,
                    OUTPUT.getAbsolutePath(),
                    HISTOGRAM_FILE.getAbsolutePath(),
                    INPUT.getName(),
                    String.valueOf( HISTOGRAM_WIDTH ) ); //HISTOGRAM_WIDTH is passed because R automatically sets histogram width to the last
                                                         //bin that has data, which may be less than HISTOGRAM_WIDTH and confuse the user.
            }

            if (rResult != 0) {
                throw new PicardException("R script " + HISTOGRAM_R_SCRIPT + " failed with return code " + rResult);
            }
        }
    }
}
