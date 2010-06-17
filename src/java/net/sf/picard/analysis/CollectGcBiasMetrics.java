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

import net.sf.picard.analysis.GcBiasDetailMetrics;
import net.sf.picard.analysis.GcBiasSummaryMetrics;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.util.SequenceUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.*;
import java.text.NumberFormat;

import net.sf.picard.util.RExecutor;
import net.sf.picard.util.QualityUtil;

/**
 * Tool to collect information about GC bias in the reads in a given BAM file. Computes
 * the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%
 * and counts the number of read starts in each GC bin.  What is output and plotted is
 * the "normalized coverage" in each bin - i.e. the number of reads per window normalized
 * to the average number of reads per window across the whole genome.
 *
 * @author Tim Fennell
 */
public class CollectGcBiasMetrics extends CommandLineProgram {
    /** The location of the R script to do the plotting. */
    private static final String R_SCRIPT = "edu/mit/broad/picard/sam/gcBias.R";

    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference sequence fasta file.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The BAM or SAM file containing aligned reads.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The text file to write the metrics table to.")
    public File OUTPUT;

    @Option(shortName="CHART", doc="The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Option(doc="The text file to write summary metrics to.", optional=true)
    public File SUMMARY_OUTPUT;

    @Option(doc="The size of windows on the genome that are used to bin reads.")
    public int WINDOW_SIZE = 100;

    @Option(doc="For summary metrics, exclude GC windows that include less than this fraction of the genome.")
    public double MINIMUM_GENOME_FRACTION = 0.00001;

    // Used to keep track of the total clusters as this is kinda important for bias
    private int totalClusters = 0;
    private int totalAlignedReads = 0;

    /** Stock main method. */
    public static void main(String[] args) {
        System.exit(new CollectGcBiasMetrics().instanceMain(args));
    }

    protected int doWork() {
        IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        IoUtil.assertFileIsWritable(CHART_OUTPUT);
        if (SUMMARY_OUTPUT != null) IoUtil.assertFileIsWritable(SUMMARY_OUTPUT);

        // Histograms to track the number of windows at each GC, and the number of read starts
        // at windows of each GC
        int[] windowsByGc = new int[101];
        int[] readsByGc   = new int[101];
        long[] basesByGc  = new long[101];
        long[] errorsByGc = new long[101];

        SAMFileReader sam = new SAMFileReader(INPUT);
        PeekableIterator<SAMRecord> iterator = new PeekableIterator<SAMRecord>(sam.iterator());
        ReferenceSequenceFile referenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);

        {
            // Check that the sequence dictionaries match if present
            SAMSequenceDictionary referenceDictionary= referenceFile.getSequenceDictionary();
            SAMSequenceDictionary samFileDictionary = sam.getFileHeader().getSequenceDictionary();
            if (referenceDictionary != null && samFileDictionary != null) {
                SequenceUtil.assertSequenceDictionariesEqual(referenceDictionary, samFileDictionary);
            }
        }

        ////////////////////////////////////////////////////////////////////////////
        // Loop over the reference and the reads and calculate the basic metrics
        ////////////////////////////////////////////////////////////////////////////
        ReferenceSequence ref = null;
        while ((ref = referenceFile.nextSequence()) != null) {
            final byte[] refBases = ref.getBases();
            StringUtil.toUpperCase(refBases);
            final int refLength = refBases.length;
            final int lastWindowStart = refLength - WINDOW_SIZE;

            byte[] gc = null;
            for (int countSeq=0 ; iterator.hasNext() && iterator.peek().getReferenceIndex() == ref.getContigIndex() ; countSeq++) {
                SAMRecord rec = iterator.next();
                if (countSeq == 0) { 
                    gc = calculateAllGcs(refBases, windowsByGc, lastWindowStart) ;
                    }
                if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++this.totalClusters;

                if (!rec.getReadUnmappedFlag()) {
                    int pos = rec.getReadNegativeStrandFlag() ? rec.getAlignmentEnd() - WINDOW_SIZE : rec.getAlignmentStart();
                    ++this.totalAlignedReads;

                    if (pos > 0) {
                        int windowGc = gc[pos];

                        if (windowGc >= 0) {
                            ++readsByGc[windowGc];
                            basesByGc[windowGc]  += rec.getReadLength();
                            errorsByGc[windowGc] += SequenceUtil.countMismatches(rec, refBases);
                        }
                    }
                }
            }

            System.out.println("Processed: " + ref.getName());
        }

        // Finish up the reads, presumably all unaligned
        while (iterator.hasNext()) {
            SAMRecord rec = iterator.next();
            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++this.totalClusters;

        }

        /////////////////////////////////////////////////////////////////////////////
        // Synthesize the normalized coverage metrics and write it all out to a file
        /////////////////////////////////////////////////////////////////////////////
        MetricsFile<GcBiasDetailMetrics,?> metricsFile = getMetricsFile();
        double totalWindows = sum(windowsByGc);
        double totalReads   = sum(readsByGc);
        double meanReadsPerWindow = totalReads / totalWindows;
        double minimumWindowsToCountInSummary = totalWindows * this.MINIMUM_GENOME_FRACTION;

        for (int i=0; i<windowsByGc.length; ++i) {
            if (windowsByGc[i] == 0) continue;

            GcBiasDetailMetrics m = new GcBiasDetailMetrics();
            m.GC = i;
            m.WINDOWS             = windowsByGc[i];
            m.READ_STARTS         = readsByGc[i];
            if (errorsByGc[i] > 0) m.MEAN_BASE_QUALITY = QualityUtil.getPhredScoreFromObsAndErrors(basesByGc[i], errorsByGc[i]);
            m.NORMALIZED_COVERAGE = (m.READ_STARTS / (double) m.WINDOWS) / meanReadsPerWindow;
            m.ERROR_BAR_WIDTH     = (Math.sqrt(m.READ_STARTS) / (double) m.WINDOWS) / meanReadsPerWindow;

            metricsFile.addMetric(m);
        }

        metricsFile.write(OUTPUT);

        // Synthesize the high level metrics
        if (SUMMARY_OUTPUT != null) {
            MetricsFile<GcBiasSummaryMetrics,?> summaryMetricsFile = getMetricsFile();
            GcBiasSummaryMetrics summary = new GcBiasSummaryMetrics();
            summary.WINDOW_SIZE    = this.WINDOW_SIZE;
            summary.TOTAL_CLUSTERS = this.totalClusters;
            summary.ALIGNED_READS  = this.totalAlignedReads;
            summary.LOW_GC_BIAS  = calculateSummaryMetric(metricsFile.getMetrics(), minimumWindowsToCountInSummary, 0,  33);
            summary.MID_GC_BIAS  = calculateSummaryMetric(metricsFile.getMetrics(), minimumWindowsToCountInSummary, 34,  66);
            summary.HIGH_GC_BIAS = calculateSummaryMetric(metricsFile.getMetrics(), minimumWindowsToCountInSummary, 67, 100);
            summary.TOTAL_BIAS   = calculateSummaryMetric(metricsFile.getMetrics(), minimumWindowsToCountInSummary,  0, 100);
            summary.JAFFE_BIAS_METRIC = calculateJaffeMetric(metricsFile.getMetrics(), minimumWindowsToCountInSummary, 0, 100);

            summaryMetricsFile.addMetric(summary);
            if (SUMMARY_OUTPUT != null) summaryMetricsFile.write(SUMMARY_OUTPUT);
        }

        // Plot the results
        NumberFormat fmt = NumberFormat.getIntegerInstance();
        fmt.setGroupingUsed(true);
        String subtitle = "Total clusters: " + fmt.format(this.totalClusters) +
                          ", Aligned reads: " + fmt.format(this.totalAlignedReads);
        RExecutor.executeFromClasspath(R_SCRIPT,
                                       OUTPUT.getAbsolutePath(),
                                       CHART_OUTPUT.getAbsolutePath(),
                                       INPUT.getName().replace(".duplicates_marked", "").replace(".aligned.bam", ""),
                                       subtitle,
                                       String.valueOf(WINDOW_SIZE));

        return 0;
    }


    /** Sums the values in an int[]. */
    private double sum(final int[] values) {
        final int length = values.length;
        double total = 0;
        for (int i=0; i<length; ++i) {
            total += values[i];
        }

        return total;
    }

    /**
     * Calculates a summary metric that is intended to confer whether or not a set of reads exhibits
     * significant bias.
     */
    private double calculateSummaryMetric(Collection<GcBiasDetailMetrics> details, double minimumWindows, int minGc, int maxGc) {
        double retval = 0;
        double totalBins = 0;

        Map<Integer, GcBiasDetailMetrics> metricsByGc = new TreeMap<Integer, GcBiasDetailMetrics>();
        for (GcBiasDetailMetrics detail : details) metricsByGc.put(detail.GC, detail);

        int[][] loopvars = new int[][] {new int[] {50,-1, -1}, new int[] {51, 101, 1}};

        for (int[] limits : loopvars) {
            double last=0;

            for (int gc=limits[0]; gc != limits[1]; gc += limits[2]) {
                GcBiasDetailMetrics detail = metricsByGc.get(gc);

                if (detail != null && detail.WINDOWS >= minimumWindows) {
                    double current = detail.NORMALIZED_COVERAGE + detail.ERROR_BAR_WIDTH;
                    if (current == 0) current = last;
                    else last = current;

                    if (detail.GC >= minGc && detail.GC <= maxGc) {
                        totalBins += 1;
                        if (current > 0 && current < 1) {
                            retval += (1d/current);
                        }
                    }
                }
            }
        }

        if (totalBins > 0) {
            return retval / totalBins;
        }
        else {
            return 0;
        }
    }

    private double calculateJaffeMetric(Collection<GcBiasDetailMetrics> details, double minimumWindows, int minGc, int maxGc) {
        double jaffeMetric = 0;
        double binsCounted = 0;
        GcBiasDetailMetrics previous = null;

        for (GcBiasDetailMetrics detail : details) {
//            if (detail.WINDOWS < minimumWindows && detail.GC >= minGc && detail.GC <= maxGc) {
                if (previous != null) {
                    double h = detail.GC - previous.GC;
                    double fx0 = Math.pow((1d - previous.NORMALIZED_COVERAGE), 2);
                    double fx1 = Math.pow((1d - detail.NORMALIZED_COVERAGE), 2);
                    jaffeMetric += (h/2.0)*(fx0 + fx1);

//                    binsCounted++;
//                }

                previous = detail;
            }
        }

        return 100d * Math.sqrt(jaffeMetric) / details.size();
    }

    /** Calculcate all the GC values for all windows. */
    private byte[] calculateAllGcs(byte [] refBases, int [] windowsByGc, int lastWindowStart) {
        final int refLength = refBases.length;
        byte[] gc = new byte[refLength + 1];
        final CalculateGcState state = new CalculateGcState();
        for (int i=1; i<lastWindowStart; ++i) {
            final int windowEnd = i + WINDOW_SIZE;
            final int windowGc = calculateGc(refBases, i, windowEnd, state) ; 
            gc[i] = (byte) windowGc;
            if (windowGc != -1) windowsByGc[windowGc]++;
        }
        return gc;
    }

    /**
     * Calculates GC as a number from 0 to 100 in the specified window. If the window includes
     * more than five no-calls then -1 is returned.
     */
    private int calculateGc(final byte[] bases, final int startIndex, final int endIndex, CalculateGcState state) {
        if (state.init) {
            state.init = false ;
            state.gcCount = 0;
            state.nCount  = 0;
            for (int i=startIndex; i<endIndex; ++i) {
                final byte base = bases[i];
                if (base == 'G' || base == 'C') ++state.gcCount;
                else if (base == 'N') ++state.nCount;
            }
        } else {
            byte newBase = bases[endIndex-1];
            if (newBase == 'G' || newBase == 'C') ++state.gcCount;
            else if (newBase == 'N') ++state.nCount;

            if (state.priorBase == 'G' || state.priorBase == 'C') --state.gcCount;
            else if (state.priorBase == 'N') --state.nCount;
        }
        state.priorBase = bases[startIndex];
        if (state.nCount > 4) return -1;
        else return (state.gcCount * 100) / (endIndex - startIndex);
    }

    /** Keeps track of current GC calculation state. */
    class CalculateGcState {
        boolean init = true ;
        int nCount ;
        int gcCount ;
        byte priorBase ;
    }
}
