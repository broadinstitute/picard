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
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.RExecutor;

import java.io.File;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.List;

/**
 * Tool to collect information about GC bias in the reads in a given BAM file. Computes
 * the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%
 * and counts the number of read starts in each GC bin.  What is output and plotted is
 * the "normalized coverage" in each bin - i.e. the number of reads per window normalized
 * to the average number of reads per window across the whole genome.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Tool to collect information about GC bias in the reads in a given BAM file. Computes" +
                " the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%" +
                " and counts the number of read starts in each GC bin.  What is output and plotted is" +
                " the \"normalized coverage\" in each bin - i.e. the number of reads per window normalized" +
                " to the average number of reads per window across the whole genome..\n",
        usageShort = "Collects information about GC bias in the reads in the provided SAM or BAM",
        programGroup = Metrics.class
)
public class CollectGcBiasMetrics extends CommandLineProgram {
    /** The location of the R script to do the plotting. */
    private static final String R_SCRIPT = "picard/analysis/gcBias.R";

    // Usage and parameters

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "The reference sequence fasta file.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The BAM or SAM file containing aligned reads.  " +
            "Must be coordinate-sorted.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The text file to write the metrics table to.")
    public File OUTPUT;

    @Option(shortName = "CHART", doc = "The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Option(shortName = "S", doc = "The text file to write summary metrics to.", optional = true)
    public File SUMMARY_OUTPUT;

    @Option(doc = "The size of windows on the genome that are used to bin reads.")
    public int WINDOW_SIZE = 100;

    @Option(doc = "For summary metrics, exclude GC windows that include less than this fraction of the genome.")
    public double MINIMUM_GENOME_FRACTION = 0.00001;

    @Option(doc = "If true, assume that the input file is coordinate sorted, even if the header says otherwise.",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = false;

    @Option(shortName = "BS", doc = "Whether the SAM or BAM file consists of bisulfite sequenced reads.  ")
    public boolean IS_BISULFITE_SEQUENCED = false;

    private static final Log log = Log.getInstance(CollectGcBiasMetrics.class);

    // Used to keep track of the total clusters as this is kinda important for bias
    private int totalClusters = 0;
    private int totalAlignedReads = 0;

    /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new CollectGcBiasMetrics().instanceMain(args));
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        if (SUMMARY_OUTPUT != null) IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);

        // Histograms to track the number of windows at each GC, and the number of read starts
        // at windows of each GC
        final int[] windowsByGc = new int[101];
        final int[] readsByGc = new int[101];
        final long[] basesByGc = new long[101];
        final long[] errorsByGc = new long[101];

        final SamReader sam = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        if (!ASSUME_SORTED && sam.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException("Header of input file " + INPUT.getAbsolutePath() + " indicates that it is not coordinate sorted.  " +
                    "If you believe the records are in coordinate order, pass option ASSUME_SORTED=true.  If not, sort the file with SortSam.");
        }
        final PeekableIterator<SAMRecord> iterator = new PeekableIterator<SAMRecord>(sam.iterator());
        final ReferenceSequenceFile referenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);

        {
            // Check that the sequence dictionaries match if present
            final SAMSequenceDictionary referenceDictionary = referenceFile.getSequenceDictionary();
            final SAMSequenceDictionary samFileDictionary = sam.getFileHeader().getSequenceDictionary();
            if (referenceDictionary != null && samFileDictionary != null) {
                SequenceUtil.assertSequenceDictionariesEqual(referenceDictionary, samFileDictionary);
            }
        }

        ////////////////////////////////////////////////////////////////////////////
        // Loop over the reference and the reads and calculate the basic metrics
        ////////////////////////////////////////////////////////////////////////////
        ReferenceSequence ref = null;
        final ProgressLogger
                progress = new ProgressLogger(log);
        while ((ref = referenceFile.nextSequence()) != null) {
            final byte[] refBases = ref.getBases();
            StringUtil.toUpperCase(refBases);
            final int refLength = refBases.length;
            final int lastWindowStart = refLength - WINDOW_SIZE;

            final byte[] gc = calculateAllGcs(refBases, windowsByGc, lastWindowStart);
            while (iterator.hasNext() && iterator.peek().getReferenceIndex() == ref.getContigIndex()) {
                final SAMRecord rec = iterator.next();

                if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++this.totalClusters;

                if (!rec.getReadUnmappedFlag()) {
                    final int pos = rec.getReadNegativeStrandFlag() ? rec.getAlignmentEnd() - WINDOW_SIZE : rec.getAlignmentStart();
                    ++this.totalAlignedReads;

                    if (pos > 0) {
                        final int windowGc = gc[pos];

                        if (windowGc >= 0) {
                            ++readsByGc[windowGc];
                            basesByGc[windowGc] += rec.getReadLength();
                            errorsByGc[windowGc] +=
                                    SequenceUtil.countMismatches(rec, refBases, IS_BISULFITE_SEQUENCED) +
                                            SequenceUtil.countInsertedBases(rec) + SequenceUtil.countDeletedBases(rec);
                        }
                    }
                }

                progress.record(rec);
            }

            log.info("Processed reference sequence: " + ref.getName());
        }

        // Finish up the reads, presumably all unaligned
        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++this.totalClusters;

        }

        /////////////////////////////////////////////////////////////////////////////
        // Synthesize the normalized coverage metrics and write it all out to a file
        /////////////////////////////////////////////////////////////////////////////
        final MetricsFile<GcBiasDetailMetrics, ?> metricsFile = getMetricsFile();
        final double totalWindows = sum(windowsByGc);
        final double totalReads = sum(readsByGc);
        final double meanReadsPerWindow = totalReads / totalWindows;
        final double minimumWindowsToCountInSummary = totalWindows * this.MINIMUM_GENOME_FRACTION;

        for (int i = 0; i < windowsByGc.length; ++i) {
            if (windowsByGc[i] == 0) continue;

            GcBiasDetailMetrics m = new GcBiasDetailMetrics();
            m.GC = i;
            m.WINDOWS = windowsByGc[i];
            m.READ_STARTS = readsByGc[i];
            if (errorsByGc[i] > 0) m.MEAN_BASE_QUALITY = QualityUtil.getPhredScoreFromObsAndErrors(basesByGc[i], errorsByGc[i]);
            m.NORMALIZED_COVERAGE = (m.READ_STARTS / (double) m.WINDOWS) / meanReadsPerWindow;
            m.ERROR_BAR_WIDTH = (Math.sqrt(m.READ_STARTS) / (double) m.WINDOWS) / meanReadsPerWindow;

            metricsFile.addMetric(m);
        }

        metricsFile.write(OUTPUT);

        // Synthesize the high level metrics
        if (SUMMARY_OUTPUT != null) {
            final MetricsFile<GcBiasSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
            final GcBiasSummaryMetrics summary = new GcBiasSummaryMetrics();
            summary.WINDOW_SIZE = this.WINDOW_SIZE;
            summary.TOTAL_CLUSTERS = this.totalClusters;
            summary.ALIGNED_READS = this.totalAlignedReads;
            calculateDropoutMetrics(metricsFile.getMetrics(), summary);

            summaryMetricsFile.addMetric(summary);
            summaryMetricsFile.write(SUMMARY_OUTPUT);
        }

        // Plot the results
        final NumberFormat fmt = NumberFormat.getIntegerInstance();
        fmt.setGroupingUsed(true);
        final String subtitle = "Total clusters: " + fmt.format(this.totalClusters) +
                ", Aligned reads: " + fmt.format(this.totalAlignedReads);
        String title = INPUT.getName().replace(".duplicates_marked", "").replace(".aligned.bam", "");

        // Qualify the title with the library name iff it's for a single sample
        final List<SAMReadGroupRecord> readGroups = sam.getFileHeader().getReadGroups();
        if (readGroups.size() == 1) {
            title += "." + readGroups.get(0).getLibrary();
        }
        RExecutor.executeFromClasspath(R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                title,
                subtitle,
                String.valueOf(WINDOW_SIZE));

        return 0;
    }

    /** Sums the values in an int[]. */
    private double sum(final int[] values) {
        final int length = values.length;
        double total = 0;
        for (int i = 0; i < length; ++i) {
            total += values[i];
        }

        return total;
    }

    /** Calculates the Illumina style AT and GC dropout numbers. */
    private void calculateDropoutMetrics(final Collection<GcBiasDetailMetrics> details,
                                         final GcBiasSummaryMetrics summary) {
        // First calculate the totals
        double totalReads = 0;
        double totalWindows = 0;

        for (final GcBiasDetailMetrics detail : details) {
            totalReads += detail.READ_STARTS;
            totalWindows += detail.WINDOWS;
        }

        double atDropout = 0;
        double gcDropout = 0;

        for (final GcBiasDetailMetrics detail : details) {
            final double relativeReads = detail.READ_STARTS / totalReads;
            final double relativeWindows = detail.WINDOWS / totalWindows;
            final double dropout = (relativeWindows - relativeReads) * 100;

            if (dropout > 0) {
                if (detail.GC <= 50) atDropout += dropout;
                if (detail.GC >= 50) gcDropout += dropout;
            }
        }

        summary.AT_DROPOUT = atDropout;
        summary.GC_DROPOUT = gcDropout;
    }

    /** Calculcate all the GC values for all windows. */
    private byte[] calculateAllGcs(final byte[] refBases, final int[] windowsByGc, final int lastWindowStart) {
        final int refLength = refBases.length;
        final byte[] gc = new byte[refLength + 1];
        final CalculateGcState state = new CalculateGcState();
        for (int i = 1; i < lastWindowStart; ++i) {
            final int windowEnd = i + WINDOW_SIZE;
            final int windowGc = calculateGc(refBases, i, windowEnd, state);
            gc[i] = (byte) windowGc;
            if (windowGc != -1) windowsByGc[windowGc]++;
        }
        return gc;
    }

    /**
     * Calculates GC as a number from 0 to 100 in the specified window. If the window includes
     * more than five no-calls then -1 is returned.
     */
    private int calculateGc(final byte[] bases, final int startIndex, final int endIndex, final CalculateGcState state) {
        if (state.init) {
            state.init = false;
            state.gcCount = 0;
            state.nCount = 0;
            for (int i = startIndex; i < endIndex; ++i) {
                final byte base = bases[i];
                if (base == 'G' || base == 'C') ++state.gcCount;
                else if (base == 'N') ++state.nCount;
            }
        } else {
            final byte newBase = bases[endIndex - 1];
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
        boolean init = true;
        int nCount;
        int gcCount;
        byte priorBase;
    }
}
