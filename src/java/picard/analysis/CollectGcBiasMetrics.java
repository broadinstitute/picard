/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.analysis.directed.GcBiasMetricsCollector;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.metrics.GcBiasMetrics;
import picard.util.RExecutor;

import java.io.File;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Tool to collect information about GC bias in the reads in a given BAM file. Computes
 * the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%
 * and counts the number of read starts in each GC bin.  What is output and plotted is
 * the "normalized coverage" in each bin - i.e. the number of reads per window normalized
 * to the average number of reads per window across the whole genome.
 *
 * @author Tim Fennell
 * edited by Kylee Bergin
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
public class CollectGcBiasMetrics extends SinglePassSamProgram {
    /** The location of the R script to do the plotting. */
    private static final String R_SCRIPT = "picard/analysis/gcBias.R";

    // Usage and parameters

    @Option(shortName = "CHART", doc = "The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Option(shortName = "S", doc = "The text file to write summary metrics to.", optional = true)
    public File SUMMARY_OUTPUT;

    @Option(doc = "The size of windows on the genome that are used to bin reads.")
    public int WINDOW_SIZE = 100;

    @Option(doc = "For summary metrics, exclude GC windows that include less than this fraction of the genome.")
    public double MINIMUM_GENOME_FRACTION = 0.00001;

    @Option(shortName = "BS", doc = "Whether the SAM or BAM file consists of bisulfite sequenced reads.")
    public boolean IS_BISULFITE_SEQUENCED = false;

    @Option(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    // Calculates GcBiasMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private GcBiasMetricsCollector multiCollector;

    //windowSize is the size of the scanning window that goes over the reference
    private final int windowSize = WINDOW_SIZE;
    final int[] windowsByGc = new int[WINDOWS];

    // Histograms to track the number of windows at each GC, and the number of read starts
    // at windows of each GC. Need 101 to get from 0-100.
    private static final int WINDOWS = 101;

    //Hash map of gc[] with reference name as key
    private final Map<String, byte[]> gcByRef = new HashMap<String, byte[]>();

    ////////////////////////////////////////////////////////////////////////////
    // Stock main method
    ////////////////////////////////////////////////////////////////////////////
    public static void main(final String[] args) {
        System.exit(new CollectGcBiasMetrics().instanceMain(args));
    }

    /////////////////////////////////////////////////////////////////////////////
    // Setup calculates gc[] for the reference. Must be done at startup to avoid
    // missing reference sequences in the case of small files that may
    // not have reads aligning to every reference sequence
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);

        if (SUMMARY_OUTPUT != null) IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);

        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
        ReferenceSequence ref;

        while ((ref = refFile.nextSequence()) != null) {
            final byte[] refBases = ref.getBases();
            final String refName = ref.getName();
            StringUtil.toUpperCase(refBases);
            final int refLength = refBases.length;
            final int lastWindowStart = refLength - windowSize;
            final byte[] gc = calculateAllGcs(refBases, windowsByGc, lastWindowStart);
            gcByRef.put(refName, gc);
        }
        //Delegate actual collection to GcBiasMetricCollector
        multiCollector = new GcBiasMetricsCollector(METRIC_ACCUMULATION_LEVEL, gcByRef, windowsByGc, header.getReadGroups(), windowSize, IS_BISULFITE_SEQUENCED);
    }

    ////////////////////////////////////////////////////////////////////////////
    // MultiCollector acceptRead
    ////////////////////////////////////////////////////////////////////////////
    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        multiCollector.acceptRecord(rec, ref);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Write out all levels of normalized coverage metrics to a file
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void finish() {
        multiCollector.finish();
        final MetricsFile<GcBiasMetrics, Integer> file = getMetricsFile();
        final MetricsFile<GcBiasDetailMetrics, ?> detailMetricsFile = getMetricsFile();
        final MetricsFile<GcBiasSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
        multiCollector.addAllLevelsToFile(file);
        final List<GcBiasMetrics> gcBiasMetricsList = file.getMetrics();
        for(final GcBiasMetrics gcbm : gcBiasMetricsList){
            final List<GcBiasDetailMetrics> gcDetailList = gcbm.DETAILS.getMetrics();
            for(final GcBiasDetailMetrics d : gcDetailList) {
                detailMetricsFile.addMetric(d);
            }
            summaryMetricsFile.addMetric(gcbm.SUMMARY);
        }
        detailMetricsFile.write(OUTPUT);
        summaryMetricsFile.write(SUMMARY_OUTPUT);

        final NumberFormat fmt = NumberFormat.getIntegerInstance();
        fmt.setGroupingUsed(true);
        RExecutor.executeFromClasspath(R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                SUMMARY_OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                String.valueOf(WINDOW_SIZE));
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculcate all the GC values for all windows
    /////////////////////////////////////////////////////////////////////////////
    private byte[] calculateAllGcs(final byte[] refBases, final int[] windowsByGc, final int lastWindowStart) {
        final CalculateGcState state = new CalculateGcState();
        final int refLength = refBases.length;
        final byte[] gc = new byte[refLength + 1];
        for (int i = 1; i < lastWindowStart; ++i) {
            final int windowEnd = i + windowSize;
            final int windowGc = calculateGc(refBases, i, windowEnd, state);
            gc[i] = (byte) windowGc;
            if (windowGc != -1) windowsByGc[windowGc]++;
        }
        return gc;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculates GC as a number from 0 to 100 in the specified window.
    // If the window includes more than five no-calls then -1 is returned.
    /////////////////////////////////////////////////////////////////////////////
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

    /////////////////////////////////////////////////////////////////////////////
    // Keeps track of current GC calculation state
    /////////////////////////////////////////////////////////////////////////////
    class CalculateGcState {
        boolean init = true;
        int nCount;
        int gcCount;
        byte priorBase;
    }
}


