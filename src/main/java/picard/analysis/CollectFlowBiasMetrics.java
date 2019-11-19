/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.io.output.NullOutputStream;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.LinkedList;
import java.util.NoSuchElementException;
import java.util.Random;

import static htsjdk.samtools.util.SequenceUtil.N;

/**
 *
 */
@CommandLineProgramProperties(
        summary = CollectFlowBiasMetrics.USAGE_SUMMARY + CollectFlowBiasMetrics.USAGE_DETAILS,
        oneLineSummary = CollectFlowBiasMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CollectFlowBiasMetrics extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Collect metrics regarding Flow. ";
    static final String USAGE_DETAILS = "";
    /**
     * The location of the R script to do the plotting.
     */
    private static final String R_SCRIPT = "picard/analysis/flowBias.R";

    // Usage and parameters

    @Argument(shortName = "CHART", doc = "The PDF file to render the chart to.", optional = true)
    public File DETAILED_OUTPUT;

    @Argument(shortName = "S", doc = "The text file to write summary metrics to.")
    public File SUMMARY_OUTPUT;

    @Argument(shortName = "FLOW", doc = "The flow order used in sequencing (will be used repeatedly)")
    public String FLOW_ORDER;

    @Argument(shortName = "SIZE", doc = "The number of flows that were used")
    public Integer FLOW_LENGTH;

    @Argument
    public Integer MAX_READ_TO_CONSIDER = 2000;

    @Argument
    public double READ_PROBABILITY = 1D;

    @Argument
    long STOP_AFTER = 0;
    @Argument
    int PROGRESS_INTERVAL = 100000;

    private final long randomSeed = 42;

    private final int[] localCoverage = new int[MAX_READ_TO_CONSIDER];
    private String FLOW_ORDER_RC;

    /////////////////////////////////////////////////////////////////////////////
    // Setup calculates windowsByGc for the entire reference. Must be done at
    // startup to avoid missing reference contigs in the case of small files
    // that may not have reads aligning to every reference contig.
    /////////////////////////////////////////////////////////////////////////////

    @Override
    protected String[] customCommandLineValidation() {
        if (DETAILED_OUTPUT != null) {
            IOUtil.assertFileIsWritable(DETAILED_OUTPUT);
        }
        IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        return super.customCommandLineValidation();
    }

    @Override
    protected boolean requiresReference() {
        return true;
    }

    final Log log = Log.getInstance(CollectFlowBiasMetrics.class);
    private ProgressLogger progressLogger;
    private Random random;

    private final Histogram<Integer> coverageHist = new Histogram<>();

    @Override
    protected int doWork() {
        try (PrintStream detailedOutput = DETAILED_OUTPUT != null ?
                new PrintStream(DETAILED_OUTPUT) :
                new PrintStream(NullOutputStream.NULL_OUTPUT_STREAM)) {

            detailedOutput.println("#CHROM\tPOS\tEND\tNAME\tCOUNT\tCOUNT_POS\tCOUNT_NEG");
            random = new Random(randomSeed);

            progressLogger = new ProgressLogger(log, PROGRESS_INTERVAL, "positions", "examined");

            FLOW_ORDER_RC = SequenceUtil.reverseComplement(FLOW_ORDER);
            coverageHist.setValueLabel("Coverage");
            coverageHist.setValueLabel("Count");
            try (ReferenceSequenceFileWalker referenceSequenceWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE)) {
                for (SAMSequenceRecord samSequenceRecord : referenceSequenceWalker.getSequenceDictionary().getSequences()) {

                    final ReferenceSequence referenceSequence = referenceSequenceWalker.get(samSequenceRecord.getSequenceIndex());
                    log.info("examining sequence " + referenceSequence.getName());
                    updateMetrics(referenceSequence, detailedOutput);
                    if (STOP_AFTER != 0 && progressLogger.getCount() > STOP_AFTER) {
                        break;
                    }
                }
            } catch (
                    IOException e) {
                throw new PicardException("Error while reading reference " + REFERENCE_SEQUENCE, e);
            }
        } catch (FileNotFoundException e) {
            throw new PicardException("Error while writing to " + DETAILED_OUTPUT, e);
        }

        MetricsFile<?, Integer> metricsFile = getMetricsFile();
        metricsFile.addHistogram(coverageHist);
        metricsFile.write(SUMMARY_OUTPUT);
        return 0;
    }

    private void updateMetrics(final ReferenceSequence referenceSequence, final PrintStream detailedOutput) {
        final int N = referenceSequence.length();
        final byte[] referenceBytesCopy = new byte[N];
        System.arraycopy(referenceSequence.getBases(), 0, referenceBytesCopy, 0, N);

        SequenceUtil.reverseComplement(referenceBytesCopy);
        final ReferenceSequence referenceSequenceRC = new ReferenceSequence(referenceSequence.getName(), referenceSequence.getContigIndex(), referenceBytesCopy);

        final byte[] baseString = referenceSequence.getBases();
        final byte[] baseStringRC = referenceSequenceRC.getBases();
        final LinkedList<Integer> readLengths = new LinkedList<>();
        final LinkedList<Integer> readLengthsRC = new LinkedList<>();
        final byte[] flowOrder = FLOW_ORDER.getBytes();
        final byte[] flowOrderRC = FLOW_ORDER_RC.getBytes();
        final byte[] referenceSlice = new byte[MAX_READ_TO_CONSIDER];
        final byte[] referenceSliceRC = new byte[MAX_READ_TO_CONSIDER];

        for (int i = MAX_READ_TO_CONSIDER; i < N - MAX_READ_TO_CONSIDER; i++) {
            progressLogger.record(referenceSequence.getName(), i);

            if (STOP_AFTER != 0 && progressLogger.getCount() > STOP_AFTER) {
                return;
            }

            System.arraycopy(baseString, i, referenceSlice, 0, MAX_READ_TO_CONSIDER);
            final int lengthGivenFlow = readLengthGivenFlow(referenceSlice, flowOrder, FLOW_LENGTH);
            if (random.nextDouble() < READ_PROBABILITY) {
                readLengths.addFirst(lengthGivenFlow);
            } else {
                readLengths.addFirst(-1); //-1 means that theres no read here.
            }

            System.arraycopy(baseStringRC, N - 1 - i - MAX_READ_TO_CONSIDER, referenceSliceRC, 0, MAX_READ_TO_CONSIDER);

            final int lengthGivenFlowRC = readLengthGivenFlow(referenceSliceRC, flowOrderRC, FLOW_LENGTH);
            if (random.nextDouble() < READ_PROBABILITY) {
                readLengthsRC.addFirst(lengthGivenFlowRC);
            } else {
                readLengthsRC.addFirst(-1); //-1 means that there's no read here.
            }
            final Integer coverageFor = countCoverage(readLengths);
            final Integer coverageRev = countCoverage(readLengthsRC);
            final int localCoverage;
            localCoverage = coverageFor + coverageRev;
            if (DETAILED_OUTPUT != null) {
                detailedOutput.print(referenceSequence.getName());
                detailedOutput.print("\t");
                detailedOutput.print(i);
                detailedOutput.print("\t");
                detailedOutput.print(i);
                detailedOutput.print("\t");
                detailedOutput.print(localCoverage);
                detailedOutput.print("\t");
                detailedOutput.print(coverageFor);
                detailedOutput.print("\t");
                detailedOutput.print(coverageRev);
                detailedOutput.println();
            }
            coverageHist.increment(localCoverage);
        }
    }

    private Integer countCoverage(final LinkedList<Integer> readLengths) {

        int i = 0;
        int depth = 0;
        int readLengthI;
        while (i < readLengths.size() && ((readLengthI = readLengths.get(i)) == -1 || i < readLengthI)) {
            if (readLengthI != -1) {
                depth++;
            }
            i++;
        }
        while (readLengths.size() > i) {
            readLengths.removeLast();
        }
        return depth;
    }

    static int readLengthGivenFlow(final byte[] referenceSlice, final byte[] flowOrder, final int flowLength) {
        int i = 0;
        for (int j = 0; j < flowLength; j++) {
            if (SequenceUtil.basesEqual(N, referenceSlice[i])) {
                return 0;
            }
            while (SequenceUtil.basesEqual(referenceSlice[i], flowOrder[j % flowOrder.length])) {
                i++;
            }
        }
        return i;
    }

//    /////////////////////////////////////////////////////////////////////////////
//    // Write out all levels of normalized coverage metrics to a file
//    /////////////////////////////////////////////////////////////////////////////
//    @Override
//    protected void finish() {
//        multiCollector.finish();
//        writeResultsToFiles();
//    }
//
//    private void writeResultsToFiles() {
//        final MetricsFile<GcBiasMetrics, Integer> file = getMetricsFile();
//        final MetricsFile<GcBiasDetailMetrics, ?> detailMetricsFile = getMetricsFile();
//        final MetricsFile<GcBiasSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
//        multiCollector.addAllLevelsToFile(file);
//        final List<GcBiasMetrics> gcBiasMetricsList = file.getMetrics();
//        for (final GcBiasMetrics gcbm : gcBiasMetricsList) {
//            final List<GcBiasDetailMetrics> gcDetailList = gcbm.DETAILS.getMetrics();
//            for (final GcBiasDetailMetrics d : gcDetailList) {
//                detailMetricsFile.addMetric(d);
//            }
//            summaryMetricsFile.addMetric(gcbm.SUMMARY);
//        }
//        detailMetricsFile.write(OUTPUT);
//        summaryMetricsFile.write(SUMMARY_OUTPUT);
//
//        final NumberFormat fmt = NumberFormat.getIntegerInstance();
//        fmt.setGroupingUsed(true);
//        RExecutor.executeFromClasspath(R_SCRIPT,
//                OUTPUT.getAbsolutePath(),
//                SUMMARY_OUTPUT.getAbsolutePath(),
//                CHART_OUTPUT.getAbsolutePath(),
//                String.valueOf(SCAN_WINDOW_SIZE));
//    }
}


