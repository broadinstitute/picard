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

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.metrics.GcBiasMetrics;
import picard.metrics.MultiLevelCollector;
import picard.metrics.PerUnitMetricCollector;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashMap;

/** Calculates GC Bias Metrics on multiple levels
 *  Created by kbergin on 3/23/15.
 */
public class GcBiasMetricsCollector extends MultiLevelCollector<GcBiasMetrics, Integer, GcBiasCollectorArgs> {
    // Histograms to track the number of windows at each GC, and the number of read starts
    // at windows of each GC. Need 101 to get from 0-100.
    private final int scanWindowSize;
    private final boolean bisulfite;
    private int[] windowsByGc = new int[BINS];
    private static final int BINS = 101;
    //Use to calculate additional results without duplicates
    private boolean ignoreDuplicates;

    //will hold the relevant gc information per contig
    private byte [] gc = null;
    private int referenceIndex = -1;
    private byte [] refBases = null;
    private static final Log log = Log.getInstance(GcBiasMetricsCollector.class);

    public GcBiasMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final int[] windowsByGc,
                                  final List<SAMReadGroupRecord> samRgRecords, final int scanWindowSize,
                                  final boolean bisulfite) {
        this(accumulationLevels, windowsByGc, samRgRecords, scanWindowSize, bisulfite, false);
    }

    public GcBiasMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final int[] windowsByGc,
                                  final List<SAMReadGroupRecord> samRgRecords, final int scanWindowSize,
                                  final boolean bisulfite, final boolean ignoreDuplicates) {
        this.scanWindowSize = scanWindowSize;
        this.bisulfite = bisulfite;
        this.windowsByGc = windowsByGc;
        this.ignoreDuplicates = ignoreDuplicates;
        setup(accumulationLevels, samRgRecords);
    }

    /////////////////////////////////////////////////////////////////////////////
    // This method is called once Per samRecord
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected GcBiasCollectorArgs makeArg(final SAMRecord rec, final ReferenceSequence ref) {
        return new GcBiasCollectorArgs(rec, ref);
    }

    /////////////////////////////////////////////////////////////////////////////
    //Make a GcBiasCollector with the given arguments
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected PerUnitMetricCollector<GcBiasMetrics, Integer, GcBiasCollectorArgs> makeChildCollector(final String sample, final String library, final String readGroup) {
        return new PerUnitGcBiasMetricsCollector(sample, library, readGroup);
    }

    /////////////////////////////////////////////////////////////////////////////
    //A collector for individual GcBiasMetrics for a given SAMPLE or SAMPLE/LIBRARY
    //or SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels)
    /////////////////////////////////////////////////////////////////////////////
    public class PerUnitGcBiasMetricsCollector implements PerUnitMetricCollector<GcBiasMetrics, Integer, GcBiasCollectorArgs> {
        private final Map<String, GcObject> gcData;
        // Additional object to store data without duplicates (null if option ALSO_IGNORE_DUPLICATES is not specified)
        private Map<String, GcObject> gcDataNonDups;
        private final String sample;
        private final String library;
        private final String readGroup;
        private static final String allReads = "All_Reads";
        final static String ACCUMULATION_LEVEL_ALL_READS = "All Reads";
        final static String ACCUMULATION_LEVEL_LIBRARY = "Library";
        final static String ACCUMULATION_LEVEL_SAMPLE = "Sample";
        final static String ACCUMULATION_LEVEL_READ_GROUP = "Read Group";
        final static String READS_USED_ALL = "ALL";
        final static String READS_USED_UNIQUE = "UNIQUE";
        private int logCounter;

        /////////////////////////////////////////////////////////////////////////////
        //Records the accumulation level for each level of collection and initializes
        // a GcObject for this accumulation level
        /////////////////////////////////////////////////////////////////////////////
        public PerUnitGcBiasMetricsCollector(final String sample, final String library, final String readGroup) {
            this.sample = sample;
            this.library = library;
            this.readGroup = readGroup;
            this.gcData = prepareGcData();
            if (ignoreDuplicates) {
                this.gcDataNonDups = prepareGcData();
            }
        }

        /////////////////////////////////////////////////////////////////////////////
        //Takes each record and sends them to addRead to calculate gc metrics for
        // that read for each accumulation level
        /////////////////////////////////////////////////////////////////////////////
        @Override
        public void acceptRecord(final GcBiasCollectorArgs args) {
            final SAMRecord rec = args.getRec();
            if (logCounter < 100 && rec.getReadBases().length == 0) {
                log.warn("Omitting read " + rec.getReadName() + " with '*' in SEQ field.");
                if (++logCounter == 100) {
                    log.warn("There are more than 100 reads with '*' in SEQ field in file.");
                }
                return;
            }
            if (!rec.getReadUnmappedFlag()) {
                if (referenceIndex != rec.getReferenceIndex() || gc == null) {
                    final ReferenceSequence ref = args.getRef();
                    refBases = ref.getBases();
                    StringUtil.toUpperCase(refBases);
                    final int refLength = refBases.length;
                    final int lastWindowStart = refLength - scanWindowSize;
                    gc = GcBiasUtils.calculateAllGcs(refBases, lastWindowStart, scanWindowSize);
                    referenceIndex = rec.getReferenceIndex();
                }

                addReadToGcData(rec, this.gcData);
                if (ignoreDuplicates && !rec.getDuplicateReadFlag()) {
                    addReadToGcData(rec, this.gcDataNonDups);
                }
            } else {
                updateTotalClusters(rec, this.gcData);
                if (ignoreDuplicates && !rec.getDuplicateReadFlag()) {
                    updateTotalClusters(rec, this.gcDataNonDups);
                }
            }
        }

        @Override
        public void finish() {}

        /////////////////////////////////////////////////////////////////////////////
        //Called to add metrics to the output file for each level of collection
        // these metrics are used for graphing gc bias in R script
        /////////////////////////////////////////////////////////////////////////////
        @Override
        public void addMetricsToFile(final MetricsFile<GcBiasMetrics, Integer> file) {
            addGcDataToFile(file, this.gcData, true);
            if (ignoreDuplicates) {
                addGcDataToFile(file, this.gcDataNonDups, false);
            }
        }

        private void updateTotalClusters(final SAMRecord rec, final Map<String, GcObject> gcData) {
            for (final Map.Entry<String, GcObject> entry : gcData.entrySet()) {
                final GcObject gcCur = entry.getValue();
                if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++gcCur.totalClusters;
            }
        }

        private Map<String, GcObject> prepareGcData() {
            final String prefix;
            final Map<String, GcObject> gcData = new HashMap<>();
            if (this.readGroup != null) {
                prefix = this.readGroup;
            } else if (this.library != null) {
                prefix = this.library;
            } else if (this.sample != null) {
                prefix = this.sample;
            } else {
                prefix = allReads;
            }
            gcData.put(prefix, new GcObject());
            return gcData;
        }

        private void addReadToGcData(final SAMRecord rec, final Map<String, GcObject> gcData) {
            final String type;
            String group;
            if (this.readGroup != null) {
                type = this.readGroup;
                group = ACCUMULATION_LEVEL_READ_GROUP;
            } else if (this.library != null) {
                type = this.library;
                group = ACCUMULATION_LEVEL_LIBRARY;
            } else if (this.sample != null) {
                type = this.sample;
                group = ACCUMULATION_LEVEL_SAMPLE;
            } else {
                type = allReads;
                group = ACCUMULATION_LEVEL_ALL_READS;
            }
            addRead(gcData.get(type), rec, group, gc, refBases);
        }

        /////////////////////////////////////////////////////////////////////////////
        // Sums the values in an int[].
        /////////////////////////////////////////////////////////////////////////////
        private double sum(final int[] values) {
            final int length = values.length;
            double total = 0;
            for (int i = 0; i < length; i++) {
                total += values[i];
            }

            return total;
        }

        private void addGcDataToFile(final MetricsFile<GcBiasMetrics, Integer> file, final Map<String, GcObject> gcData,
                                     final boolean includeDuplicates) {
            for (final Map.Entry<String, GcObject> entry : gcData.entrySet()) {
                final GcObject gcCur = entry.getValue();
                final String gcType = entry.getKey();

                final int[] readsByGc = gcCur.readsByGc;
                final long[] errorsByGc = gcCur.errorsByGc;
                final long[] basesByGc = gcCur.basesByGc;
                final long totalClusters = gcCur.totalClusters;
                final long totalAlignedReads = gcCur.totalAlignedReads;
                final String group = gcCur.group;

                final GcBiasMetrics metrics = new GcBiasMetrics();

                final double totalWindows = sum(windowsByGc);
                final double totalReads = sum(readsByGc);
                final double meanReadsPerWindow = totalReads / totalWindows;

                if (totalAlignedReads > 0) {
                    for (int i = 0; i < windowsByGc.length; ++i) {
                        final GcBiasDetailMetrics detail = new GcBiasDetailMetrics();
                        detail.GC = i;
                        detail.WINDOWS = windowsByGc[i];
                        detail.READ_STARTS = readsByGc[i];
                        if (errorsByGc[i] > 0) {
                            detail.MEAN_BASE_QUALITY = QualityUtil.getPhredScoreFromObsAndErrors(basesByGc[i], errorsByGc[i]);
                        }
                        if (windowsByGc[i] != 0) {
                            detail.NORMALIZED_COVERAGE = (detail.READ_STARTS / (double) detail.WINDOWS) / meanReadsPerWindow;
                            detail.ERROR_BAR_WIDTH = (Math.sqrt(detail.READ_STARTS) / (double) detail.WINDOWS) / meanReadsPerWindow;
                        } else {
                            detail.NORMALIZED_COVERAGE = 0;
                            detail.ERROR_BAR_WIDTH = 0;
                        }
                        detail.ACCUMULATION_LEVEL = group;
                        if (group.equals(ACCUMULATION_LEVEL_READ_GROUP)) {detail.READ_GROUP = gcType;}
                        else if (group.equals(ACCUMULATION_LEVEL_SAMPLE)) {detail.SAMPLE = gcType;}
                        else if (group.equals(ACCUMULATION_LEVEL_LIBRARY)) {detail.LIBRARY = gcType;}

                        detail.READS_USED = includeDuplicates ? READS_USED_ALL : READS_USED_UNIQUE;

                        metrics.DETAILS.addMetric(detail);
                    }

                    // Synthesize the high level summary metrics
                    final GcBiasSummaryMetrics summary = new GcBiasSummaryMetrics();
                    if (group.equals(ACCUMULATION_LEVEL_READ_GROUP)) {summary.READ_GROUP = gcType;}
                    else if (group.equals(ACCUMULATION_LEVEL_SAMPLE)) {summary.SAMPLE = gcType;}
                    else if (group.equals(ACCUMULATION_LEVEL_LIBRARY)) {summary.LIBRARY = gcType;}

                    summary.READS_USED = includeDuplicates ? READS_USED_ALL : READS_USED_UNIQUE;

                    summary.ACCUMULATION_LEVEL = group;
                    summary.WINDOW_SIZE = scanWindowSize;
                    summary.TOTAL_CLUSTERS = totalClusters;
                    summary.ALIGNED_READS = totalAlignedReads;
                    summary.GC_NC_0_19 = calculateGcNormCoverage(meanReadsPerWindow, readsByGc, 0, 19);
                    summary.GC_NC_20_39 = calculateGcNormCoverage(meanReadsPerWindow, readsByGc, 20, 39);
                    summary.GC_NC_40_59 = calculateGcNormCoverage(meanReadsPerWindow, readsByGc, 40, 59);
                    summary.GC_NC_60_79 = calculateGcNormCoverage(meanReadsPerWindow, readsByGc, 60, 79);
                    summary.GC_NC_80_100 = calculateGcNormCoverage(meanReadsPerWindow, readsByGc, 80, 100);

                    calculateDropoutMetrics(metrics.DETAILS.getMetrics(), summary);

                    metrics.SUMMARY = summary;

                    file.addMetric(metrics);
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculates the normalized coverage over a given gc content region
    /////////////////////////////////////////////////////////////////////////////
    private double calculateGcNormCoverage(final double meanReadsPerWindow, final int[] readsByGc,
                                           final int start, final int end) {
        int windowsTotal = 0;
        double sum = 0.0;
        for (int i = start; i <= end; i++) {
            if (windowsByGc[i] != 0) {
                sum += (double) readsByGc[i];
                windowsTotal += windowsByGc[i];
            }
        }

        if (windowsTotal == 0) {
            return 0.0;
        } else {
            return (sum / (windowsTotal * meanReadsPerWindow));
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculates the Illumina style AT and GC dropout numbers
    /////////////////////////////////////////////////////////////////////////////
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
                else{ gcDropout += dropout; }
            }
        }

        summary.AT_DROPOUT = atDropout;
        summary.GC_DROPOUT = gcDropout;
    }

    /////////////////////////////////////////////////////////////////////////////
    //Keeps track of each level of GcCalculation
    /////////////////////////////////////////////////////////////////////////////
    class GcObject {
        long totalClusters = 0;
        long totalAlignedReads = 0;
        int[] readsByGc = new int[BINS];
        long[] basesByGc = new long[BINS];
        long[] errorsByGc = new long[BINS];
        String group = null;
    }

    /////////////////////////////////////////////////////////////////////////////
    //Adds each read to the appropriate gcObj which is determined in acceptRecord above
    //Also calculates values for calculating GC Bias at each level
    /////////////////////////////////////////////////////////////////////////////
     private void addRead(final GcObject gcObj, final SAMRecord rec, final String group, final byte[] gc, final byte[] refBases) {
        if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++gcObj.totalClusters;
        final int pos = rec.getReadNegativeStrandFlag() ? rec.getAlignmentEnd() - scanWindowSize : rec.getAlignmentStart();
        ++gcObj.totalAlignedReads;
        if (pos > 0) {
            final int windowGc = gc[pos];
            if (windowGc >= 0) {
                ++gcObj.readsByGc[windowGc];
                gcObj.basesByGc[windowGc] += rec.getReadLength();
                gcObj.errorsByGc[windowGc] +=
                        SequenceUtil.countMismatches(rec, refBases, bisulfite) +
                                SequenceUtil.countInsertedBases(rec) + SequenceUtil.countDeletedBases(rec);
            }
        }
        if (gcObj.group == null) {
            gcObj.group = group;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Arguments that need to be passed to each PerUnitMetricCollector
// for the given record
/////////////////////////////////////////////////////////////////////////////
class GcBiasCollectorArgs {
    private final SAMRecord rec;
    private final ReferenceSequence ref;
    public SAMRecord getRec() {return rec;}
    public ReferenceSequence getRef() {return ref;}
    public GcBiasCollectorArgs(final SAMRecord rec, final ReferenceSequence ref) {
        this.rec = rec;
        this.ref = ref;
    }
}
