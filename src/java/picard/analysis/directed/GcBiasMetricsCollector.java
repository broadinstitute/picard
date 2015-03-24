package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SequenceUtil;
import picard.analysis.GcBiasDetailMetrics;
import picard.analysis.GcBiasSummaryMetrics;
import picard.metrics.GcBiasMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.metrics.MultiLevelCollector;
import picard.metrics.PerUnitMetricCollector;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashMap;

/** Calculates GC Bias Metrics on multiple levels
 * Created by kbergin on 3/23/15.
 */
public class GcBiasMetricsCollector extends MultiLevelCollector<GcBiasMetrics, Integer, GcBiasCollectorArgs> {
    // Histograms to track the number of windows at each GC, and the number of read starts
    // at windows of each GC. Need 101 to get from 0-100.
    public final int windowSize;
    public final boolean bisulfite;
    public Map<String, byte[]> gcByRef;
    public int extraClusters = 0;
    public int[] windowsByGc = new int[WINDOWS];
    public static final int WINDOWS = 101;

    public GcBiasMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final Map<String, byte[]> gcByRef, final int[] windowsByGc, final List<SAMReadGroupRecord> samRgRecords, final int windowSize, final boolean bisulfite) {
        this.windowSize = windowSize;
        this.bisulfite = bisulfite;
        this.gcByRef = gcByRef;
        this.windowsByGc = windowsByGc;
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

    @Override
    public void acceptRecord(final SAMRecord rec, final ReferenceSequence ref) {super.acceptRecord(rec, ref);}

    /////////////////////////////////////////////////////////////////////////////
    //A collector for individual GcBiasMetrics for a given SAMPLE or SAMPLE/LIBRARY
    //or SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels)
    /////////////////////////////////////////////////////////////////////////////
    public class PerUnitGcBiasMetricsCollector implements PerUnitMetricCollector<GcBiasMetrics, Integer, GcBiasCollectorArgs> {
        Map<String, GcObject> gcData = new HashMap<String, GcObject>();
        String sample = null;
        String library = null;
        String readGroup = null;

        /////////////////////////////////////////////////////////////////////////////
        //Records the accumulation level for each level of collection and initializes
        // a GcObject for this accumulation level
        /////////////////////////////////////////////////////////////////////////////
        public PerUnitGcBiasMetricsCollector(final String sample, final String library, final String readGroup) {
            this.sample = sample;
            this.library = library;
            this.readGroup = readGroup;
            final String prefix;
            if (this.readGroup != null) {
                prefix = this.readGroup;
                gcData.put(prefix, new GcObject());
            } else if (this.library != null) {
                prefix = this.library;
                gcData.put(prefix, new GcObject());
            } else if (this.sample != null) {
                prefix = this.sample;
                gcData.put(prefix, new GcObject());
            } else {
                prefix = "All_Reads";
                gcData.put(prefix, new GcObject());
            }
        }

        /////////////////////////////////////////////////////////////////////////////
        //Takes each record and sends them to addRead to calculate gc metrics for
        // that read for each accumulation level
        /////////////////////////////////////////////////////////////////////////////
        public void acceptRecord(final GcBiasCollectorArgs args) {
            final SAMRecord rec = args.getRec();
            final String type;
            if(!rec.getReadUnmappedFlag()) {
                final ReferenceSequence ref = args.getRef();
                final byte[] refBases = ref.getBases();
                final String refName = ref.getName();
                final byte[] gc = gcByRef.get(refName);
                final String group;
                if (this.readGroup != null) {
                    type = this.readGroup;
                    group = "Read Group";
                    addRead(gcData.get(type), rec, group, gc, refBases);
                } else if (this.library != null) {
                    type = this.library;
                    group = "Library";
                    addRead(gcData.get(type), rec, group, gc, refBases);
                } else if (this.sample != null) {
                    type = this.sample;
                    group = "Sample";
                    addRead(gcData.get(type), rec, group, gc, refBases);
                } else {
                    type = "All_Reads";
                    group = "All Reads";
                    addRead(gcData.get(type), rec, group, gc, refBases);
                }
            }
            else {
                for (final Map.Entry<String, GcObject> entry : gcData.entrySet()) {
                    final GcObject gcCur = entry.getValue();
                    if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++gcCur.totalClusters;
                }
            }
        }

        public void finish() {}

        /////////////////////////////////////////////////////////////////////////////
        // Sums the values in an int[].
        /////////////////////////////////////////////////////////////////////////////
        private double sum(final int[] values) {
            final int length = values.length;
            double total = 0;
            for (int i = 0; i < length; ++i) {
                total += values[i];
            }

            return total;
        }

        /////////////////////////////////////////////////////////////////////////////
        //Called to add metrics to the output file for each level of collection
        // these metrics are used for graphing gc bias in R script
        /////////////////////////////////////////////////////////////////////////////
        public void addMetricsToFile(final MetricsFile<GcBiasMetrics, Integer> file) {
            for (final Map.Entry<String, GcObject> entry : gcData.entrySet()) {
                final GcObject gcCur = entry.getValue();
                final String gcType = entry.getKey();
                final int[] readsByGc = gcCur.readsByGc;
                final long[] errorsByGc = gcCur.errorsByGc;
                final long[] basesByGc = gcCur.basesByGc;
                final int totalClusters = gcCur.totalClusters;
                final int totalAlignedReads = gcCur.totalAlignedReads;
                final String group = gcCur.group;
                final GcBiasMetrics metrics = new GcBiasMetrics();
                final double totalWindows = sum(windowsByGc);
                final double totalReads = sum(readsByGc);
                final double meanReadsPerWindow = totalReads / totalWindows;
                if(totalAlignedReads>0) {
                    for (int i = 0; i < windowsByGc.length; ++i) {
                        final GcBiasDetailMetrics m = new GcBiasDetailMetrics();
                        m.GC = i;
                        m.WINDOWS = windowsByGc[i];
                        m.READ_STARTS = readsByGc[i];
                        if (errorsByGc[i] > 0) m.MEAN_BASE_QUALITY = QualityUtil.getPhredScoreFromObsAndErrors(basesByGc[i], errorsByGc[i]);
                        if(windowsByGc[i]!=0) {
                            m.NORMALIZED_COVERAGE = (m.READ_STARTS / (double) m.WINDOWS) / meanReadsPerWindow;
                            m.ERROR_BAR_WIDTH = (Math.sqrt(m.READ_STARTS) / (double) m.WINDOWS) / meanReadsPerWindow;
                        }
                        else{
                            m.NORMALIZED_COVERAGE = 0;
                            m.ERROR_BAR_WIDTH = 0;
                        }
                        m.ACCUMULATION_LEVEL = group;
                        if(group.equals("Read Group")) {m.READ_GROUP = gcType;}
                        else if(group.equals("Sample")) {m.SAMPLE = gcType;}
                        else if(group.equals("Library")) {m.LIBRARY = gcType;}
                        metrics.DETAILS.addMetric(m);
                    }
                    // Synthesize the high level summary metrics
                    final GcBiasSummaryMetrics s = new GcBiasSummaryMetrics();
                    if(group.equals("Read Group")) {s.READ_GROUP = gcType;}
                    else if(group.equals("Sample")) {s.SAMPLE = gcType;}
                    else if(group.equals("Library")) {s.LIBRARY = gcType;}
                    s.ACCUMULATION_LEVEL = group;
                    s.WINDOW_SIZE = windowSize;
                    s.TOTAL_CLUSTERS = totalClusters;
                    s.ALIGNED_READS = totalAlignedReads;
                    calculateDropoutMetrics(metrics.DETAILS.getMetrics(), s);
                    metrics.SUMMARY = s;
                    file.addMetric(metrics);
                }
            }
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
                if (detail.GC >= 50) gcDropout += dropout;
            }
        }

        summary.AT_DROPOUT = atDropout;
        summary.GC_DROPOUT = gcDropout;
    }

    /////////////////////////////////////////////////////////////////////////////
    //Keeps track of each level of GcCalculation
    /////////////////////////////////////////////////////////////////////////////
    class GcObject{
        int totalClusters = 0;
        int totalAlignedReads = 0;
        int[] readsByGc = new int[WINDOWS];
        long[] basesByGc = new long[WINDOWS];
        long[] errorsByGc = new long[WINDOWS];
        String group = null;
    }

    /////////////////////////////////////////////////////////////////////////////
    //Adds each read to the appropriate gcObj which is determined in acceptRecord above
    //Also calculates values for calculating GC Bias at each level
    /////////////////////////////////////////////////////////////////////////////
     private void addRead(final GcObject gcObj, final SAMRecord rec, final String group, final byte[] gc, final byte[] refBases) {
        if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++gcObj.totalClusters;
        final int pos = rec.getReadNegativeStrandFlag() ? rec.getAlignmentEnd() - windowSize : rec.getAlignmentStart();
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
