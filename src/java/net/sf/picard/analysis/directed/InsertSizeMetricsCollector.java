package net.sf.picard.analysis.directed;

import net.sf.picard.analysis.InsertSizeMetrics;
import net.sf.picard.analysis.MetricAccumulationLevel;
import net.sf.picard.metrics.MultiLevelCollector;
import net.sf.picard.metrics.PerUnitMetricCollector;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.sam.SamPairUtil;
import net.sf.picard.util.Histogram;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Collects InserSizeMetrics on the specified accumulationLevels using
 */
public class InsertSizeMetricsCollector extends MultiLevelCollector<InsertSizeMetrics, Integer, InsertSizeCollectorArgs> {
    // When generating the histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this
    // percentage of overall reads. (Range: 0 to 1)
    private final double minimumPct;

    // Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION.
    // This is done because insert size data typically includes enough anomolous values from chimeras and other
    // artifacts to make the mean and sd grossly misleading regarding the real distribution.
    private final double deviations;

    //Explicitly sets the histogram width, overriding automatic truncation of histogram tail.
    //Also, when calculating mean and stdev, only bins <= HISTOGRAM_WIDTH will be included.
    private Integer histogramWidth;

    public InsertSizeMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords,
                                      final double minimumPct, final Integer histogramWidth, final double deviations) {
        this.minimumPct = minimumPct;
        this.histogramWidth = histogramWidth;
        this.deviations = deviations;
        setup(accumulationLevels, samRgRecords);
    }

    // We will pass insertSize and PairOrientation with the DefaultPerRecordCollectorArgs passed to the record collectors
    // This method is called once Per samRecord
    @Override
    protected InsertSizeCollectorArgs makeArg(SAMRecord samRecord, ReferenceSequence refSeq) {
        final int insertSize = Math.abs(samRecord.getInferredInsertSize());
        final SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(samRecord);

        return new InsertSizeCollectorArgs(insertSize, orientation);
    }

    /** Make an InsertSizeCollector with the given arguments */
    @Override
    protected PerUnitMetricCollector<InsertSizeMetrics, Integer, InsertSizeCollectorArgs> makeChildCollector(final String sample, final String library, final String readGroup) {
        return new PerUnitInsertSizeMetricsCollector(sample, library, readGroup);
    }

    @Override
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
        if (!record.getReadPairedFlag() ||
                record.getReadUnmappedFlag() ||
                record.getMateUnmappedFlag() ||
                record.getFirstOfPairFlag() ||
                record.getNotPrimaryAlignmentFlag() ||
                record.getDuplicateReadFlag() ||
                record.getInferredInsertSize() == 0) {
            return;
        }

        super.acceptRecord(record, refSeq);
    }

    /** A Collector for individual InsertSizeMetrics for a given SAMPLE or SAMPLE/LIBRARY or SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels) */
    public class PerUnitInsertSizeMetricsCollector implements PerUnitMetricCollector<InsertSizeMetrics, Integer, InsertSizeCollectorArgs> {
        final EnumMap<SamPairUtil.PairOrientation, Histogram<Integer>> histograms = new EnumMap<SamPairUtil.PairOrientation, Histogram<Integer>>(SamPairUtil.PairOrientation.class);
        final String sample;
        final String library;
        final String readGroup;
        private double totalInserts = 0;

        public PerUnitInsertSizeMetricsCollector(final String sample, final String library, final String readGroup) {
            this.sample = sample;
            this.library = library;
            this.readGroup = readGroup;
            String prefix = null;
            if (this.readGroup != null) {
                prefix = this.readGroup + ".";
            }
            else if (this.library != null) {
                prefix = this.library + ".";
            }
            else if (this.sample != null) {
                prefix = this.sample + ".";
            }
            else {
                prefix = "All_Reads.";
            }
            histograms.put(SamPairUtil.PairOrientation.FR,     new Histogram<Integer>("insert_size", prefix + "fr_count"));
            histograms.put(SamPairUtil.PairOrientation.TANDEM, new Histogram<Integer>("insert_size", prefix + "tandem_count"));
            histograms.put(SamPairUtil.PairOrientation.RF,     new Histogram<Integer>("insert_size", prefix + "rf_count"));
        }

        public void acceptRecord(final InsertSizeCollectorArgs args) {
            histograms.get(args.getPairOrientation()).increment(args.getInsertSize());
        }

        public void finish() { }

        public double getTotalInserts() {
            return totalInserts;
        }

        public void addMetricsToFile(final MetricsFile<InsertSizeMetrics,Integer> file) {
            for (final Histogram<Integer> h : this.histograms.values()) totalInserts += h.getCount();

            for(final Map.Entry<SamPairUtil.PairOrientation, Histogram<Integer>> entry : histograms.entrySet()) {
                final SamPairUtil.PairOrientation pairOrientation = entry.getKey();
                final Histogram<Integer> histogram = entry.getValue();
                final double total = histogram.getCount();

                // Only include a category if it has a sufficient percentage of the data in it
                if( total > totalInserts * minimumPct ) {
                    final InsertSizeMetrics metrics = new InsertSizeMetrics();
                    metrics.SAMPLE             = this.sample;
                    metrics.LIBRARY            = this.library;
                    metrics.READ_GROUP         = this.readGroup;
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
                    if (histogramWidth == null) {
                        histogramWidth = (int) (metrics.MEDIAN_INSERT_SIZE + (deviations * metrics.MEDIAN_ABSOLUTE_DEVIATION));
                    }

                    trimmedHisto.trimByWidth(histogramWidth);

                    metrics.MEAN_INSERT_SIZE = trimmedHisto.getMean();
                    metrics.STANDARD_DEVIATION = trimmedHisto.getStandardDeviation();

                    file.addHistogram(trimmedHisto);
                    file.addMetric(metrics);
                }
            }
        }
    }
}

// Arguments that need to be calculated once per SAMRecord that are then passed to each PerUnitMetricCollector
// for the given record
class InsertSizeCollectorArgs {
    private final int insertSize;
    private final SamPairUtil.PairOrientation po;


    public int getInsertSize() {
        return insertSize;
    }

    public SamPairUtil.PairOrientation getPairOrientation() {
        return po;
    }

    public InsertSizeCollectorArgs(final int insertSize, final SamPairUtil.PairOrientation po) {
        this.insertSize = insertSize;
        this.po = po;
    }
}