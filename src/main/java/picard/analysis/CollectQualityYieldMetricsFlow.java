/*
 * The MIT License
 *
 * Copyright (c) 2024 The Broad Institute
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
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.flow.FlowBasedArgumentCollection;
import picard.flow.FlowBasedRead;
import picard.flow.FlowBasedKeyCodec;
import picard.flow.FlowReadGroupInfo;
import picard.util.SeriesStats;
import picard.util.help.HelpConstants;

import java.io.File;
import java.util.Vector;

/**
 * Command line program to calculate quality yield metrics for flow based read files
 *
 * @author Dror Kessler
 */


@CommandLineProgramProperties(
        summary = CollectQualityYieldMetricsFlow.USAGE_SUMMARY + CollectQualityYieldMetricsFlow.USAGE_DETAILS,
        oneLineSummary = CollectQualityYieldMetricsFlow.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@ExperimentalFeature
public class CollectQualityYieldMetricsFlow extends SinglePassSamProgram {

    private static final byte MIN_QUAL = 0;
    private static final byte MAX_QUAL = 100;
    private final Log log = Log.getInstance(CollectQualityYieldMetricsFlow.class);

    private static final int CYCLE_LENGTH = 4;
    private QualityYieldMetricsCollectorFlow collector = null;
    public Histogram<Integer> qualityHistogram = new Histogram<>("KEY", "QUAL_COUNT");
    private Vector<SeriesStats> flowQualityStats = new Vector<>();

    static final String USAGE_SUMMARY = "Collect metrics about reads that pass quality thresholds from flow based read files.  ";
    static final String USAGE_DETAILS = "This tool evaluates the overall quality of reads within a bam file containing one read group. " +
            "The output indicates the total numbers of flows within a read group that pass a minimum base quality score threshold " +
            "<h4>Usage Example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectQualityYieldMetricsFlow \\<br /> " +
            "      I=input.bam \\<br /> " +
            "      O=quality_yield_metrics.txt \\<br />" +
            "</pre>" +
            "<hr />";

    @Argument(doc = "If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting " +
            "of bases if there are secondary alignments in the input file.")
    public boolean INCLUDE_SECONDARY_ALIGNMENTS = false;

    @Argument(doc = "If true, include bases from supplemental alignments in metrics. Setting to true may cause double-counting " +
            "of bases if there are supplemental alignments in the input file.")
    public boolean INCLUDE_SUPPLEMENTAL_ALIGNMENTS = false;

    @Argument(doc = "Determines whether to include the flow quality histogram in the metrics file.")
    public boolean INCLUDE_BQ_HISTOGRAM = false;

    @ArgumentCollection(doc = "flow based args")
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    /**
     * Ensure that we get all reads regardless of alignment status.
     */
    @Override
    protected boolean usesNoRefReads() {
        return true;
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        this.collector = new QualityYieldMetricsCollectorFlow(INCLUDE_SECONDARY_ALIGNMENTS, INCLUDE_SUPPLEMENTAL_ALIGNMENTS, fbargs);
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        this.collector.acceptRecord(rec, ref);
    }

    @Override
    protected void finish() {
        final MetricsFile<QualityYieldMetricsFlow, Integer> metricsFile = getMetricsFile();
        this.collector.finish();
        this.collector.addMetricsToFile(metricsFile);
        if ( INCLUDE_BQ_HISTOGRAM ) {
            metricsFile.addHistogram(qualityHistogram);
            this.collector.addHistograms(metricsFile);
        }
        metricsFile.write(OUTPUT);
    }

    public class QualityYieldMetricsCollectorFlow {
        // If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting
        // of bases if there are secondary alignments in the input file.
        private final boolean includeSecondaryAlignments;

        // If true, include bases from supplemental alignments in metrics. Setting to true may cause double-counting
        // of bases if there are supplemental alignments in the input file.
        public final boolean includeSupplementalAlignments;

        // The metrics to be accumulated
        private final QualityYieldMetricsFlow metrics;

        // flow based args
        private final FlowBasedArgumentCollection fbargs;

        public QualityYieldMetricsCollectorFlow(final boolean includeSecondaryAlignments,
                                                final boolean includeSupplementalAlignments,
                                                final FlowBasedArgumentCollection fbargs) {
            this.includeSecondaryAlignments = includeSecondaryAlignments;
            this.includeSupplementalAlignments = includeSupplementalAlignments;
            this.fbargs = fbargs;
            this.metrics = new QualityYieldMetricsFlow();
        }

        public void acceptRecord(final SAMRecord rec, final ReferenceSequence ref) {
            if (rec.getReadLength() == 0) return;
            if (!this.includeSecondaryAlignments && rec.isSecondaryAlignment()) return;
            if (!this.includeSupplementalAlignments && rec.getSupplementaryAlignmentFlag()) return;
            metrics.TOTAL_READS++;

            final boolean isPfRead = !rec.getReadFailsVendorQualityCheckFlag();
            if ( !isPfRead ) {
                return;
            }

            // NOTE: code below runs only isPfRead reads

            // convert to a flow based read
            FlowReadGroupInfo info = FlowBasedKeyCodec.getReadGroupInfo(rec.getHeader(), rec);
            if ( !info.isFlowPlatform ) {
                throw new PicardException("Reads should originate from a flow based platform");
            }
            FlowBasedRead fread = new FlowBasedRead(rec, info.flowOrder, info.maxClass, fbargs);
            metrics.PF_READS++;
            metrics.PF_FLOWS += fread.getKey().length;

            // get flow quals
            final byte[] quals = getFlowQualities(fread);

            // allocate cycle qual accounting
            int cycleCount = INCLUDE_BQ_HISTOGRAM ? (int)Math.ceil((float)quals.length / CYCLE_LENGTH) : 0;
            int cycleQualCount[] = INCLUDE_BQ_HISTOGRAM ? new int[cycleCount] : null;
            int cycleQualSum[] = INCLUDE_BQ_HISTOGRAM ? new int[cycleCount] : null;

            // add up quals, and quals >= 20
            int flow = 0;
            for (final int qual : quals) {
                metrics.PF_Q20_EQUIVALENT_YIELD += qual;
                if (qual >= 30) {
                    metrics.PF_Q20_FLOWS++;
                    metrics.PF_Q30_FLOWS++;
                } else if (qual >= 20) {
                    metrics.PF_Q20_FLOWS++;
                }

                if ( INCLUDE_BQ_HISTOGRAM ) {

                    // enter quality into histograms
                    qualityHistogram.increment(qual);

                    // enter quality into cycle stats
                    final int cycle = flow / CYCLE_LENGTH;
                    cycleQualCount[cycle]++;
                    cycleQualSum[cycle] += qual;
                }

                // advance
                flow++;
            }

            // make sure flowQualityCount/Sum are large enough and enter accounted values
            if ( INCLUDE_BQ_HISTOGRAM ) {
                while (flowQualityStats.size() < cycleCount)
                    flowQualityStats.add(new SeriesStats());
                for (int cycle = 0; cycle < cycleCount; cycle++) {
                    int id = !rec.getReadNegativeStrandFlag() ? cycle : (cycleCount - 1 - cycle);
                    flowQualityStats.get(id).add((double) cycleQualSum[cycle] / cycleQualCount[cycle]);
                }
            }
        }

        private byte[] getFlowQualities(FlowBasedRead fread) {

            double[] errorProbs = computeErrorProb(fread);
            byte[] quals = new byte[errorProbs.length];
            for ( int i = 0 ; i < errorProbs.length ; i++ ) {
                if ( errorProbs[i] == 0.0 ) {
                    // this is a special case that should not happen
                    log.warn(fread.getReadName() + ": zero errorProb on flow: " + i);
                    quals[i] = MAX_QUAL;
                } else {
                    long q = Math.round(-10 * Math.log10(errorProbs[i]));
                    if ( q < MIN_QUAL || q > MAX_QUAL ) {
                        // this is an out-of-range condition. should not happen as well
                        log.warn(fread.getReadName() + ": qual " + q + " is out of range on flow: " + i);
                        q = Math.max(MIN_QUAL, Math.min(MAX_QUAL, q));
                    }
                    quals[i] = (byte)q;
                }
            }
            return quals;
        }

        /*
         * compute error probability vector for a read
         *
         * The vector has one element for each flow key, representing the probability complementing the call-probability to 1
         * This is further complicated by the optional presence of a genome-prior database, which provides factoring for
         * each hmer length (on a base basis)
         */
        private double[] computeErrorProb(final FlowBasedRead flowRead) {

            final int[] key = flowRead.getKey();
            final double[] probCol = new double[flowRead.getMaxHmer() + 1];
            double[] result = new double[key.length];

            for ( int i = 0 ; i < key.length ; i++ ) {

                // step 1 - extract column & sum
                double  sum = 0;
                for ( int j = 0 ; j < probCol.length ; j++ ) {
                    sum += (probCol[j] = flowRead.getProb(i, j));
                }

                // step 2 - normalize column
                for ( int j = 0 ; j < probCol.length ; j++ ) {
                    probCol[j] /= sum;
                }

                // assign normalized result
                result[i] = 1 - probCol[Math.min(key[i], flowRead.getMaxHmer())];
            }

            return result;
        }

        public void finish() {
            metrics.PF_Q20_EQUIVALENT_YIELD = metrics.PF_Q20_EQUIVALENT_YIELD / 20;
            metrics.calculateDerivedFields();
        }

        public void addMetricsToFile(final MetricsFile<QualityYieldMetricsFlow, Integer> metricsFile) {
            metricsFile.addMetric(metrics);
        }

        public void addHistograms(MetricsFile<QualityYieldMetricsFlow, Integer> metricsFile) {

            Histogram<Integer> meanHist = new Histogram<>("KEY", "MEAN_CYCLE_QUAL");
            Histogram<Integer> medianHist = new Histogram<>("KEY", "MEDIAN_CYCLE_QUAL");
            Histogram<Integer> q25Hist = new Histogram<>("KEY", "Q25_CYCLE_QUAL");
            Histogram<Integer> q75Hist = new Histogram<>("KEY", "Q75_CYCLE_QUAL");
            for (int i = 0; i < flowQualityStats.size() ; i++ ) {
                SeriesStats ss = flowQualityStats.get(i);
                meanHist.increment(i, ss.getMean());
                medianHist.increment(i, ss.getMedian());
                q25Hist.increment(i, ss.getPercentile(25));
                q75Hist.increment(i, ss.getPercentile(75));
            }
            metricsFile.addHistogram(meanHist);
            metricsFile.addHistogram(medianHist);
            metricsFile.addHistogram(q25Hist);
            metricsFile.addHistogram(q75Hist);
        }
    }

    /**
     * A set of metrics used to describe the general quality of a BAM file
     */
    @DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
    public static class QualityYieldMetricsFlow extends MergeableMetricBase {

        public QualityYieldMetricsFlow() {
        }

        /**
         * The total number of reads in the input file
         */
        @MergeByAdding
        public long TOTAL_READS = 0;

        /**
         * The number of reads that are PF - pass filter
         */
        @MergeByAdding
        public long PF_READS = 0;

        /**
         * The average number of flows in PF reads
         */
        @NoMergingIsDerived
        public int MEAN_PF_READ_NUMBER_OF_FLOWS = 0;

        /**
         * The total number of flows in all PF reads
         */
        @MergeByAdding
        public long PF_FLOWS = 0;

        /**
         * The number of flows in PF reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long PF_Q20_FLOWS = 0;

        /**
         * The percentage of flows in all reads that achieve quality score 20 or higher
         */
        @MergingIsManual
        public double PCT_PF_Q20_FLOWS = 0;

        /**
         * The number of flows in PF reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long PF_Q30_FLOWS = 0;

        /**
         * The percentage of flows in all reads that achieve quality score 30 or higher
         */
        @MergingIsManual
        public double PCT_PF_Q30_FLOWS = 0;

        /**
         * The sum of quality scores of all flows in PF reads divided by 20
         */
        @MergeByAdding
        public long PF_Q20_EQUIVALENT_YIELD = 0;

        @Override
        public void calculateDerivedFields() {
            super.calculateDerivedFields();
            this.MEAN_PF_READ_NUMBER_OF_FLOWS = this.PF_READS == 0 ? 0 : (int) (this.PF_FLOWS / this.PF_READS);
            this.PCT_PF_Q20_FLOWS = this.PF_FLOWS == 0 ? 0 : (double)this.PF_Q20_FLOWS / this.PF_FLOWS;
            this.PCT_PF_Q30_FLOWS = this.PF_FLOWS == 0 ? 0 : (double)this.PF_Q30_FLOWS / this.PF_FLOWS;
        }

        @Override
        public MergeableMetricBase merge(final MergeableMetricBase other) {
            if (!(other instanceof QualityYieldMetricsFlow)){
                throw new PicardException("Only objects of the same type can be merged");
            }

            final QualityYieldMetricsFlow otherMetric = (QualityYieldMetricsFlow) other;

            super.merge(otherMetric);
            calculateDerivedFields();
            return this;
        }
    }
}
