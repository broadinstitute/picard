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
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.util.SeriesStats;
import picard.util.help.HelpConstants;

import java.io.File;
import java.util.Vector;

/**
 * Command line program to calculate SNV quality yield metrics for read files
 *
 * @author Dror Kessler
 */

@CommandLineProgramProperties(
        summary = CollectQualityYieldMetricsSNVQ.USAGE_SUMMARY + CollectQualityYieldMetricsSNVQ.USAGE_DETAILS,
        oneLineSummary = CollectQualityYieldMetricsSNVQ.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CollectQualityYieldMetricsSNVQ extends SinglePassSamProgram {
    public static final String SNVQ_BASE_ORDER = "ACGT";

    private QualityYieldMetricsCollector collector = null;
    public Histogram<Integer> qualityHistogram = new Histogram<>("KEY", "BQ_COUNT");
    private Vector<SeriesStats> readPositionQualityStats = new Vector<>();
    public Histogram<Integer> snvqHistogram = new Histogram<>("KEY", "SNVQ_COUNT");
    private Vector<SeriesStats> readPositionSnvqStats = new Vector<>();
    final byte[] baseOrder = SNVQ_BASE_ORDER.getBytes();

    static final String USAGE_SUMMARY = "Collect SNVQ metrics about reads that pass quality thresholds and other filters (such as vendor fail, etc).  ";
    static final String USAGE_DETAILS = "This tool evaluates the overall SNVQ quality of reads within a bam file containing one read group. " +
            "The output indicates the total numbers of bases within a read group that pass a minimum base quality score threshold and " +
            "(in the case of Illumina data) pass Illumina quality filters as described in the " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=6329'>GATK Dictionary entry</a>. " +
            "<br />" +

            "<h4>Usage Example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectSNVQualityYieldMetrics \\<br /> " +
            "      I=input.bam \\<br /> " +
            "      O=quality_yield_metrics.txt \\<br />" +
            "</pre>" +
            "Please see " +
            "<a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectSNVQualityYieldMetrics.QualityYieldMetrics'>" +
            "the QualityYieldMetrics documentation</a> for details and explanations of the output metrics." +
            "<hr />";

    @Argument(shortName = StandardOptionDefinitions.ALTERNATE_QUALITY_ATTRIBUTE_SHORT_NAME,
            doc = "Use an alternative  tag instead of base quality (QUAL) scores ", optional = true)
    public String ALTERNATE_QUALITY_ATTRIBUTE = null;

    @Argument(doc = "If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting " +
            "of bases if there are secondary alignments in the input file.")
    public boolean INCLUDE_SECONDARY_ALIGNMENTS = false;

    @Argument(doc = "If true, include bases from supplemental alignments in metrics. Setting to true may cause double-counting " +
            "of bases if there are supplemental alignments in the input file.")
    public boolean INCLUDE_SUPPLEMENTAL_ALIGNMENTS = false;

    @Argument(doc = "Determines whether to include the base quality histogram in the metrics file.")
    public boolean INCLUDE_BQ_HISTOGRAM = false;

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
        this.collector = new QualityYieldMetricsCollector(ALTERNATE_QUALITY_ATTRIBUTE, INCLUDE_SECONDARY_ALIGNMENTS, INCLUDE_SUPPLEMENTAL_ALIGNMENTS);
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        this.collector.acceptRecord(rec, ref);
    }

    @Override
    protected void finish() {
        final MetricsFile<QualityYieldMetrics, Integer> metricsFile = getMetricsFile();
        this.collector.finish();
        this.collector.addMetricsToFile(metricsFile);
        if ( INCLUDE_BQ_HISTOGRAM ) {
            metricsFile.addHistogram(qualityHistogram);
            metricsFile.addHistogram(snvqHistogram);
            this.collector.addHistograms(metricsFile);
        }
        metricsFile.write(OUTPUT);
    }

    public class QualityYieldMetricsCollector {
        // optional attribute name to retrieve quality from (instead of standard QUAL string)
        private final String alternateQualityAttribute;

        // If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting
        // of bases if there are secondary alignments in the input file.
        private final boolean includeSecondaryAlignments;

        // If true, include bases from supplemental alignments in metrics. Setting to true may cause double-counting
        // of bases if there are supplemental alignments in the input file.
        public final boolean includeSupplementalAlignments;

        // The metrics to be accumulated
        private final QualityYieldMetrics metrics;

        public QualityYieldMetricsCollector(final String alternateQualityAttribute,
                                            final boolean includeSecondaryAlignments,
                                            final boolean includeSupplementalAlignments) {
            this.alternateQualityAttribute = alternateQualityAttribute;
            this.includeSecondaryAlignments = includeSecondaryAlignments;
            this.includeSupplementalAlignments = includeSupplementalAlignments;
            this.metrics = new QualityYieldMetrics(alternateQualityAttribute);
        }

        public void acceptRecord(final SAMRecord rec, final ReferenceSequence ref) {
            if (!this.includeSecondaryAlignments && rec.isSecondaryAlignment()) return;
            if (!this.includeSupplementalAlignments && rec.getSupplementaryAlignmentFlag()) return;

            final int length = rec.getReadLength();
            metrics.TOTAL_READS++;
            metrics.TOTAL_BASES += length;

            final boolean isPfRead = !rec.getReadFailsVendorQualityCheckFlag();
            if (isPfRead) {
                metrics.PF_READS++;
                metrics.PF_BASES += length;
            }

            // access regular quality
            final byte[] quals;
            if (alternateQualityAttribute != null) {
                byte[] tmp = rec.getStringAttribute(alternateQualityAttribute).getBytes();
                if (tmp == null) {
                    // fall back on base queslities
                    tmp = rec.getBaseQualities();
                } else {
                    SAMUtils.fastqToPhred(tmp);
                }
                quals = tmp;
            } else {
                quals = rec.getBaseQualities();
            }
            if ( quals.length != rec.getReadBases().length ) {
                throw new PicardException("quality string length does not match bases string");
            }

            // access snv qualities
            final byte[][] snvq = new byte[4][];
            for ( int i = 0 ; i < snvq.length ; i++ )
                snvq[i] = rec.getStringAttribute(String.format("q%c", Character.toLowerCase(baseOrder[i]))).getBytes();

            // access bases
            final byte[] bases = rec.getReadBases();

            // add up quals, and quals >= 20
            int readPosition = 0;
            for (final int qual : quals) {

                if (qual >= 40) {
                    metrics.Q20_BASES++;
                    metrics.Q30_BASES++;
                    metrics.Q40_BASES++;
                } else if (qual >= 30) {
                    metrics.Q20_BASES++;
                    metrics.Q30_BASES++;
                } else if (qual >= 20) {
                    metrics.Q20_BASES++;
                }

                if (isPfRead) {
                    if (qual >= 40) {
                        metrics.PF_Q20_BASES++;
                        metrics.PF_Q30_BASES++;
                        metrics.PF_Q40_BASES++;
                    } else if (qual >= 30) {
                        metrics.PF_Q20_BASES++;
                        metrics.PF_Q30_BASES++;
                    } else if (qual >= 20) {
                        metrics.PF_Q20_BASES++;
                    }
                }

                // SNVQ counters
                final byte base = bases[readPosition];
                for ( int i = 0 ; i < baseOrder.length ; i++ ) {
                    if ( base != baseOrder[i] ) {

                        int q = SAMUtils.fastqToPhred((char)snvq[i][readPosition]);
                        metrics.TOTAL_SNVQ++;
                        if ( isPfRead )
                            metrics.PF_SNVQ++;

                        if (q >= 40) {
                            metrics.Q20_SNVQ++;
                            metrics.Q30_SNVQ++;
                            metrics.Q40_SNVQ++;
                            if ( isPfRead ) {
                                metrics.PF_Q20_SNVQ++;
                                metrics.PF_Q30_SNVQ++;
                                metrics.PF_Q40_SNVQ++;
                            }
                        } else if (q >= 30) {
                            metrics.Q20_SNVQ++;
                            metrics.Q30_SNVQ++;
                            if ( isPfRead ) {
                                metrics.PF_Q20_SNVQ++;
                                metrics.PF_Q30_SNVQ++;
                            }
                        } else if (q >= 20) {
                            metrics.Q20_SNVQ++;
                            if ( isPfRead ) {
                                metrics.PF_Q20_SNVQ++;
                            }
                        }

                        if (INCLUDE_BQ_HISTOGRAM) {
                            snvqHistogram.increment(q);
                            while (readPositionSnvqStats.size() <= readPosition)
                                readPositionSnvqStats.add(new SeriesStats());
                            readPositionSnvqStats.get(readPosition).add(q);
                        }
                    }
                }

                if (INCLUDE_BQ_HISTOGRAM) {

                    // enter quality into histograms
                    qualityHistogram.increment(qual);

                    // collect read position quality stats
                    while (readPositionQualityStats.size() <= readPosition)
                        readPositionQualityStats.add(new SeriesStats());
                    readPositionQualityStats.get(readPosition).add(qual);
                }

                readPosition++;
            }
        }

        public void finish() {
            metrics.calculateDerivedFields();
        }

        public void addMetricsToFile(final MetricsFile<QualityYieldMetrics, Integer> metricsFile) {
            metricsFile.addMetric(metrics);
        }

        public void addHistograms(MetricsFile<QualityYieldMetrics, Integer> metricsFile) {

            final Histogram<Integer> h1 = new Histogram<>("KEY", "READ_INDEX_MEAN_BQ");
            for ( int i = 0; i < readPositionQualityStats.size() ; i++ ) {
                SeriesStats ss = readPositionQualityStats.get(i);
                h1.increment(i, ss.getMean());
            }
            metricsFile.addHistogram(h1);

            final Histogram<Integer> h2 = new Histogram<>("KEY", "READ_INDEX_MEAN_SNVQ");
            for ( int i = 0; i < readPositionSnvqStats.size() ; i++ ) {
                SeriesStats ss = readPositionSnvqStats.get(i);
                h2.increment(i, ss.getMean());
            }
            metricsFile.addHistogram(h2);
        }
    }

    /**
     * A set of metrics used to describe the general quality of a BAM file
     */
    @DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
    public static class QualityYieldMetrics extends MergeableMetricBase {

        public QualityYieldMetrics() {
            this(null);
        }

        public QualityYieldMetrics(final String alternateQualityAttribute) {
            super();
            this.alternateQualityAttribute = alternateQualityAttribute;
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
         * The average read length of all the reads
         */
        @NoMergingIsDerived
        public int READ_LENGTH = 0;

        /**
         * The total number of bases in all reads
         */
        @MergeByAdding
        public long TOTAL_BASES = 0;

        /**
         * The total number of bases in all PF reads
         */
        @MergeByAdding
        public long PF_BASES = 0;

        /**
         * The number of bases in all reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long Q20_BASES = 0;

        /**
         * The number of bases in PF reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long PF_Q20_BASES = 0;

        /**
         * The number of bases in all reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long Q30_BASES = 0;

        /**
         * The number of bases in PF reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long PF_Q30_BASES = 0;

        /**
         * The number of bases in all reads that achieve quality score 40 or higher
         */
        @MergeByAdding
        public long Q40_BASES = 0;

        /**
         * The number of bases in PF reads that achieve quality score 40 or higher
         */
        @MergeByAdding
        public long PF_Q40_BASES = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 20 or higher
         */
        @MergingIsManual
        public double PCT_Q20_BASES = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 30 or higher
         */
        @MergingIsManual
        public double PCT_Q30_BASES = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 40 or higher
         */
        @MergingIsManual
        public double PCT_Q40_BASES = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 20 or higher and pass filter
         */
        @MergingIsManual
        public double PCT_PF_Q20_BASES = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 30 or higher and pass filter
         */
        @MergingIsManual
        public double PCT_PF_Q30_BASES = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 40 or higher and pass filter
         */
        @MergingIsManual
        public double PCT_PF_Q40_BASES = 0;

        /**
         * The total number of SNVQ values in all reads
         */
        @MergeByAdding
        public long TOTAL_SNVQ = 0;

        /**
         * The total number of SNVQ values in all PF reads
         */
        @MergeByAdding
        public long PF_SNVQ = 0;

        /**
         * The number of SNVQ values in all reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long Q20_SNVQ = 0;

        /**
         * The number of SNVQ values in PF reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long PF_Q20_SNVQ = 0;

        /**
         * The number of SNVQ values in all reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long Q30_SNVQ = 0;

        /**
         * The number of SNVQ values in PF reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long PF_Q30_SNVQ = 0;

        /**
         * The number of SNVQ values in all reads that achieve quality score 40 or higher
         */
        @MergeByAdding
        public long Q40_SNVQ = 0;

        /**
         * The number of SNVQ values in PF reads that achieve quality score 40 or higher
         */
        @MergeByAdding
        public long PF_Q40_SNVQ = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 20 or higher
         */
        @MergingIsManual
        public double PCT_Q20_SNVQ = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 30 or higher
         */
        @MergingIsManual
        public double PCT_Q30_SNVQ = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 40 or higher
         */
        @MergingIsManual
        public double PCT_Q40_SNVQ = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 20 or higher and pass filter
         */
        @MergingIsManual
        public double PCT_PF_Q20_SNVQ = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 30 or higher and pass filter
         */
        @MergingIsManual
        public double PCT_PF_Q30_SNVQ = 0;

        /**
         * The percentage of SNVQ values in all reads that achieve quality score 40 or higher and pass filter
         */
        @MergingIsManual
        public double PCT_PF_Q40_SNVQ = 0;

        @MergeByAssertEquals
        protected final String  alternateQualityAttribute;

        @Override
        public void calculateDerivedFields() {
            super.calculateDerivedFields();
            this.READ_LENGTH = this.TOTAL_READS == 0 ? 0 : (int) (this.TOTAL_BASES / this.TOTAL_READS);

            this.PCT_Q20_BASES = this.TOTAL_BASES == 0 ? 0 : (double)this.Q20_BASES / this.TOTAL_BASES;
            this.PCT_Q30_BASES = this.TOTAL_BASES == 0 ? 0 : (double)this.Q30_BASES / this.TOTAL_BASES;
            this.PCT_Q40_BASES = this.TOTAL_BASES == 0 ? 0 : (double)this.Q40_BASES / this.TOTAL_BASES;

            this.PCT_PF_Q20_BASES = this.PF_BASES == 0 ? 0 : (double)this.PF_Q20_BASES / this.PF_BASES;
            this.PCT_PF_Q30_BASES = this.PF_BASES == 0 ? 0 : (double)this.PF_Q30_BASES / this.PF_BASES;
            this.PCT_PF_Q40_BASES = this.PF_BASES == 0 ? 0 : (double)this.PF_Q40_BASES / this.PF_BASES;


            this.PCT_Q20_SNVQ = this.TOTAL_SNVQ == 0 ? 0 : (double)this.Q20_SNVQ / this.TOTAL_SNVQ;
            this.PCT_Q30_SNVQ = this.TOTAL_SNVQ == 0 ? 0 : (double)this.Q30_SNVQ / this.TOTAL_SNVQ;
            this.PCT_Q40_SNVQ = this.TOTAL_SNVQ == 0 ? 0 : (double)this.Q40_SNVQ / this.TOTAL_SNVQ;

            this.PCT_PF_Q20_SNVQ = this.PF_SNVQ == 0 ? 0 : (double)this.PF_Q20_SNVQ / this.PF_SNVQ;
            this.PCT_PF_Q30_SNVQ = this.PF_SNVQ == 0 ? 0 : (double)this.PF_Q30_SNVQ / this.PF_SNVQ;
            this.PCT_PF_Q40_SNVQ = this.PF_SNVQ == 0 ? 0 : (double)this.PF_Q40_SNVQ / this.PF_SNVQ;
        }

        @Override
        public MergeableMetricBase merge(final MergeableMetricBase other) {
            if (!(other instanceof QualityYieldMetrics)){
                throw new PicardException("Only objects of the same type can be merged");
            }

            final QualityYieldMetrics otherMetric = (QualityYieldMetrics) other;

            super.merge(otherMetric);
            calculateDerivedFields();
            return this;
        }
    }
}
