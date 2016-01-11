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
package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.DbSnpBitSetUtil;
import picard.vcf.processor.VariantProcessor;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/** Collects summary and per-sample metrics about variant calls in a VCF file. */
@CommandLineProgramProperties(
        usage = "Collects per-sample and aggregate (spanning all samples) metrics from the provided VCF file.",
        usageShort = "Collects per-sample and aggregate (spanning all samples) metrics from the provided VCF file",
        programGroup = Metrics.class
)
public class CollectVariantCallingMetrics extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input vcf file for analysis")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Path (except for the file extension) of output metrics files " +
            "to write.")
    public File OUTPUT;

    @Option(doc = "Reference dbSNP file in dbSNP or VCF format.")
    public File DBSNP;

    @Option(shortName = "TI", doc = "Target intervals to restrict analysis to.", optional = true)
    public File TARGET_INTERVALS;

    @Option(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, optional = true,
            doc = "If present, speeds loading of dbSNP file, will look for dictionary in vcf if not present here.")
    public File SEQUENCE_DICTIONARY = null;

    @Option(doc = "Set to true if running on a single-sample gvcf.", optional = true)
    public boolean GVCF_INPUT = false;

    @Option
    public int THREAD_COUNT = 1;

    private final Log log = Log.getInstance(CollectVariantCallingMetrics.class);

    public static void main(final String[] args) {
        new CollectVariantCallingMetrics().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);

        final boolean requiresIndex = this.TARGET_INTERVALS != null || this.THREAD_COUNT > 1;
        final VCFFileReader variantReader = new VCFFileReader(INPUT, requiresIndex);
        final VCFHeader vcfHeader = variantReader.getFileHeader();
        CloserUtil.close(variantReader);

        final SAMSequenceDictionary sequenceDictionary =
                SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY == null ? INPUT : SEQUENCE_DICTIONARY);

        log.info("Loading dbSNP file ...");
        final DbSnpBitSetUtil.DbSnpBitSets dbsnp = DbSnpBitSetUtil.createSnpAndIndelBitSets(DBSNP, sequenceDictionary);

        log.info("Starting iteration of variants.");

        final VariantProcessor.Builder<CallingMetricAccumulator, CallingMetricAccumulator.Result> builder =
                VariantProcessor.Builder
                        .generatingAccumulatorsBy(() -> {
                            CallingMetricAccumulator accumulator = GVCF_INPUT ? new GvcfMetricAccumulator(dbsnp) : new CallingMetricAccumulator(dbsnp);
                            accumulator.setup(vcfHeader);
                            return accumulator;
                        })
                        .combiningResultsBy(CallingMetricAccumulator.Result::merge)
                        .withInput(INPUT)
                        .multithreadingBy(THREAD_COUNT);

        if (TARGET_INTERVALS != null) {
            builder.limitingProcessedRegionsTo(IntervalList.fromFile(TARGET_INTERVALS).uniqued());
        }

        final CallingMetricAccumulator.Result result = builder.build().process();

        // Fetch and write the metrics.
        final MetricsFile<CollectVariantCallingMetrics.VariantCallingDetailMetrics, Integer> detail = getMetricsFile();
        final MetricsFile<CollectVariantCallingMetrics.VariantCallingSummaryMetrics, Integer> summary = getMetricsFile();
        summary.addMetric(result.summary);
        result.details.forEach(detail::addMetric);

        final String outputPrefix = OUTPUT.getAbsolutePath() + ".";
        detail.write(new File(outputPrefix + CollectVariantCallingMetrics.VariantCallingDetailMetrics.getFileExtension()));
        summary.write(new File(outputPrefix + CollectVariantCallingMetrics.VariantCallingSummaryMetrics.getFileExtension()));

        return 0;
    }

    /** A collection of metrics relating to snps and indels within a variant-calling file (VCF). */
    public static class VariantCallingSummaryMetrics extends MetricBase {
        /** The number of high confidence SNPs calls (i.e. non-reference genotypes) that were examined */
        public long TOTAL_SNPS;

        /** The number of high confidence SNPs found in dbSNP */
        public long NUM_IN_DB_SNP;

        /** The number of high confidence SNPS called that were not found in dbSNP */
        public long NOVEL_SNPS;

        /** The number of SNPs that are also filtered */
        public long FILTERED_SNPS;

        /** The percentage of high confidence SNPs in dbSNP */
        public float PCT_DBSNP;

        /** The Transition/Transversion ratio of the SNP calls made at dbSNP sites */
        public double DBSNP_TITV;

        /** The Transition/Transversion ratio of the SNP calls made at non-dbSNP sites */
        public double NOVEL_TITV;

        /** The number of high confidence Indel calls that were examined */
        public long TOTAL_INDELS;

        /** The number of high confidence Indels called that were not found in dbSNP */
        public long NOVEL_INDELS;

        /** The number of indels that are also filtered */
        public long FILTERED_INDELS;

        /** The percentage of high confidence Indels in dbSNP */
        public float PCT_DBSNP_INDELS;

        /** The number of high confidence Indels found in dbSNP */
        public long NUM_IN_DB_SNP_INDELS;

        /** The Insertion/Deletion ratio of the Indel calls made at dbSNP sites */
        public double DBSNP_INS_DEL_RATIO;

        /** The Insertion/Deletion ratio of the Indel calls made at non-dbSNP sites */
        public double NOVEL_INS_DEL_RATIO;

        /** The number of high confidence multiallelic SNP calls that were examined */
        public double TOTAL_MULTIALLELIC_SNPS;

        /** The number of high confidence multiallelic SNPs found in dbSNP */
        public double NUM_IN_DB_SNP_MULTIALLELIC;

        /** The number of high confidence complex Indel calls that were examined */
        public double TOTAL_COMPLEX_INDELS;

        /** The number of high confidence complex Indels found in dbSNP */
        public double NUM_IN_DB_SNP_COMPLEX_INDELS;

        /** The rate at which reference bases are observed at ref/alt heterozygous SNP sites. */
        public double SNP_REFERENCE_BIAS;

        /**
         * For summary metrics, the number of variants that appear in only one sample.
         * For detail metrics, the number of variants that appear only in the current sample.
         */
        public long NUM_SINGLETONS;

        /** Hidden fields that are not propagated to the metrics file, but are used to house temporary values. */
        long refAlleleObs, altAlleleObs, novelDeletions, novelInsertions, novelTransitions, novelTransversions, dbSnpDeletions,
                dbSnpInsertions, dbSnpTransitions, dbSnpTransversions;

        public static String getFileExtension() {
            return "variant_calling_summary_metrics";
        }

        public void updateDerivedValuesInPlace() {
            this.PCT_DBSNP = this.NUM_IN_DB_SNP / (float) this.TOTAL_SNPS;
            this.NOVEL_SNPS = this.TOTAL_SNPS - this.NUM_IN_DB_SNP;
            this.SNP_REFERENCE_BIAS = (this.refAlleleObs) / (double) (this.refAlleleObs + this.altAlleleObs);

            if (this.dbSnpTransversions > 0)
                this.DBSNP_TITV = (double) this.dbSnpTransitions / (double) this.dbSnpTransversions;

            if (this.novelTransversions > 0)
                this.NOVEL_TITV = this.novelTransitions / (double) this.novelTransversions;

            this.PCT_DBSNP_INDELS = this.NUM_IN_DB_SNP_INDELS / (float) this.TOTAL_INDELS;
            this.NOVEL_INDELS = this.TOTAL_INDELS - this.NUM_IN_DB_SNP_INDELS;

            if (this.dbSnpDeletions > 0)
                this.DBSNP_INS_DEL_RATIO = this.dbSnpInsertions / (double) this.dbSnpDeletions;

            if (this.novelDeletions > 0)
                this.NOVEL_INS_DEL_RATIO = this.novelInsertions / (double) this.novelDeletions;
        }

        public static <T extends VariantCallingSummaryMetrics> void foldInto(final T target, final Collection<T> metrics) {
            for (final T metric : metrics) {
                target.TOTAL_SNPS += metric.TOTAL_SNPS;
                target.NUM_IN_DB_SNP += metric.NUM_IN_DB_SNP;
                target.FILTERED_SNPS += metric.FILTERED_SNPS;
                target.TOTAL_INDELS += metric.TOTAL_INDELS;
                target.FILTERED_INDELS += metric.FILTERED_INDELS;
                target.NUM_IN_DB_SNP_INDELS += metric.NUM_IN_DB_SNP_INDELS;
                target.TOTAL_MULTIALLELIC_SNPS += metric.TOTAL_MULTIALLELIC_SNPS;
                target.NUM_IN_DB_SNP_MULTIALLELIC += metric.NUM_IN_DB_SNP_MULTIALLELIC;
                target.TOTAL_COMPLEX_INDELS += metric.TOTAL_COMPLEX_INDELS;
                target.NUM_IN_DB_SNP_COMPLEX_INDELS += metric.NUM_IN_DB_SNP_COMPLEX_INDELS;
                target.NUM_SINGLETONS += metric.NUM_SINGLETONS;
                target.refAlleleObs += metric.refAlleleObs;
                target.altAlleleObs += metric.altAlleleObs;
                target.novelDeletions += metric.novelDeletions;
                target.novelInsertions += metric.novelInsertions;
                target.novelTransitions += metric.novelTransitions;
                target.novelTransversions += metric.novelTransversions;
                target.dbSnpDeletions += metric.dbSnpDeletions;
                target.dbSnpInsertions += metric.dbSnpInsertions;
                target.dbSnpTransitions += metric.dbSnpTransitions;
                target.dbSnpTransversions += metric.dbSnpTransversions;
            }
            target.updateDerivedValuesInPlace();
        }
    }

    /** A collection of metrics relating to snps and indels within a variant-calling file (VCF) for a given sample. */
    public static class VariantCallingDetailMetrics extends CollectVariantCallingMetrics.VariantCallingSummaryMetrics {
        /** The name of the sample being assayed */
        public String SAMPLE_ALIAS;

        /**
         * (count of hets)/(count of homozygous non-ref) for this sample
         */
        public double HET_HOMVAR_RATIO;

        /**
         * Hidden fields not propagated to the metrics file.
         */
        long numHets, numHomVar;

        public static String getFileExtension() {
            return "variant_calling_detail_metrics";
        }

        public static void foldInto(final VariantCallingDetailMetrics target,
                                    final Collection<VariantCallingDetailMetrics> metrics) {
            VariantCallingSummaryMetrics.foldInto(target, metrics);
            final Set<String> sampleAliases = new HashSet<String>();
            for (final VariantCallingDetailMetrics metric : metrics) {
                target.numHets += metric.numHets;
                target.numHomVar += metric.numHomVar;
                sampleAliases.add(metric.SAMPLE_ALIAS);
            }
            target.updateDerivedValuesInPlace();
            if (sampleAliases.size() != 1) {
                throw new IllegalArgumentException("Provided metrics do not have the same sample name.");
            }
            target.SAMPLE_ALIAS = sampleAliases.iterator().next();
        }

        @Override
        public void updateDerivedValuesInPlace() {
            super.updateDerivedValuesInPlace();
            // Divide by zero should be OK -- NaN should get propagated to metrics file and to DB.
            HET_HOMVAR_RATIO = numHets / (double) numHomVar;
        }
    }
}
