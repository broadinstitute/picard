package picard.vcf;

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

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFHeader;
import picard.util.DbSnpBitSetUtil;
import picard.vcf.processor.VariantProcessor;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static picard.vcf.CollectVariantCallingMetrics.VariantCallingDetailMetrics;
import static picard.vcf.CollectVariantCallingMetrics.VariantCallingSummaryMetrics;

/**
 * Collects variants and generates metrics about them.  To use, construct, call {@link #setup(VCFHeader)} once, then
 * {@link #accumulate(htsjdk.variant.variantcontext.VariantContext)} as desired, then call {@link #result()}.
 *
 * @author mccowan
 */
public class CallingMetricAccumulator implements VariantProcessor.Accumulator<CallingMetricAccumulator.Result> {
    public static class Result {
        final VariantCallingSummaryMetrics summary;
        final Collection<VariantCallingDetailMetrics> details;

        Result(final VariantCallingSummaryMetrics summary, final Collection<VariantCallingDetailMetrics> details) {
            this.summary = summary;
            this.details = details;
        }

        public static Result merge(final Collection<Result> results) {
            final Collection<VariantCallingDetailMetrics> details = new ArrayList<>();
            final Collection<VariantCallingSummaryMetrics> summaries = new ArrayList<>();
            results.stream().forEach(result -> {
                summaries.add(result.summary);
                details.addAll(result.details);
            });

            final Map<String, List<VariantCallingDetailMetrics>> sampleDetailsMap =
                    details.stream().collect(Collectors.groupingBy(vcDetailMetrics -> vcDetailMetrics.SAMPLE_ALIAS));

            final Collection<CollectVariantCallingMetrics.VariantCallingDetailMetrics> collapsedDetails = new ArrayList<>();

            sampleDetailsMap.values().stream().forEach(sampleDetails -> {
                final VariantCallingDetailMetrics collapsed = new VariantCallingDetailMetrics();
                VariantCallingDetailMetrics.foldInto(collapsed, sampleDetails);
                collapsedDetails.add(collapsed);
            });

            final VariantCallingSummaryMetrics collapsedSummary = new VariantCallingSummaryMetrics();
            VariantCallingSummaryMetrics.foldInto(collapsedSummary, summaries);

            return new Result(collapsedSummary, collapsedDetails);
        }
    }

    private static final Log LOG = Log.getInstance(CallingMetricAccumulator.class);
    private static final ProgressLogger progress = new ProgressLogger(LOG, 10000);

    private final DbSnpBitSetUtil.DbSnpBitSets dbsnp;
    private final VariantCallingSummaryMetrics summaryMetric = new VariantCallingSummaryMetrics();
    /**
     * A map of sample names to metrics.  If .get() for a not-yet-existing sample name, a metric is generated, inserted into the map,
     * then returned.
     */
    private final CollectionUtil.DefaultingMap<String, VariantCallingDetailMetrics> sampleMetricsMap =
            new CollectionUtil.DefaultingMap<>(
                    sampleName -> {
                        final VariantCallingDetailMetrics detail = new VariantCallingDetailMetrics();
                        detail.SAMPLE_ALIAS = sampleName;
                        return detail;
                    }, true);

    public CallingMetricAccumulator(final DbSnpBitSetUtil.DbSnpBitSets dbsnp) {
        this.dbsnp = dbsnp;
    }

    public void setup(final VCFHeader vcfHeader) {
        //noop.
    }

    /** Incorporates the provided variant's data into the metric analysis. */
    @Override
    public void accumulate(final VariantContext vc) {
        progress.record(vc.getContig(), vc.getStart());
        if (!isVariantExcluded(vc)) {
            final String singletonSample = getSingletonSample(vc);
            updateSummaryMetric(summaryMetric, null, vc, singletonSample != null); // The summary metric has no genotype.

            vc.getSampleNames().stream()
                    .filter(sampleName -> !vc.getGenotype(sampleName).isHomRef())
                    .forEach(sampleName ->
                            updateDetailMetric(sampleMetricsMap.get(sampleName), vc.getGenotype(sampleName), vc,
                                    sampleName.equals(singletonSample)));
        }
    }

    /**
     * @return Sample name if there is only one sample that contains alternate allele(s), else null if either multiple samples that
     * are not homref, or no samples that are not homref.
     */
    protected static String getSingletonSample(final VariantContext vc) {

        // peek can only change effectively final variables...workaround
        final String[] sampleName = new String[1];

        if (vc.getGenotypes()
                .stream()
                        // look at het or homVar genotypes
                .filter(genotype -> genotype.isHet() || genotype.isHomVar())
                        // two such genotypes will be enough
                .limit(2)
                        //get any of the sample names
                .peek(genotype -> sampleName[0] = genotype.getSampleName())
                        //map to the number of variant chromosomes
                .mapToInt(genotype -> genotype.isHet() ? 1 : 2)
                        //add them up
                .reduce(Integer::sum)
                        // compare to 1 with 0 as default
                .orElse(0) == 1) {
            return sampleName[0];
        } else {
            return null;
        }
    }

    public Result result() {
        final Collection<VariantCallingDetailMetrics> values = sampleMetricsMap.values();
        values.forEach(CollectVariantCallingMetrics.VariantCallingDetailMetrics::updateDerivedValuesInPlace);

        summaryMetric.updateDerivedValuesInPlace();
        return new Result(summaryMetric, values);
    }

    /** Returns true if the variant is --NOT-- interesting enough to be included in metrics calculations. */
    static private boolean isVariantExcluded(final VariantContext vc) {

        // If the entire record is not a variant, exclude it
        return !vc.isVariant() || vc.getGenotypes().stream().allMatch(Genotype::isHomRef);
    }

    private void updateDetailMetric(final VariantCallingDetailMetrics metric,
                                    final Genotype genotype,
                                    final VariantContext vc,
                                    final boolean hasSingletonSample) {
        updateSummaryMetric(metric, genotype, vc, hasSingletonSample);

        if (genotype != null && !vc.isFiltered()) {
            if (genotype.isHet()) {
                ++metric.numHets;
            } else if (genotype.isHomVar()) {
                ++metric.numHomVar;
            }
        }
    }

    /** Amends the provided metric with the data in the provided variant.  Also amends the summary metric re: reference bias. */
    private void updateSummaryMetric(final VariantCallingSummaryMetrics metric,
                                     final Genotype genotype,
                                     final VariantContext vc,
                                     final boolean hasSingletonSample) {

        // If this sample's genotype doesn't have any variation, exclude it
        if (genotype != null && genotype.isNoCall()) return;

        // Tally up the filtered SNPs & indels, then exit. The other metrics shouldn't be
        // computed on low-confidence calls.
        if (vc.isFiltered()) {
            if (vc.isSNP()) metric.FILTERED_SNPS++;
            else if (vc.isIndel()) metric.FILTERED_INDELS++;
            return;
        }

        if (hasSingletonSample) {
            ++metric.NUM_SINGLETONS;
        }

        if (vc.isBiallelic() && vc.isSNP()) {
            // Biallelic SNPs
            final boolean isInDbSnp = dbsnp.snps.isDbSnpSite(vc.getContig(), vc.getStart());
            final boolean isTransition = VariantContextUtils.isTransition(vc);

            metric.TOTAL_SNPS++;

            if (isInDbSnp) {
                metric.NUM_IN_DB_SNP++;
                if (isTransition) metric.dbSnpTransitions++;
                else metric.dbSnpTransversions++;
            } else {
                if (isTransition) metric.novelTransitions++;
                else metric.novelTransversions++;
            }

            // Calculate reference bias numbers.  Note: genotype == null for summary metric, so this block won't be called when metric ==
            // summaryMetric.
            if (genotype != null && genotype.isHet()) { //
                final int[] alleleDepths = genotype.getAD();
                /*
                 * Null check: work around GATK issue in which some biallelic sites are missing allele depth.  This should affect only ~1%
                 * of samples and should not have a significant impact on the reference bias calculation.
                 */
                if (alleleDepths != null) {
                    final int indexOfRef = vc.getAlleleIndex(vc.getReference());
                    final int indexOfAlt = (indexOfRef + 1) % 2;

                    metric.refAlleleObs += alleleDepths[indexOfRef];
                    metric.altAlleleObs += alleleDepths[indexOfAlt];

                    // Always count these values for summary metrics.
                    summaryMetric.refAlleleObs += alleleDepths[indexOfRef];
                    summaryMetric.altAlleleObs += alleleDepths[indexOfAlt];
                } else {
                    LOG.debug("Skipping aggregation of genotype due to missing allele depth data: ", genotype, ".");
                }
            }
        } else if (vc.isSNP() && vc.getAlternateAlleles().size() > 1) {
            // Multiallelic SNPs
            metric.TOTAL_MULTIALLELIC_SNPS++;
            if (dbsnp.snps.isDbSnpSite(vc.getContig(), vc.getStart())) metric.NUM_IN_DB_SNP_MULTIALLELIC++;
        } else if (vc.isIndel() && !vc.isComplexIndel()) {
            // Simple Indels
            final boolean isInDbSnp = dbsnp.indels.isDbSnpSite(vc.getContig(), vc.getStart());
            final boolean isInsertion = vc.isSimpleInsertion();

            metric.TOTAL_INDELS++;

            if (isInDbSnp) {
                metric.NUM_IN_DB_SNP_INDELS++;
                if (isInsertion) metric.dbSnpInsertions++;
                else metric.dbSnpDeletions++;
            } else {
                if (isInsertion) metric.novelInsertions++;
                else {
                    metric.novelDeletions++;
                }
            }
        } else if (vc.isComplexIndel()) {
            // Complex Indels
            metric.TOTAL_COMPLEX_INDELS++;
            if (dbsnp.indels.isDbSnpSite(vc.getContig(), vc.getStart())) metric.NUM_IN_DB_SNP_COMPLEX_INDELS++;
        }
    }
}
