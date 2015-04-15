package picard.vcf;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import picard.util.DbSnpBitSetUtil;
import picard.vcf.processor.VariantProcessor;

import static picard.vcf.CollectVariantCallingMetrics.VariantCallingSummaryMetrics;
import static picard.vcf.CollectVariantCallingMetrics.VariantCallingDetailMetrics;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

/**
 * Collects variants and generates metrics about them.  To use, construct, call
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
            final Collection<VariantCallingDetailMetrics> details = new ArrayList<VariantCallingDetailMetrics>();
            final Collection<VariantCallingSummaryMetrics> summaries = new ArrayList<VariantCallingSummaryMetrics>();
            for (final Result result : results) {
                summaries.add(result.summary);
                details.addAll(result.details);
            }


            final Map<String, Collection<CollectVariantCallingMetrics.VariantCallingDetailMetrics>> sampleDetailsMap = CollectionUtil.partition(details,
                    new CollectionUtil.Partitioner<CollectVariantCallingMetrics.VariantCallingDetailMetrics, String>() {
                        @Override
                        public String getPartition(final CollectVariantCallingMetrics.VariantCallingDetailMetrics variantCallingDetailMetrics) {
                            return variantCallingDetailMetrics.SAMPLE_ALIAS;
                        }
                    });
            final Collection<CollectVariantCallingMetrics.VariantCallingDetailMetrics> collapsedDetails = new ArrayList<VariantCallingDetailMetrics>();
            for (final Collection<VariantCallingDetailMetrics> sampleDetails : sampleDetailsMap.values()) {
                final VariantCallingDetailMetrics collapsed = new VariantCallingDetailMetrics();
                VariantCallingDetailMetrics.foldInto(collapsed, sampleDetails);
                collapsedDetails.add(collapsed);
            }


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
            new CollectionUtil.DefaultingMap<String, VariantCallingDetailMetrics>(
                    new CollectionUtil.DefaultingMap.Factory<VariantCallingDetailMetrics, String>() {
                        @Override
                        public VariantCallingDetailMetrics make(final String sampleName) {
                            final VariantCallingDetailMetrics detail = new VariantCallingDetailMetrics();
                            detail.SAMPLE_ALIAS = sampleName;
                            return detail;
                        }
                    }, true);

    public CallingMetricAccumulator(final DbSnpBitSetUtil.DbSnpBitSets dbsnp) {
        this.dbsnp = dbsnp;
    }

    /** Incorporates the provided variant's data into the metric analysis. */
    @Override
    public void accumulate(final VariantContext vc) {
        progress.record(vc.getChr(), vc.getStart());
        if (!isVariantExcluded(vc)) {
            final String singletonSample = getSingletonSample(vc);
            updateSummaryMetric(summaryMetric, null, vc, singletonSample != null); // The summary metric has no genotype.
            for (final String sampleName : vc.getSampleNames()) {
                // Skip homozygous reference calls.
                if (!vc.getGenotype(sampleName).isHomRef()) {
                    updateDetailMetric(sampleMetricsMap.get(sampleName), vc.getGenotype(sampleName), vc,
                            sampleName.equals(singletonSample));
                }
            }
        }
    }

    /**
     * @return Sample name if there is only one sample that contains alternate allele(s), else null if either multiple samples that
     * are not homref, or no samples that are not homref.
     */
    private String getSingletonSample(final VariantContext vc) {
        String singletonSample = null;
        for (final String sampleName : vc.getSampleNames()) {
            final Genotype genotype = vc.getGenotype(sampleName);
            if (genotype.isHomVar()) {
                return null;
            } else if (genotype.isHet()) {
                if (singletonSample != null) {
                    // second sample with non-reference allele, so not a singleton
                    return null;
                } else {
                    singletonSample = sampleName;
                }
            }
        }
        return singletonSample;
    }

    public Result result() {
        final Collection<VariantCallingDetailMetrics> values = sampleMetricsMap.values();
        for (final VariantCallingDetailMetrics value : values) {
            value.updateDerivedValuesInPlace();
        }

        summaryMetric.updateDerivedValuesInPlace();
        return new Result(summaryMetric, values);
    }

    /** Returns true if the variant is --NOT-- interesting enough to be included in metrics calculations. */
    private boolean isVariantExcluded(final VariantContext vc) {

        // If the entire record is not a variant, exclude it
        if (!vc.isVariant()) {
            return true;
        }

        // Skip calls which are homozygous reference for all samples.
        for (final String sample : vc.getSampleNames()) {
            if (!vc.getGenotype(sample).isHomRef()) {
                return false; // TODO: Is this correct?
            }
        }

        return true;
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
            final boolean isInDbSnp = dbsnp.snps.isDbSnpSite(vc.getChr(), vc.getStart());
            final boolean isTransition = isTransition(vc);

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
            if (dbsnp.snps.isDbSnpSite(vc.getChr(), vc.getStart())) metric.NUM_IN_DB_SNP_MULTIALLELIC++;
        } else if (vc.isIndel() && !vc.isComplexIndel()) {
            // Simple Indels
            final boolean isInDbSnp = dbsnp.indels.isDbSnpSite(vc.getChr(), vc.getStart());
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
            if (dbsnp.indels.isDbSnpSite(vc.getChr(), vc.getStart())) metric.NUM_IN_DB_SNP_COMPLEX_INDELS++;
        }

    }

    /**
     * Answers if the provided variant is transitional (otherwise, it's transversional).
     * Transitions:
     * A->G
     * G->A
     * C->T
     * T->C
     * <p/>
     * Transversions:
     * A->C
     * A->T
     * C->A
     * C->G
     * G->C
     * G->T
     * T->A
     * T->G
     */
    static private boolean isTransition(final VariantContext vc) {
        final byte refAllele = vc.getReference().getBases()[0];
        final Collection<Allele> altAlleles = vc.getAlternateAlleles();

        Byte altAllele = null;
        for (final Allele a : altAlleles) {
            if (a.getBases()[0] != refAllele) {
                altAllele = a.getBases()[0];
                break;
            }
        }
        if (altAllele == null) {
            // This should never happen
            throw new IllegalArgumentException("All alternate alleles match the reference base " + (char) refAllele);
        }

        return (refAllele == 'A' && altAllele == 'G')
                || (refAllele == 'G' && altAllele == 'A')
                || (refAllele == 'C' && altAllele == 'T')
                || (refAllele == 'T' && altAllele == 'C');
    }

}
