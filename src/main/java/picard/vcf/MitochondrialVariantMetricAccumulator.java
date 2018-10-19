package picard.vcf;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import picard.util.DbSnpBitSetUtil;

/**
 * Created by gauthier on 9/20/18.
 */
public class MitochondrialVariantMetricAccumulator extends CallingMetricAccumulator {

    public static final String filenameTag = "mito";
    public static final double HIGH_AF_THRESHOLD = 0.9;
    public static final double LOW_AF_THRESHOLD = 0.1;

    public MitochondrialVariantMetricAccumulator(final DbSnpBitSetUtil.DbSnpBitSets dbsnp) {
        super(dbsnp);
        sampleMetricsMap =
                new CollectionUtil.DefaultingMap<>(
                        sampleName -> {
                            final CollectVariantCallingMetrics.VariantCallingDetailMetrics detail = new CollectVariantCallingMetrics.MitochondrialVariantCallingDetailMetrics();
                            detail.SAMPLE_ALIAS = sampleName;
                            return detail;
                        }, true);
    }

    @Override
    protected void updateDetailMetric(final CollectVariantCallingMetrics.VariantCallingDetailMetrics metric,
                                      final Genotype genotype,
                                      final VariantContext vc,
                                      final boolean hasSingletonSample) {
        super.updateDetailMetric(metric, genotype, vc, hasSingletonSample);
        updateSummaryMetric(metric, genotype, vc, hasSingletonSample);

        final CollectVariantCallingMetrics.MitochondrialVariantCallingDetailMetrics mtMetric = (CollectVariantCallingMetrics.MitochondrialVariantCallingDetailMetrics)metric;

        if (genotype != null && !vc.isFiltered()) {
            if (genotype.hasExtendedAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
                if (vc.isBiallelic()) {
                    final Double alleleFraction = Double.parseDouble(genotype.getExtendedAttribute(VCFConstants.ALLELE_FREQUENCY_KEY).toString());
                    if (alleleFraction <= LOW_AF_THRESHOLD) {
                        ++mtMetric.TOTAL_LOW_AF_VARIANTS;
                    } else if (alleleFraction >= HIGH_AF_THRESHOLD) {
                        ++mtMetric.TOTAL_HOMOPLASMIC_VARIANTS;
                    } else {
                        ++mtMetric.TOTAL_MID_AF_VARIANTS;
                    }
                }
            }
        }
    }


    @Override
    public String getFilenameTag() {
        return filenameTag;
    }
}
