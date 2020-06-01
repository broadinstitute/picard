package picard.arrays;

import picard.arrays.illumina.InfiniumVcfFields;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import picard.pedigree.Sex;
import picard.util.DbSnpBitSetUtil;
import picard.vcf.CallingMetricAccumulator;
import picard.vcf.processor.VariantProcessor;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;
import java.util.stream.Collectors;

class ArraysCallingMetricAccumulator implements VariantProcessor.Accumulator<ArraysCallingMetricAccumulator.Result> {

    private static final Log LOG = Log.getInstance(ArraysCallingMetricAccumulator.class);
    private static final ProgressLogger progress = new ProgressLogger(LOG, 10000);

    private final DbSnpBitSetUtil.DbSnpBitSets dbsnp;
    private final CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics summaryMetric
            = new CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics();

    private String sampleAlias;
    private Integer analysisVersionNumber;

    private String chipTypeName;
    private String reportedGender;
    private String fingerprintGender;
    private String autocallGender;
    private Iso8601Date autocallDate;
    private Iso8601Date imagingDate;
    private String autocallVersion;
    private String extendedIlluminaManifestVersion;
    private String zcallVersion;
    private String zcallThresholdsFile;
    private String clusterFileName;
    private Integer p95Green;
    private Integer p95Red;
    private String scannerName;
    private String pipelineVersion;

    /**
     * A map of sample names to metrics.  If .get() for a not-yet-existing sample name, a metric is generated, inserted into the map,
     * then returned.
     */
    private final CollectionUtil.DefaultingMap<String, CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics> sampleMetricsMap =
            new CollectionUtil.DefaultingMap<>(
                    sampleName -> {
                        final CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics detail = new CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics();
                        detail.CHIP_WELL_BARCODE = sampleName;
                        detail.SAMPLE_ALIAS = sampleAlias;
                        detail.ANALYSIS_VERSION = analysisVersionNumber;

                        detail.CHIP_TYPE = chipTypeName;
                        detail.REPORTED_GENDER = reportedGender;
                        detail.FP_GENDER = fingerprintGender;
                        detail.AUTOCALL_GENDER = autocallGender;

                        detail.AUTOCALL_VERSION = autocallVersion;
                        detail.AUTOCALL_DATE = new Iso8601Date(autocallDate);
                        detail.IMAGING_DATE = new Iso8601Date(imagingDate);
                        detail.EXTENDED_MANIFEST_VERSION = extendedIlluminaManifestVersion;
                        detail.ZCALL_VERSION = zcallVersion;
                        detail.zcallThresholdsFile = zcallThresholdsFile;
                        detail.CLUSTER_FILE_NAME = clusterFileName;
                        detail.P95_GREEN = p95Green;
                        detail.P95_RED = p95Red;
                        detail.SCANNER_NAME = scannerName;
                        detail.PIPELINE_VERSION = pipelineVersion;
                        return detail;
                    }, true);

    ArraysCallingMetricAccumulator(DbSnpBitSetUtil.DbSnpBitSets dbsnp) {
        this.dbsnp = dbsnp;
    }

    public void setup(final VCFHeader vcfHeader) {
        this.sampleAlias = InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.SAMPLE_ALIAS);
        this.analysisVersionNumber = InfiniumVcfFields.getOptionalIntegerFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.ANALYSIS_VERSION_NUMBER);
        this.chipTypeName = InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.ARRAY_TYPE);
        this.reportedGender = getOptionalGenderStringFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.EXPECTED_GENDER);
        this.fingerprintGender = getOptionalGenderStringFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.FINGERPRINT_GENDER);
        this.autocallGender = Sex.fromString(InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.AUTOCALL_GENDER)).toSymbol();
        this.autocallVersion = InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.AUTOCALL_VERSION);
        final SimpleDateFormat autocallDateFormat = new SimpleDateFormat("MM/dd/yyyy HH:mm");         // of the form '09/21/2016 20:40'
        this.autocallDate = InfiniumVcfFields.getDateFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.AUTOCALL_DATE, autocallDateFormat);
        final SimpleDateFormat imagingDateFormat = new SimpleDateFormat("MM/dd/yyyy hh:mm:ss a");     // of the form '8/15/2015 7:28:52 AM'
        imagingDateFormat.setTimeZone(TimeZone.getTimeZone("America/New_York"));
        this.imagingDate = InfiniumVcfFields.getDateFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.IMAGING_DATE, imagingDateFormat);
        this.extendedIlluminaManifestVersion = InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.EXTENDED_ILLUMINA_MANIFEST_VERSION);
        this.zcallVersion = InfiniumVcfFields.getOptionalValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.ZCALL_VERSION);
        this.zcallThresholdsFile = InfiniumVcfFields.getOptionalValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.ZCALL_THRESHOLDS);
        this.clusterFileName = InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.CLUSTER_FILE);
        this.p95Green = InfiniumVcfFields.getIntegerFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.P_95_GREEN);
        this.p95Red = InfiniumVcfFields.getIntegerFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.P_95_RED);
        this.scannerName = InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.SCANNER_NAME);
        this.pipelineVersion = InfiniumVcfFields.getOptionalValueFromVcfOtherHeaderLine(vcfHeader, InfiniumVcfFields.PIPELINE_VERSION);

        vcfHeader.getGenotypeSamples().forEach(sampleName -> sampleMetricsMap.get(sampleName));
    }

    private String getOptionalGenderStringFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String fieldName) {
        String genderString = InfiniumVcfFields.getOptionalValueFromVcfOtherHeaderLine(vcfHeader, fieldName);
        if (genderString != null) {
            return Sex.fromString(genderString).toSymbol();
        }
        return Sex.NotReported.toSymbol();
    }

    @Override
    public void accumulate(VariantContext vc) {
        progress.record(vc.getContig(), vc.getStart());
        final String singletonSample = CallingMetricAccumulator.getSingletonSample(vc);

        vc.getSampleNames().forEach(sampleName ->
            updateDetailMetric(sampleMetricsMap.get(sampleName), vc.getGenotype(sampleName), vc,
                    sampleName.equals(singletonSample)));
    }

    private void updateDetailMetric(final CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics metric,
                                    final Genotype genotype,
                                    final VariantContext vc,
                                    final boolean hasSingletonSample) {
        metric.NUM_ASSAYS++;
        if (!vc.isFiltered() || vc.getCommonInfo().getFilters().contains(InfiniumVcfFields.DUPE)) {
            metric.NUM_NON_FILTERED_ASSAYS++;
            if (genotype.isCalled()) {
                metric.NUM_CALLS++;
                String gtA = (String) genotype.getExtendedAttribute(InfiniumVcfFields.GTA, genotype.getGenotypeString());
                if (!gtA.equals(VCFConstants.EMPTY_GENOTYPE)) {
                    metric.NUM_AUTOCALL_CALLS++;
                }
            } else {
                metric.NUM_NO_CALLS++;
            }

            if (vc.isSNP()) {
                // Biallelic SNPs
                final boolean isInDbSnp = dbsnp.snps.isDbSnpSite(vc.getContig(), vc.getStart());

                metric.NUM_SNPS++;

                if (isInDbSnp) {
                    metric.NUM_IN_DB_SNP++;
                }
            } else if (vc.isIndel()) {
                metric.NUM_INDELS++;
            }

            if (hasSingletonSample) {
                metric.NUM_SINGLETONS++;
            }

            if (genotype.isHet()) {
                metric.numHets++;
            } else if (genotype.isHomVar()) {
                metric.numHomVar++;
            }
        } else {
            metric.NUM_FILTERED_ASSAYS++;
            if (vc.getCommonInfo().getFilters().contains(InfiniumVcfFields.ZEROED_OUT_ASSAY)) {
                // A "zeroed-out SNP".  Marked as unusable/uncallable
                metric.NUM_ZEROED_OUT_ASSAYS++;
            }
        }
    }

    @Override
    public Result result() {
        return new Result(summaryMetric, sampleMetricsMap.values());
    }

    public static class Result {
        final CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics summary;
        final Collection<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics> details;

        Result(final CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics summary, final Collection<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics> details) {
            this.summary = summary;
            this.details = details;
        }

        /**
         * Combines results generated across 1 or more threads
         */
        public static Result merge(final Collection<Result> results) {
            final Collection<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics> details = new ArrayList<>();
            results.forEach(result -> {
                details.addAll(result.details);
            });

            final Map<String, List<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics>> sampleDetailsMap =
                    details.stream().collect(Collectors.groupingBy(vcDetailMetrics -> vcDetailMetrics.CHIP_WELL_BARCODE));

            final Collection<CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics> collapsedDetails = new ArrayList<>();


            final CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics collapsedSummary = new CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics();
            sampleDetailsMap.values().forEach(sampleDetails -> {
                final CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics collapsed = new CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics();
                CollectArraysVariantCallingMetrics.ArraysVariantCallingDetailMetrics.foldInto(collapsed, sampleDetails);
                CollectArraysVariantCallingMetrics.ArraysVariantCallingSummaryMetrics.foldInto(collapsedSummary, sampleDetails);
                collapsedDetails.add(collapsed);
                collapsed.calculateDerivedFields();
            });
            collapsedSummary.calculateDerivedFields();

            return new Result(collapsedSummary, collapsedDetails);
        }
    }
}
