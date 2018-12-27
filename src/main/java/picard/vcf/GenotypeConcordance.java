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
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import htsjdk.tribble.Tribble;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import picard.vcf.GenotypeConcordanceStates.CallState;
import picard.vcf.GenotypeConcordanceStates.ContingencyState;
import picard.vcf.GenotypeConcordanceStates.TruthAndCallStates;
import picard.vcf.GenotypeConcordanceStates.TruthState;
import picard.vcf.PairedVariantSubContextIterator.VcfTuple;

import java.io.File;
import java.util.*;

import static htsjdk.variant.variantcontext.VariantContext.Type.*;
import static htsjdk.variant.vcf.VCFConstants.MISSING_VALUE_v4;
import static picard.vcf.GenotypeConcordanceStateCodes.*;

/**
 * <h3>Summary</h3>
 * Calculates the concordance between genotype data of one sample in each of two VCFs - one being considered the truth (or reference)
 * the other being the call.  The concordance is broken into separate results sections for SNPs and indels.  Satistics
 * are reported in three different files.
 *
 * <h3>Details</h3>
 * This tool evaluates the concordance between genotype calls for a sample in different callsets where one is being considered
 * as the \"truth\" (aka standard, or reference) and the other as the \"call\" that is being evaluated for accuracy. The
 * Comparison can be restricted to a confidence interval which is typically used in order to enable proper assessment of
 * False Positives and the False-Positive Rate (FPR).
 * <br />
 * <h3>Usage example</h3>
 * <h4>Compare two VCFs within a confidence region</h4>
 * <pre>
 * java -jar picard.jar GenotypeConcordance \\
 *       CALL_VCF=input.vcf \\
 *       CALL_SAMPLE=sample_name \\
 *       O=gc_concordance.vcf \\
 *       TRUTH_VCF=truth_set.vcf \\
 *       TRUTH_SAMPLE=sample_in_truth \\
 *       INTERVALS=confident.interval_list \\
 *       MISSING_SITES_HOM_REF = true
 * </pre>
 *
 * <h3>Output Metrics:</h3>
 * Output metrics consists of GenotypeConcordanceContingencyMetrics, GenotypeConcordanceSummaryMetrics, and
 * GenotypeConcordanceDetailMetrics.  For each set of metrics, the data is broken into separate sections for
 * SNPs and INDELs.  Note that only SNP and INDEL variants are considered, MNP, Symbolic, and Mixed classes
 *  of variants are not included.
 * <ul>
 * <li>{@link GenotypeConcordanceContingencyMetrics} enumerate the constituents of each contingent in a callset including true-positive
 * (TP), true-negative (TN), false-positive (FP), and false-negative (FN) calls.</li>
 * <li>{@link GenotypeConcordanceDetailMetrics} include the numbers of SNPs and INDELs for each contingent genotype as well as the
 * number of validated genotypes.</li>
 * <li>{@link GenotypeConcordanceSummaryMetrics} provide specific details for the variant caller performance on a callset including
 * values for sensitivity, specificity, and positive predictive values.</li>
 * </ul>
 * <br />
 * <br />
 * Useful definitions applicable to alleles and genotypes:
 * <ul>
 * <li>Truthset - A callset (typically in VCF format) containing variant calls and genotypes that have been cross-validated
 * with multiple technologies e.g. Genome In A Bottle Consortium (GIAB) (https://sites.stanford.edu/abms/giab)</li>
 * <li>TP - True-positives are variant sites that match against the truth-set</li>
 * <li>FP - False-positives are reference sites miscalled as variant</li>
 * <li>FN - False-negatives are variant sites miscalled as reference</li>
 * <li>TN - True-negatives are correctly called as reference</li>
 * <li>Validated genotypes - are TP sites where the exact genotype (HET or HOM-VAR) appears in the truth-set</li>
 * </ul>
 * <h3>VCF Output:</h3>
 * <ul>
 * <li>The concordance state will be stored in the CONC_ST tag in the INFO field</li>
 * <li>The truth sample name will be \"truth\" and call sample name will be \"call\"</li>
 * </ul>
 *
 *
 * @author Tim Fennell
 * @author George Grant
 */
@CommandLineProgramProperties(
        summary = GenotypeConcordance.USAGE_SUMMARY + GenotypeConcordance.USAGE_DETAILS,
        oneLineSummary =  GenotypeConcordance.USAGE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class)
@DocumentedFeature
public class GenotypeConcordance extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Calculates the concordance between genotype data of one sample in each of two VCFs -" +
            " truth (or reference) vs. calls.";
    static final String USAGE_DETAILS =
            "<h3>Summary</h3>" +
            "Calculates the concordance between genotype data of one sample in each of two VCFs - one being considered the truth (or reference) " +
            "the other being the call.  The concordance is broken into separate results sections for SNPs and indels.  Summary and detailed statistics " +
            "are reported.\n" +
            "\n" +
            "<h3>Details</h3>\n" +
            "This tool evaluates the concordance between genotype calls for a sample in different callsets where one is being considered " +
            "as the \"truth\" (aka standard, or reference) and the other as the \"call\" that is being evaluated for accuracy. The " +
            "Comparison can be restricted to a confidence interval which is typically used in order to enable proper assessment of " +
            "False Positives and the False-Positive Rate (FPR).\n" +
            "\n" +
            "<h3>Usage example</h3>\n" +
            "<h4>Compare two VCFs within a confidence region</h4>\n" +
            "\n" +
            "java -jar picard.jar GenotypeConcordance \\\n" +
            "      CALL_VCF=input.vcf \\\n" +
            "      CALL_SAMPLE=sample_name \\\n" +
            "      O=gc_concordance.vcf \\\n" +
            "      TRUTH_VCF=truth_set.vcf \\\n" +
            "      TRUTH_SAMPLE=sample_in_truth \\\n" +
            "      INTERVALS=confident.interval_list \\\n" +
            "      MISSING_SITES_HOM_REF = true\n" +
            "\n" +
            "<h3>Output Metrics:</h3>\n" +
            "Output metrics consists of GenotypeConcordanceContingencyMetrics, GenotypeConcordanceSummaryMetrics, and " +
            "GenotypeConcordanceDetailMetrics.  For each set of metrics, the data is broken into separate sections for " +
            "SNPs and INDELs.  Note that only SNP and INDEL variants are considered, MNP, Symbolic, and Mixed classes " +
            " of variants are not included.\n" +
            "\n" +
            "- GenotypeConcordanceContingencyMetrics enumerate the constituents of each contingent in a callset including true-positive" +
            "(TP), true-negative (TN), false-positive (FP), and false-negative (FN) calls. See" +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#GenotypeConcordanceContingencyMetrics for more details.\n" +
            "- GenotypeConcordanceDetailMetrics include the numbers of SNPs and INDELs for each contingent genotype as well as the" +
            "number of validated genotypes. See" +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#GenotypeConcordanceDetailMetrics for more details." +
            "- GenotypeConcordanceSummaryMetrics provide specific details for the variant caller performance on a callset including:" +
            "values for sensitivity, specificity, and positive predictive values. See" +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#GenotypeConcordanceSummaryMetrics for more details.\n" +
            "\n" +
            "Useful definitions applicable to alleles and genotypes:\n" +
            "\n" +
            "Truthset - A callset (typically in VCF format) containing variant calls and genotypes that have been cross-validated" +
            "with multiple technologies e.g. Genome In A Bottle Consortium (GIAB) (https://sites.stanford.edu/abms/giab)\n" +
            "TP - True-positives are variant sites that match against the truth-set\n" +
            "FP - False-positives are reference sites miscalled as variant\n" +
            "FN - False-negatives are variant sites miscalled as reference\n" +
            "TN - True-negatives are correctly called as reference\n" +
            "Validated genotypes - are TP sites where the exact genotype (HET or HOM-VAR) appears in the truth-set\n" +
            "\n" +
            "<h3>VCF Output:</h3>\n" +
            "- The concordance state will be stored in the CONC_ST tag in the INFO field\n" +
            "- The truth sample name will be \"truth\" and call sample name will be \"call\"";
    @Argument(shortName = "TV", doc="The VCF containing the truth sample")
    public File TRUTH_VCF;

    @Argument(shortName = "CV", doc="The VCF containing the call sample")
    public File CALL_VCF;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Basename for the three metrics files that are to be written." +
            " Resulting files will be <OUTPUT>" + SUMMARY_METRICS_FILE_EXTENSION + ", <OUTPUT>" + DETAILED_METRICS_FILE_EXTENSION + ", and <OUTPUT>" + CONTINGENCY_METRICS_FILE_EXTENSION + ".")
    public File OUTPUT;

    @Argument(doc = "Output a VCF annotated with concordance information.")
    public boolean OUTPUT_VCF = false;

    @Argument(shortName = "TS", doc="The name of the truth sample within the truth VCF. Not required if only one sample exists.", optional = true)
    public String TRUTH_SAMPLE = null;

    @Argument(shortName = "CS", doc="The name of the call sample within the call VCF. Not required if only one sample exists.", optional = true)
    public String CALL_SAMPLE = null;

    @Argument(doc="One or more interval list files that will be used to limit the genotype concordance.  Note - if intervals are specified, the VCF files must be indexed.", optional = true)
    public List<File> INTERVALS;

    @Argument(doc="If true, multiple interval lists will be intersected. If false multiple lists will be unioned.")
    public boolean INTERSECT_INTERVALS = true;

    @Argument(doc="Genotypes below this genotype quality will have genotypes classified as LowGq.")
    public int MIN_GQ = 0;

    @Argument(doc="Genotypes below this depth will have genotypes classified as LowDp.")
    public int MIN_DP = 0;

    @Argument(doc="If true, output all rows in detailed statistics even when count == 0.  When false only output rows with non-zero counts.")
    public boolean OUTPUT_ALL_ROWS = false;

    @Argument(doc="If true, use the VCF index, else iterate over the entire VCF.", optional = true)
    public boolean USE_VCF_INDEX = false;

    @Argument(shortName = "MISSING_HOM", doc="Default is false, which follows the GA4GH Scheme. If true, missing sites in the truth set will be " +
            "treated as HOM_REF sites and sites missing in both the truth and call sets will be true negatives. Useful when hom ref sites are left out of the truth set. " +
            "This flag can only be used with a high confidence interval list.")
    public boolean MISSING_SITES_HOM_REF = false;

    @Argument(doc="Default is false. If true, filter status of sites will be ignored so that we include filtered sites when calculating genotype concordance. ", optional = true)
    public boolean IGNORE_FILTER_STATUS = false;

    private final Log log = Log.getInstance(GenotypeConcordance.class);
    private final ProgressLogger progress = new ProgressLogger(log, 10000, "checked", "variants");

    public static final String SUMMARY_METRICS_FILE_EXTENSION     = ".genotype_concordance_summary_metrics";
    public static final String DETAILED_METRICS_FILE_EXTENSION    = ".genotype_concordance_detail_metrics";
    public static final String CONTINGENCY_METRICS_FILE_EXTENSION = ".genotype_concordance_contingency_metrics";
    public static final String OUTPUT_VCF_FILE_EXTENSION          = ".genotype_concordance" + IOUtil.COMPRESSED_VCF_FILE_EXTENSION;

    protected GenotypeConcordanceCounts snpCounter;
    public GenotypeConcordanceCounts getSnpCounter() { return snpCounter; }

    protected GenotypeConcordanceCounts indelCounter;
    public GenotypeConcordanceCounts getIndelCounter() { return indelCounter; }

    public static final String CONTINGENCY_STATE_TAG = "CONC_ST";
    public static final VCFHeaderLine CONTINGENCY_STATE_HEADER_LINE = new VCFInfoHeaderLine(CONTINGENCY_STATE_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "The genotype concordance contingency state(s)");
    public static final String OUTPUT_VCF_TRUTH_SAMPLE_NAME = "truth";
    public static final String OUTPUT_VCF_CALL_SAMPLE_NAME = "call";

    public static void main(final String[] args) {
        new GenotypeConcordance().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        // Note - If the user specifies to use INTERVALS, the code will fail if the vcfs are not indexed, so we set USE_VCF_INDEX to true and check that the vcfs are indexed.
        IOUtil.assertFileIsReadable(TRUTH_VCF);
        IOUtil.assertFileIsReadable(CALL_VCF);
        final boolean usingIntervals = this.INTERVALS != null && !this.INTERVALS.isEmpty();
        final List<String> errors = new ArrayList<String>();
        if (usingIntervals) {
            USE_VCF_INDEX = true;
        }
        if (USE_VCF_INDEX) {
            // Index file is required either because we are using intervals, or because user-set parameter
            if (!indexExists(TRUTH_VCF)) {
                errors.add("The index file was not found for the TRUTH VCF.  Note that if intervals are specified, the VCF files must be indexed.");
            }
            if (!indexExists(CALL_VCF)) {
                errors.add("The index file was not found for the CALL VCF.  Note that if intervals are specified, the VCF files must be indexed.");
            }
        }
        if (MISSING_SITES_HOM_REF) {
            //If you are using this flag you must include a high confidence interval list where missing sites are hom_ref.
            if (!usingIntervals) {
                errors.add("You cannot use the MISSING_HOM option without also supplying an interval list over which missing " +
                        "sites are considered confident homozygous reference calls.");
            }
        }

        if (errors.isEmpty()) {
            return null;
        } else {
            return errors.toArray(new String[errors.size()]);
        }
    }

    /**
     * Determines whether an index file exists for the given vcf file using standard extension suffixes
     *
     * @param vcf  the vcf file to investigate
     * @return true if an index file exists, false otherwise
     */
    private boolean indexExists(final File vcf) {
        return Tribble.indexFile(vcf).exists() || Tribble.tabixIndexFile(vcf).exists();
    }

    @Override protected int doWork() {
        final File summaryMetricsFile     = new File(OUTPUT + SUMMARY_METRICS_FILE_EXTENSION);
        final File detailedMetricsFile    = new File(OUTPUT + DETAILED_METRICS_FILE_EXTENSION);
        final File contingencyMetricsFile = new File(OUTPUT + CONTINGENCY_METRICS_FILE_EXTENSION);
        IOUtil.assertFileIsWritable(summaryMetricsFile);
        IOUtil.assertFileIsWritable(detailedMetricsFile);
        IOUtil.assertFileIsWritable(contingencyMetricsFile);

        final GenotypeConcordanceSchemeFactory schemeFactory = new GenotypeConcordanceSchemeFactory();
        final GenotypeConcordanceScheme scheme = schemeFactory.getScheme(MISSING_SITES_HOM_REF);

        final boolean usingIntervals = this.INTERVALS != null && !this.INTERVALS.isEmpty();
        IntervalList intervals = null;
        SAMSequenceDictionary intervalsSamSequenceDictionary = null;
        if (usingIntervals) {
            log.info("Starting to load intervals list(s).");
            long genomeBaseCount = 0;
            for (final File f : INTERVALS) {
                IOUtil.assertFileIsReadable(f);
                final IntervalList tmpIntervalList = IntervalList.fromFile(f);
                if (genomeBaseCount == 0) {         // Don't count the reference length more than once.
                    intervalsSamSequenceDictionary = tmpIntervalList.getHeader().getSequenceDictionary();
                    genomeBaseCount = intervalsSamSequenceDictionary.getReferenceLength();
                }

                if (intervals == null)        intervals = tmpIntervalList;
                else if (INTERSECT_INTERVALS) intervals = IntervalList.intersection(intervals, tmpIntervalList);
                else intervals =              IntervalList.union(intervals, tmpIntervalList);
            }
            if (intervals != null) {
                intervals = intervals.uniqued();
            }
            log.info("Finished loading up intervals list(s).");
        }
        final VCFFileReader truthReader = new VCFFileReader(TRUTH_VCF, USE_VCF_INDEX);
        final VCFFileReader callReader = new VCFFileReader(CALL_VCF, USE_VCF_INDEX);

        if (TRUTH_SAMPLE == null) {
            if (truthReader.getFileHeader().getNGenotypeSamples() > 1) {
                throw new PicardException("TRUTH_SAMPLE is required when the TRUTH_VCF has more than one sample");
            }
            TRUTH_SAMPLE = truthReader.getFileHeader().getGenotypeSamples().get(0);
        }

        if (CALL_SAMPLE == null) {
            if (callReader.getFileHeader().getNGenotypeSamples() > 1) {
                throw new PicardException("CALL_SAMPLE is required when the CALL_VCF has more than one sample");
            }
            CALL_SAMPLE = callReader.getFileHeader().getGenotypeSamples().get(0);
        }

        // Check that the samples actually exist in the files!
        if (!truthReader.getFileHeader().getGenotypeSamples().contains(TRUTH_SAMPLE)) {
            throw new PicardException("File " + TRUTH_VCF.getAbsolutePath() + " does not contain genotypes for sample " + TRUTH_SAMPLE);
        }
        if (!callReader.getFileHeader().getGenotypeSamples().contains(CALL_SAMPLE)) {
            throw new PicardException("File " + CALL_VCF.getAbsolutePath() + " does not contain genotypes for sample " + CALL_SAMPLE);
        }

        // Verify that both VCFs have the same Sequence Dictionary
        SequenceUtil.assertSequenceDictionariesEqual(truthReader.getFileHeader().getSequenceDictionary(), callReader.getFileHeader().getSequenceDictionary());

        final Optional<VariantContextWriter> writer = getVariantContextWriter(truthReader, callReader);

        if (usingIntervals) {
            // If using intervals, verify that the sequence dictionaries agree with those of the VCFs
            SequenceUtil.assertSequenceDictionariesEqual(intervalsSamSequenceDictionary, truthReader.getFileHeader().getSequenceDictionary());
        }

        // Build the pair of iterators over the regions of interest
        final Iterator<VariantContext> truthIterator, callIterator;
        if (usingIntervals) {
            truthIterator = new ByIntervalListVariantContextIterator(truthReader, intervals);
            callIterator = new ByIntervalListVariantContextIterator(callReader, intervals);
        }
        else {
            truthIterator = truthReader.iterator();
            callIterator = callReader.iterator();
        }

        // Now do the iteration and count things up
        final PairedVariantSubContextIterator pairedIterator = new PairedVariantSubContextIterator(truthIterator, TRUTH_SAMPLE, callIterator, CALL_SAMPLE, truthReader.getFileHeader().getSequenceDictionary());
        snpCounter   = new GenotypeConcordanceCounts();
        indelCounter = new GenotypeConcordanceCounts();

        // A map to keep track of the count of Truth/Call States which we could not successfully classify
        final Map<String, Integer> unClassifiedStatesMap = new HashMap<String, Integer>();

        log.info("Starting iteration over variants.");
        while (pairedIterator.hasNext()) {
            final VcfTuple tuple = pairedIterator.next();
            final VariantContext.Type truthVariantContextType = tuple.leftVariantContext.map(VariantContext::getType).orElse(NO_VARIATION);
            final VariantContext.Type callVariantContextType  = tuple.rightVariantContext.map(VariantContext::getType).orElse(NO_VARIATION);

            final boolean stateClassified = classifyVariants(tuple.leftVariantContext, TRUTH_SAMPLE,
                    tuple.rightVariantContext, CALL_SAMPLE,
                    Optional.of(snpCounter), Optional.of(indelCounter),
                    MIN_GQ, MIN_DP, IGNORE_FILTER_STATUS);

            if (!stateClassified) {
                final String condition = truthVariantContextType + " " + callVariantContextType;
                final Integer count = unClassifiedStatesMap.getOrDefault(condition, 0) + 1;
                unClassifiedStatesMap.put(condition, count);
            }

            // write to the output VCF
            writer.ifPresent(w -> writeVcfTuple(tuple, w, scheme));

            //final VariantContext variantContextForLogging = tuple.leftVariantContext.orElseGet(tuple.rightVariantContext::get); // FIXME
            final VariantContext variantContextForLogging = tuple.leftVariantContext.isPresent() ? tuple.leftVariantContext.get() : tuple.rightVariantContext.get();
            progress.record(variantContextForLogging.getContig(), variantContextForLogging.getStart());
        }

        //snp counter add in X number of missing-missing hom ref's (truth and call state)
        //missing missing is total interval size minus number of iterations in while loop
        if (MISSING_SITES_HOM_REF) {
            // need to know size of region called over (intervals or whole genome) to add missing-missing sites for NIST schema.
            final long baseCount = (intervals != null) ? intervals.getBaseCount() : truthReader.getFileHeader().getSequenceDictionary().getReferenceLength();
            addMissingTruthAndMissingCallStates(snpCounter.getCounterSize(), baseCount, snpCounter);
            addMissingTruthAndMissingCallStates(indelCounter.getCounterSize(), baseCount, indelCounter);
        }

        // Calculate and store the summary-level metrics
        final MetricsFile<GenotypeConcordanceSummaryMetrics,?> genotypeConcordanceSummaryMetricsFile = getMetricsFile();
        GenotypeConcordanceSummaryMetrics summaryMetrics = new GenotypeConcordanceSummaryMetrics(SNP, snpCounter, TRUTH_SAMPLE, CALL_SAMPLE, MISSING_SITES_HOM_REF);
        genotypeConcordanceSummaryMetricsFile.addMetric(summaryMetrics);
        summaryMetrics = new GenotypeConcordanceSummaryMetrics(INDEL, indelCounter, TRUTH_SAMPLE, CALL_SAMPLE, MISSING_SITES_HOM_REF);
        genotypeConcordanceSummaryMetricsFile.addMetric(summaryMetrics);
        genotypeConcordanceSummaryMetricsFile.write(summaryMetricsFile);

        // Calculate and store the detailed metrics for both SNP and indels
        final MetricsFile<GenotypeConcordanceDetailMetrics,?> genotypeConcordanceDetailMetrics = getMetricsFile();
        outputDetailMetricsFile(SNP, genotypeConcordanceDetailMetrics, snpCounter, TRUTH_SAMPLE, CALL_SAMPLE, MISSING_SITES_HOM_REF, OUTPUT_ALL_ROWS);
        outputDetailMetricsFile(INDEL, genotypeConcordanceDetailMetrics, indelCounter, TRUTH_SAMPLE, CALL_SAMPLE, MISSING_SITES_HOM_REF, OUTPUT_ALL_ROWS);
        genotypeConcordanceDetailMetrics.write(detailedMetricsFile);

        // Calculate and score the contingency metrics
        final MetricsFile<GenotypeConcordanceContingencyMetrics,?> genotypeConcordanceContingencyMetricsFile = getMetricsFile();
        GenotypeConcordanceContingencyMetrics contingencyMetrics = new GenotypeConcordanceContingencyMetrics(SNP, snpCounter, TRUTH_SAMPLE, CALL_SAMPLE, MISSING_SITES_HOM_REF);
        genotypeConcordanceContingencyMetricsFile.addMetric(contingencyMetrics);
        contingencyMetrics = new GenotypeConcordanceContingencyMetrics(INDEL, indelCounter, TRUTH_SAMPLE, CALL_SAMPLE, MISSING_SITES_HOM_REF);
        genotypeConcordanceContingencyMetricsFile.addMetric(contingencyMetrics);
        genotypeConcordanceContingencyMetricsFile.write(contingencyMetricsFile);

        for (final String condition : unClassifiedStatesMap.keySet()) {
            log.info("Uncovered truth/call Variant Context Type Counts: " + condition + " " + unClassifiedStatesMap.get(condition));
        }

        CloserUtil.close(callReader);
        CloserUtil.close(truthReader);
        writer.ifPresent(VariantContextWriter::close);

        return 0;
    }

    /** Gets the variant context writer if the output VCF is to be written, otherwise empty. */
    private Optional<VariantContextWriter> getVariantContextWriter(final VCFFileReader truthReader, final VCFFileReader callReader) {
        if (OUTPUT_VCF) {
            final File outputVcfFile = new File(OUTPUT + OUTPUT_VCF_FILE_EXTENSION);
            final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                    .setOutputFile(outputVcfFile)
                    .setReferenceDictionary(callReader.getFileHeader().getSequenceDictionary())
                    .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .setOption(Options.INDEX_ON_THE_FLY);
            final VariantContextWriter writer = builder.build();

            // create the output header
            final List<String> sampleNames = Arrays.asList(OUTPUT_VCF_CALL_SAMPLE_NAME, OUTPUT_VCF_TRUTH_SAMPLE_NAME);
            final Set<VCFHeaderLine> headerLines = new HashSet<>();
            headerLines.addAll(callReader.getFileHeader().getMetaDataInInputOrder());
            headerLines.addAll(truthReader.getFileHeader().getMetaDataInInputOrder());
            headerLines.add(CONTINGENCY_STATE_HEADER_LINE);
            writer.writeHeader(new VCFHeader(headerLines, sampleNames));
            return Optional.of(writer);
        }
        else {
            return Optional.empty();
        }
    }

    private void writeVcfTuple(final VcfTuple tuple, final VariantContextWriter writer, final GenotypeConcordanceScheme scheme) {
        VariantContext truthContext = null, callContext = null;
        final List<Genotype> genotypes = new ArrayList<>(2);

        // get the variant contexts and initialize the output variant context builder
        if (tuple.leftVariantContext.isPresent()) {
            truthContext = tuple.leftVariantContext.get();
        }
        if (tuple.rightVariantContext.isPresent()) {
            callContext = tuple.rightVariantContext.get();
        }

        //Don't write symbolic alleles to output VCF
        if (truthContext != null && truthContext.isSymbolic() || callContext != null && callContext.isSymbolic()) {
            return;
        }

        // Get the alleles for each genotype.  No alleles will be extracted for a genotype if the genotype is
        // mixed, filtered, or missing.
        final Alleles alleles = normalizeAlleles(truthContext, TRUTH_SAMPLE, callContext, CALL_SAMPLE, IGNORE_FILTER_STATUS);

        // There will be no alleles if both genotypes are one of mixed, filtered, or missing.  Do not output any
        // variant context in this case.
        if (!alleles.allAlleles.isEmpty()) {
            if (truthContext == null && callContext == null) {
                throw new IllegalStateException("Both truth and call contexts are null!");
            }

            final VariantContextBuilder builder;
            final List<Allele> allAlleles   = alleles.asList();
            final List<Allele> truthAlleles = alleles.truthAlleles();
            final List<Allele> callAlleles  = alleles.callAlleles();

            // Get the alleles present at this site for both samples to use for the output variant context, but remove no calls.
            final Set<Allele> siteAlleles = new HashSet<>();
            siteAlleles.addAll(allAlleles);
            siteAlleles.remove(Allele.NO_CALL);

            // Initialize the variant context builder
            final VariantContext initialContext = (callContext == null) ? truthContext : callContext;
            builder = new VariantContextBuilder(initialContext.getSource(), initialContext.getContig(), initialContext.getStart(), initialContext.getEnd(), siteAlleles);
            builder.computeEndFromAlleles(allAlleles, initialContext.getStart());
            builder.log10PError(initialContext.getLog10PError());

            // Add the genotypes
            addToGenotypes(genotypes, truthContext, TRUTH_SAMPLE, OUTPUT_VCF_TRUTH_SAMPLE_NAME, allAlleles, truthAlleles, MISSING_SITES_HOM_REF);
            addToGenotypes(genotypes, callContext, CALL_SAMPLE, OUTPUT_VCF_CALL_SAMPLE_NAME, allAlleles, callAlleles, false);

            // set the alleles and genotypes
            builder.genotypes(genotypes);

            // set the concordance state attribute
            final TruthAndCallStates state = GenotypeConcordance.determineState(truthContext, TRUTH_SAMPLE, callContext, CALL_SAMPLE, MIN_GQ, MIN_DP, IGNORE_FILTER_STATUS);
            final ContingencyState[] stateArray = scheme.getConcordanceStateArray(state.truthState, state.callState);
            builder.attribute(CONTINGENCY_STATE_TAG, Arrays.asList(stateArray));

            // write it
            writer.add(builder.make());
        }
    }

    /** Adds a new genotype to the provided list of genotypes.  If the given variant context is null or has no alleles
     * (ex. a no-call), then add a no-call, otherwise make a copy of the genotype for the sample with the input
     * sample name within the provided context.  In the former case, if missingSitesHomRef is true, a homozygous
     * reference genotype is created instead of a no-call.  In the latter case, the new genotype will have the output
     * sample name. */
    private void addToGenotypes(final List<Genotype> genotypes,
                                final VariantContext ctx,
                                final String inputSampleName,
                                final String outputSampleName,
                                final List<Allele> allAlleles,
                                final List<Allele> ctxAlleles,
                                final boolean missingSitesHomRef) {
        if (ctx != null && !ctxAlleles.isEmpty()) {
            final Genotype genotype = ctx.getGenotype(inputSampleName);
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            genotypeBuilder.name(outputSampleName);
            genotypeBuilder.alleles(ctxAlleles);
            if (!genotype.hasAnyAttribute(VCFConstants.GENOTYPE_KEY)) {
                genotypeBuilder.attribute(VCFConstants.GENOTYPE_KEY, MISSING_VALUE_v4);
            }
            genotypes.add(genotypeBuilder.make());
        }
        else {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(outputSampleName);
            if (missingSitesHomRef) {
                genotypeBuilder.alleles(Arrays.asList(allAlleles.get(0), allAlleles.get(0)));
            }
            else {
                genotypeBuilder.alleles(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL));
            }
            genotypes.add(genotypeBuilder.make());
        }
    }

    public static boolean classifyVariants(final Optional<VariantContext> truthContext,
                                           final String truthSample,
                                           final Optional<VariantContext> callContext,
                                           final String callSample,
                                           final int minGq, final int minDp,
                                           final boolean ignoreFilterStatus) {
        return classifyVariants(truthContext, truthSample, callContext, callSample, Optional.empty(), Optional.empty(),
                minGq, minDp, ignoreFilterStatus);
    }

    /**
     * Attempts to determine the concordance state given the truth and all variant context and optionally increments the genotype concordance
     * count for the given variant type (SNP or INDEL).  This will ignore cases where an indel was found in the truth and a SNP was found in
     * the call, and vice versa.  We typically fail to classify Mixed, Symbolic variants, or MNPs.
     *
     * @param truthContext A variant context representing truth
     * @param truthSample The name of the truth sample
     * @param callContext A variant context representing the call
     * @param callSample The name of the call sample
     * @param snpCounter optionally a place to increment the counts for SNP truth/call states
     * @param indelCounter optionally a place to increment the counts for INDEL truth/call states
     * @param minGq Threshold for filtering by genotype attribute GQ
     * @param minDp Threshold for filtering by genotype attribute DP
     * @return true if the concordance state could be classified.
     */
    public static boolean classifyVariants(final Optional<VariantContext> truthContext,
                                           final String truthSample,
                                           final Optional<VariantContext> callContext,
                                           final String callSample,
                                           final Optional<GenotypeConcordanceCounts> snpCounter,
                                           final Optional<GenotypeConcordanceCounts> indelCounter,
                                           final int minGq, final int minDp, final boolean ignoreFilteredStatus) {
        final VariantContext.Type truthVariantContextType = truthContext.map(VariantContext::getType).orElse(NO_VARIATION);
        final VariantContext.Type callVariantContextType  = callContext.map(VariantContext::getType).orElse(NO_VARIATION);

        // A flag to keep track of whether we have been able to successfully classify the Truth/Call States.
        // Unclassified include MIXED/MNP/Symbolic...
        final TruthAndCallStates truthAndCallStates = determineState(truthContext.orElse(null), truthSample, callContext.orElse(null), callSample, minGq, minDp, ignoreFilteredStatus);
        if (truthVariantContextType == SNP) {
            if ((callVariantContextType == SNP) || (callVariantContextType == MIXED) || (callVariantContextType == NO_VARIATION)) {
                // Note.  If truth is SNP and call is MIXED, the event will be logged in the snpCounter, with row = MIXED
                snpCounter.ifPresent(counter -> counter.increment(truthAndCallStates));
                return true;
            }
        }
        else if (truthVariantContextType == INDEL) {
            // Note.  If truth is Indel and call is MIXED, the event will be logged in the indelCounter, with row = MIXED
            if ((callVariantContextType == INDEL) || (callVariantContextType == MIXED) || (callVariantContextType == NO_VARIATION)) {
                indelCounter.ifPresent(counter -> counter.increment(truthAndCallStates));
                return true;
            }
        }
        else if (truthVariantContextType == MIXED) {
            // Note.  If truth is MIXED and call is SNP, the event will be logged in the snpCounter, with column = MIXED
            if (callVariantContextType == SNP) {
                snpCounter.ifPresent(counter -> counter.increment(truthAndCallStates));
                return true;
            }
            // Note.  If truth is MIXED and call is INDEL, the event will be logged in the indelCounter, with column = MIXED
            else if (callVariantContextType == INDEL) {
                indelCounter.ifPresent(counter -> counter.increment(truthAndCallStates));
                return true;
            }
        }
        else if (truthVariantContextType == NO_VARIATION) {
            if (callVariantContextType == SNP) {
                snpCounter.ifPresent(counter -> counter.increment(truthAndCallStates));
                return true;
            }
            else if (callVariantContextType == INDEL) {
                indelCounter.ifPresent(counter -> counter.increment(truthAndCallStates));
                return true;
            }
        }
        return false;
    }


    /**
     * Method to add missing sites that are KNOWN to be HOM_REF in the case of the NIST truth data set.
     */
    public static void addMissingTruthAndMissingCallStates(final double numVariants, final long intervalBaseCount, final GenotypeConcordanceCounts counter) {
        final double countMissingMissing = intervalBaseCount-numVariants;
        final TruthAndCallStates missingMissing = new TruthAndCallStates(TruthState.MISSING, CallState.MISSING);
        counter.increment(missingMissing, countMissingMissing);
    }

    /**
    * Outputs the detailed statistics tables for SNP and Indel match categories.
    **/
    public static void outputDetailMetricsFile(final VariantContext.Type variantType, final MetricsFile<GenotypeConcordanceDetailMetrics,?> genotypeConcordanceDetailMetricsFile,
                                               final GenotypeConcordanceCounts counter, final String truthSampleName, final String callSampleName,
                                               final boolean missingSitesHomRef, final boolean outputAllRows) {
        final GenotypeConcordanceSchemeFactory schemeFactory = new GenotypeConcordanceSchemeFactory();
        final GenotypeConcordanceScheme scheme = schemeFactory.getScheme(missingSitesHomRef);
        scheme.validateScheme();
        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                final long count = counter.getCount(truthState, callState);
                final String contingencyValues = scheme.getContingencyStateString(truthState, callState);
                if (count > 0 || outputAllRows) {
                    final GenotypeConcordanceDetailMetrics detailMetrics = new GenotypeConcordanceDetailMetrics();
                    detailMetrics.VARIANT_TYPE = variantType;
                    detailMetrics.TRUTH_SAMPLE = truthSampleName;
                    detailMetrics.CALL_SAMPLE = callSampleName;
                    detailMetrics.TRUTH_STATE = truthState;
                    detailMetrics.CALL_STATE = callState;
                    detailMetrics.COUNT = count;
                    detailMetrics.CONTINGENCY_VALUES = contingencyValues;
                    genotypeConcordanceDetailMetricsFile.addMetric(detailMetrics);
                }
            }
        }
    }

    /** A simple structure to return the results of getAlleles.  NB: truthAllele1/2 and callAllele1/2 may be null if
     * no first or second allele was present respectively. */
    protected static class Alleles {
        public final OrderedSet<String> allAlleles;
        public final String truthAllele1;
        public final String truthAllele2;
        public final String callAllele1;
        public final String callAllele2;

        public Alleles(final OrderedSet<String> allAlleles,
                       final String truthAllele1,
                       final String truthAllele2,
                       final String callAllele1,
                       final String callAllele2) {

            if (truthAllele1 == null && truthAllele2 != null) {
                throw new IllegalStateException("truthAllele2 should be null if truthAllele1 is null.");
            }
            if (callAllele1 == null && callAllele2 != null) {
                throw new IllegalStateException("callAllele2 should be null if callAllele1 is null.");
            }

            this.allAlleles   = allAlleles;
            this.truthAllele1 = (truthAllele1 == null) ? null : this.allAlleles.get(allAlleles.indexOf(truthAllele1));
            this.truthAllele2 = (truthAllele2 == null) ? null : this.allAlleles.get(allAlleles.indexOf(truthAllele2));
            this.callAllele1  = (callAllele1 == null) ? null : this.allAlleles.get(allAlleles.indexOf(callAllele1));
            this.callAllele2  = (callAllele2 == null) ? null : this.allAlleles.get(allAlleles.indexOf(callAllele2));
        }

        public List<Allele> asList() {
            if (allAlleles.isEmpty()) {
                return Collections.emptyList();
            }
            else {
                final List<Allele> alleles = new ArrayList<>(this.allAlleles.size());
                for (int idx = 0; idx < this.allAlleles.size(); idx++) {
                    alleles.add(Allele.create(this.allAlleles.get(idx), idx == 0));
                }
                return alleles;
            }
        }

        public List<Allele> truthAlleles() {
            if (allAlleles.isEmpty() || this.truthAllele1 == null) {
                return Collections.emptyList();
            }
            else {
                final Allele truthAllele1 = Allele.create(this.truthAllele1, this.allAlleles.indexOf(this.truthAllele1) == 0);
                final Allele truthAllele2 = Allele.create(this.truthAllele2, this.allAlleles.indexOf(this.truthAllele2) == 0);
                return Arrays.asList(truthAllele1, truthAllele2);
            }
        }

        public List<Allele> callAlleles() {
            if (allAlleles.isEmpty() || this.callAllele1 == null) {
                return Collections.emptyList();
            }
            else {
                final Allele callAllele1 = Allele.create(this.callAllele1, this.allAlleles.indexOf(this.callAllele1) == 0);
                final Allele callAllele2 = Allele.create(this.callAllele2, this.allAlleles.indexOf(this.callAllele2) == 0);
                return Arrays.asList(callAllele1, callAllele2);
            }
        }
    }

    /** Inserts the given string into the destination string at the given index.  If the index is past the end of the
     * destination string, the given string is appended to the destination.  If the destination string is the
     * spanning deletion allele it will be returned unchanged.
     */
    static String spliceOrAppendString(final String destination, final String toInsert, final int insertIdx) {
        if (destination.equals(Allele.SPAN_DEL_STRING)) {
            return destination;
        }
        if (insertIdx <= destination.length()) {
            return destination.substring(0, insertIdx) + toInsert + destination.substring(insertIdx);
        }
        return destination + toInsert;
    }

    /** Gets the alleles for the truth and call genotypes.  In particular, this handles the case where indels can have different
     * reference alleles. */
    final protected static Alleles normalizeAlleles(final VariantContext truthContext,
                                                    final String truthSample,
                                                    final VariantContext callContext,
                                                    final String callSample,
                                                    final Boolean ignoreFilteredStatus) {

        final Genotype truthGenotype, callGenotype;

        if (truthContext == null || truthContext.isMixed() || truthContext.isFiltered()) truthGenotype = null;
        else truthGenotype = truthContext.getGenotype(truthSample);

        if (callContext == null || callContext.isMixed() || (!ignoreFilteredStatus && callContext.isFiltered())) callGenotype = null;
        else callGenotype = callContext.getGenotype(callSample);

        // initialize the reference
        String truthRef = (truthGenotype != null) ? truthContext.getReference().getBaseString() : null;
        String callRef  = (callGenotype != null) ? callContext.getReference().getBaseString() : null;

        String truthAllele1 = null;
        String truthAllele2 = null;
        if (null != truthGenotype) {
            if (truthGenotype.getAlleles().size() != 2) {
                throw new IllegalStateException("Genotype for Variant Context: " + truthContext + " does not have exactly 2 alleles");
            }
            truthAllele1 = truthGenotype.getAllele(0).getBaseString();
            truthAllele2 = truthGenotype.getAllele(1).getBaseString();
        }

        String callAllele1 = null;
        String callAllele2 = null;
        if (null != callGenotype) {
            if (callGenotype.getAlleles().size() != 2) {
                throw new IllegalStateException("Genotype for Variant Context: " + callContext + " does not have exactly 2 alleles");
            }
            callAllele1 = callGenotype.getAllele(0).getBaseString();
            callAllele2 = callGenotype.getAllele(1).getBaseString();
        }

        if ((truthRef != null && callRef != null) && (!truthRef.equals(callRef))) {
            // We must handle different representations for indels based on their reference alleles.  For example:
            // - Deletion:  TCAA*/T versus TCAACAA*/TCAA
            // - Insertion: TCAA*/TCAAGG versus TCAACAA*/TCAACAAGG
            // We do the following:
            // 1. Verify that the shorter reference allele is a prefix of the longer reference allele.
            // 2. Add the extra reference bases to the variant allele.
            // 3. Update the shorter reference allele to be the longer reference allele.
            if (truthRef.length() < callRef.length()) {
                // Truth reference is shorter than call reference
                final String suffix = getStringSuffix(callRef, truthRef, "Ref alleles mismatch between: " + truthContext + " and " + callContext);
                final int insertIdx = truthRef.length();
                truthAllele1 = truthAllele1.equals(Allele.NO_CALL_STRING) ? truthAllele1 : spliceOrAppendString(truthAllele1, suffix, insertIdx);
                truthAllele2 = truthAllele2.equals(Allele.NO_CALL_STRING) ? truthAllele2 : spliceOrAppendString(truthAllele2, suffix, insertIdx);
                truthRef = truthRef + suffix;

            }
            else if (truthRef.length() > callRef.length()) {
                // call reference is shorter than truth:
                final String suffix = getStringSuffix(truthRef, callRef, "Ref alleles mismatch between: " + truthContext + " and " + callContext);
                final int insertIdx = callRef.length();
                callAllele1 = callAllele1.equals(Allele.NO_CALL_STRING) ? callAllele1 : spliceOrAppendString(callAllele1, suffix, insertIdx);
                callAllele2 = callAllele2.equals(Allele.NO_CALL_STRING) ? callAllele2 : spliceOrAppendString(callAllele2, suffix, insertIdx);
                callRef = callRef + suffix;
            }
            else {
                // Same length - so they must just disagree...
                throw new IllegalStateException("Ref alleles mismatch between: " + truthContext + " and " + callContext);
            }
        }

        final OrderedSet<String> allAlleles = new OrderedSet<String>();

        if (truthGenotype != null || callGenotype != null) {
            // Store the refAllele as the first (0th index) allele in allAlleles (only can do if at least one genotype is non-null)
            allAlleles.smartAdd(truthGenotype == null ? callRef : truthRef); // zeroth allele;
        }

        if (null != truthGenotype) {
            allAlleles.smartAdd(truthAllele1);
            allAlleles.smartAdd(truthAllele2);
        }

        /**
         *  if either of the call alleles is in allAlleles, with index > 1, we need to make sure that allele has index 1.
         *  this is because of the following situations:
         *
         *      REF TRUTH   CALL-GT TRUTH-STATE     CALL-STATE
         *      A   C/G     C/A     HET_VAR1_VAR2   HET_REF_VAR1
         *      A   G/C     C/A     HET_VAR1_VAR2   HET_REF_VAR1
         *      A   G/C     G/A     HET_VAR1_VAR2   HET_REF_VAR1
         *      A   G/C     G/A     HET_VAR1_VAR2   HET_REF_VAR1
         *
         *  so, in effect, the order of the alleles in the TRUTH doesn't determine the assignment of allele to Var1 and Var2,
         *  only once the call is known can this assignment be made.
         */

        if (null != callGenotype) {
            if (allAlleles.indexOf(callAllele1) > 1 || allAlleles.indexOf(callAllele2) > 1) {
                allAlleles.remove(2);
                allAlleles.remove(1);
                allAlleles.smartAdd(truthAllele2);
                allAlleles.smartAdd(truthAllele1);
            }

            allAlleles.smartAdd(callAllele1);
            allAlleles.smartAdd(callAllele2);
        }

        return new Alleles(allAlleles, truthAllele1, truthAllele2, callAllele1, callAllele2);
    }

    private static GenotypeConcordanceStateCodes getStateCode(final VariantContext ctx, final String sample, final int minGq, final int minDp) {
        // Site level checks
        // Added potential to include missing sites as hom ref.
        if (ctx == null) return MISSING_CODE;
        else if (ctx.isMixed()) return IS_MIXED_CODE;
        else if (ctx.isFiltered()) return VC_FILTERED_CODE;
        else {
            // Genotype level checks
            final Genotype genotype = ctx.getGenotype(sample);
            if (genotype.isNoCall())           return NO_CALL_CODE;
            else if (genotype.isFiltered())    return GT_FILTERED_CODE;
            else if ((genotype.getGQ() != -1) && (genotype.getGQ() < minGq)) return LOW_GQ_CODE;
            else if ((genotype.getDP() != -1) && (genotype.getDP() < minDp)) return LOW_DP_CODE;
                // Note.  Genotype.isMixed means that it is called on one chromosome and NOT on the other
            else if ((genotype.isMixed())) return NO_CALL_CODE;
        }
        return null;
    }

    /**
     * A method to determine the truth and call states for a pair of variant contexts representing truth and call.
     * A variety of variant and genotype-level checks are first used to determine if either of the the variant contexts
     * are filtered and after that a comparison of the called genotype alleles to determine appropriate truth and call state
     *
     * Note that this method does NOT check for SNP versus Indel.  It is assumed that that check is done by the caller and the results
     * of this method are stored by SNP/Indel.
     * Note that if a variant context has BOTH GQ and DP less than the specified threshold, then it will be of Truth/Call State LOW_GQ
     *
     * @param truthContext A variant context representing truth
     * @param truthSample The name of the truth sample
     * @param callContext A variant context representing the call
     * @param callSample The name of the call sample
     * @param minGq Threshold for filtering by genotype attribute GQ
     * @param minDp Threshold for filtering by genotype attribute DP
     * @return TruthAndCallStates object containing the TruthState and CallState determined here.
     */
    final public static TruthAndCallStates determineState(final VariantContext truthContext, final String truthSample, final VariantContext callContext, final String callSample, final int minGq, final int minDp, final Boolean ignoreFilteredStatus) {
        TruthState truthState = null;
        CallState callState = null;

        // TODO: what about getPloidy()

        // Get truth and call states if they are filtered or are not going to be compared (ex. depth is less than minDP).
        final GenotypeConcordanceStateCodes truthStateCode = getStateCode(truthContext, truthSample, minGq, minDp);

        if (null != truthStateCode) {
            truthState = GenotypeConcordanceStates.truthMap.get(truthStateCode.ordinal());
        }
        GenotypeConcordanceStateCodes callStateCode = getStateCode(callContext, callSample, minGq, minDp);
        if (ignoreFilteredStatus && callStateCode == GenotypeConcordanceStateCodes.VC_FILTERED_CODE) {
            callStateCode = null;
        }
        if (null != callStateCode) {
            callState = GenotypeConcordanceStates.callMap.get(callStateCode.ordinal());
        }

        final Alleles alleles = normalizeAlleles(
                truthState == null ? truthContext : null,
                truthSample,
                callState == null ? callContext : null,
                callSample, ignoreFilteredStatus);
        final OrderedSet<String> allAlleles = alleles.allAlleles;
        final String truthAllele1           = alleles.truthAllele1;
        final String truthAllele2           = alleles.truthAllele2;
        final String callAllele1            = alleles.callAllele1;
        final String callAllele2            = alleles.callAllele2;

        // Truth
        if (null == truthState) {
            final int allele0idx = allAlleles.indexOf(truthAllele1);
            final int allele1idx = allAlleles.indexOf(truthAllele2);

            if (allele0idx == allele1idx) { //HOM
                truthState = TruthState.getHom(allele0idx);
            } else { //HET
                truthState = TruthState.getVar(allele0idx, allele1idx);
            }
        }

        // Call
        if (null == callState) {
            final int allele0idx = allAlleles.indexOf(callAllele1);
            final int allele1idx = allAlleles.indexOf(callAllele2);

            if (allele0idx == allele1idx) { //HOM
                callState = CallState.getHom(allele0idx);
            } else { //HET
                callState = CallState.getHet(allele0idx, allele1idx);
            }

            if (null == callState) {
                throw new IllegalStateException("This should never happen...  Could not classify the call variant: " + callContext.getGenotype(callSample));
            }
        }

        return new TruthAndCallStates(truthState, callState);
    }

    final static String getStringSuffix(final String longerString, final String shorterString, final String errorMsg) {
        // Truth reference is shorter than call reference
        if (!longerString.startsWith(shorterString)) {
            throw new IllegalStateException(errorMsg);
        }
        return longerString.substring(shorterString.length());
    }
}

/** like a list, but if you ask for an index of an item, it will first add that item.
also, same item cannot be added more than once (like a set)
 */
class OrderedSet<T> extends ArrayList<T> {

    public int smartIndexOf(final T o) {
        smartAdd(o);
        return super.indexOf(o);
    }

    public boolean smartAdd(final T o) {
        if (!this.contains(o)) {
            return add(o);
        }
        return false;
    }
}




