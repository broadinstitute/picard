
/*
 * The MIT License
 *
 * Copyright (c) 2010, 2016 The Broad Institute
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

package picard.fingerprint;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Fingerprinting;
import picard.fingerprint.CrosscheckMetric.FingerprintResult;

import java.io.*;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Checks that all data in the set of input files appear to come from the same
 * individual. Can be used to compare according to readgroups, libraries, samples, or files.
 * Operates on bams/sams and vcfs (including gvcfs).
 *
 * <h3>Summary</h3>
 * Checks if all the genetic data within a set of files appear to come from the same individual.
 * It quickly determines whether a "group's" genotype matches that of an input SAM/BAM/VCF by selective sampling,
 * and has been designed to work well even for low-depth SAM/BAMs.
 * <br/>
 * The tool collects "fingerprints" (essentially genotype information from different parts of the genome)
 * at the finest level available in the data (readgroup for SAM files
 * and sample for VCF files) and then optionally aggregates it by library, sample or file, to increase power and provide
 * results at the desired resolution. Output is in a "Moltenized" format, one row per comparison. The results will
 * be emitted into a metric file for the class {@link CrosscheckMetric}.
 * In this format the output will include the LOD score and also tumor-aware LOD score which can
 * help assess identity even in the presence of a severe loss of heterozygosity with high purity (which could
 * otherwise fail to notice that samples are from the same individual.)
 * A matrix output is also available to facilitate visual inspection of crosscheck results.
 * <br/>
 * Since there can be many rows of output in the metric file, we recommend the use of {@link ClusterCrosscheckMetrics}
 * as a follow-up step to running CrosscheckFingerprints.
 * <br/>
 * There are cases where one would like to identify a few groups out of a collection of many possible groups (say
 * to link a bam to it's correct sample in a multi-sample vcf. In this case one would not case for the cross-checking
 * of the various samples in the VCF against each other, but only in checking the identity of the bam against the various
 * samples in the vcf. The {@link #SECOND_INPUT} is provided for this use-case. With {@link #SECOND_INPUT} provided, CrosscheckFingerprints
 * does the following:
 * <il>
 * <li>aggregation of data happens independently for the input files in {@link #INPUT} and {@link #SECOND_INPUT}.</li>
 * <li>aggregation of data happens at the SAMPLE level.</li>
 * <li>each samples from {@link #INPUT} will only be compared to that same sample in {@link #INPUT}.</li>
 * <li>{@link #MATRIX_OUTPUT} is disabled.</li>
 * </il>
 * <br/>
 * <h3>Examples</h3>
 * <h4>Check that all the readgroups from a sample match each other:</h4>
 * <pre>
 *     java -jar picard.jar CrosscheckFingerprints \
 *          INPUT=sample.with.many.readgroups.bam \
 *          HAPLOTYPE_DATABASE=fingerprinting_haplotype_database.txt \
 *          LOD_THRESHOLD=-5 \
 *          OUTPUT=sample.crosscheck_metrics
 * </pre>
 *
 *
 * <h4>Check that all the readgroups match as expected when providing reads from two samples from the same individual:</h4>
 * <pre>
 *     java -jar picard.jar CrosscheckFingerprints \
 *          INPUT=sample.one.with.many.readgroups.bam \
 *          INPUT=sample.two.with.many.readgroups.bam \
 *          HAPLOTYPE_DATABASE=fingerprinting_haplotype_database.txt \
 *          LOD_THRESHOLD=-5 \
 *          EXPECT_ALL_GROUPS_TO_MATCH=true \
 *          OUTPUT=sample.crosscheck_metrics
 * </pre>
 *
 *
 * <h3> Detailed Explanation</h3>
 *
 * This tool calculates the LOD score for identity check between "groups" of data in the INPUT files as defined by
 * the CROSSCHECK_BY argument. A positive value indicates that the data seems to have come from the same individual
 * or, in other words the identity checks out. The scale is logarithmic (base 10), so a LOD of 6 indicates
 * that it is 1,000,000 more likely that the data matches the genotypes than not. A negative value indicates
 * that the data do not match. A score that is near zero is inconclusive and can result from low coverage
 * or non-informative genotypes. Each group is assigned a sample identifier (for SAM this is taken from the SM tag in
 * the appropriate readgroup header line, for VCF this is taken from the column label in the file-header.
 * After combining all the data from the same "group" together, an all-against-all comparison is performed. Results are
 * categorized a {@link FingerprintResult} enum: EXPECTED_MATCH, EXPECTED_MISMATCH, UNEXPECTED_MATCH, UNEXPECTED_MISMATCH,
 * or AMBIGUOUS depending on the LOD score and on whether the sample identifiers of the groups agree: LOD scores that are
 * less than LOD_THRESHOLD are considered mismatches, and those greater than -LOD_THRESHOLD are matches (between is ambiguous).
 * If the sample identifiers are equal, the groups are expected to match. They are expected to mismatch otherwise.
 * <br/>
 *
 * The identity check makes use of haplotype blocks defined in the HAPLOTYPE_MAP file to enable it to have higher
 * statistical power for detecting identity or swap by aggregating data from several SNPs in the haplotype block. This
 * enables an identity check of samples with very low coverage (e.g. ~1x mean coverage).
 * <br/>
 * When provided a VCF, the identity check looks at the PL, GL and GT fields (in that order) and uses the first one that
 * it finds.
 *
 *
 *
 * @author Tim Fennell
 * @author Yossi Farjoun
 */


@CommandLineProgramProperties(
        summary =
                "Checks that all data in the set of input files appear to come from the same " +
                        "individual. Can be used to cross-check readgroups, libraries, samples, or files. " +
                        "Operates on bams/sams and vcfs (including gvcfs). " +
                        "\n" +
                        "<h3>Summary</h3>\n" +
                        "Checks if all the genetic data within a set of files appear to come from the same individual. " +
                        "It quickly determines whether a group's genotype matches that of an input SAM/BAM/VCF by selective sampling, " +
                        "and has been designed to work well for low-depth SAM/BAMs (as well as high depth ones and VCFs.) " +
                        "The tool collects fingerprints (essentially, genotype information from different parts of the genome)" +
                        " at the finest level available in the data (readgroup for SAM files " +
                        "and sample for VCF files) and then optionally aggregates it by library, sample or file, to increase power and provide " +
                        "results at the desired resolution. Output is in a \"Moltenized\" format, one row per comparison. The results are " +
                        "emitted into a CrosscheckMetric metric file. " +
                        "In this format the output will include the LOD score and also tumor-aware LOD score which can " +
                        "help assess identity even in the presence of a severe loss of heterozygosity with high purity (which could cause it to " +
                        "otherwise fail to notice that samples are from the same individual.) " +
                        "A matrix output is also available to facilitate visual inspection of crosscheck results.\n " +
                        "\n" +
                        "Since there can be many rows of output in the metric file, we recommend the use of ClusterCrosscheckMetrics " +
                        "as a follow-up step to running CrosscheckFingerprints.\n " +
                        "\n" +
                        "There are cases where one would like to identify a few groups out of a collection of many possible groups (say " +
                        "to link a bam to it's correct sample in a multi-sample vcf. In this case one would not case for the cross-checking " +
                        "of the various samples in the VCF against each other, but only in checking the identity of the bam against the various " +
                        "samples in the vcf. The SECOND_INPUT is provided for this use-case. With SECOND_INPUT provided, CrosscheckFingerprints " +
                        "does the following:\n" +
                        " - aggregation of data happens independently for the input files in INPUT and SECOND_INPUT. \n" +
                        " - aggregation of data happens at the SAMPLE level \n" +
                        " - each samples from INPUT will only be compared to that same sample in SECOND_INPUT. \n" +
                        " - MATRIX_OUTPUT is disabled. " +
                        "\n" +
                        "<hr/>" +
                        "<h3>Examples</h3>" +
                        "<h4>Check that all the readgroups from a sample match each other:</h4>" +
                        "<pre>" +
                        "    java -jar picard.jar CrosscheckFingerprints \\\n" +
                        "          INPUT=sample.with.many.readgroups.bam \\\n" +
                        "          HAPLOTYPE_DATABASE=fingerprinting_haplotype_database.txt \\\n" +
                        "          LOD_THRESHOLD=-5 \\\n" +
                        "          OUTPUT=sample.crosscheck_metrics" +
                        " </pre>" +
                        "\n" +
                        " <h4>Check that all the readgroups match as expected when providing reads from two samples from the same individual:</h4>" +
                        " <pre>" +
                        "     java -jar picard.jar CrosscheckFingerprints \\\n" +
                        "           INPUT=sample.one.with.many.readgroups.bam \\\n" +
                        "           INPUT=sample.two.with.many.readgroups.bam \\\n" +
                        "           HAPLOTYPE_DATABASE=fingerprinting_haplotype_database.txt \\\n" +
                        "           LOD_THRESHOLD=-5 \\\n" +
                        "           EXPECT_ALL_GROUPS_TO_MATCH=true \\\n" +
                        "           OUTPUT=sample.crosscheck_metrics" +
                        " </pre>" +
                        "\n" +
                        "\n" +
                        "<h4>Detailed Explanation</h4>" +
                        "\n" +
                        "This tool calculates the LOD score for identity check between \"groups\" of data in the INPUT files as defined by " +
                        "the CROSSCHECK_BY argument. A positive value indicates that the data seems to have come from the same individual " +
                        "or, in other words the identity checks out. The scale is logarithmic (base 10), so a LOD of 6 indicates " +
                        "that it is 1,000,000 more likely that the data matches the genotypes than not. A negative value indicates " +
                        "that the data do not match. A score that is near zero is inconclusive and can result from low coverage " +
                        "or non-informative genotypes. Each group is assigned a sample identifier (for SAM this is taken from the SM tag in " +
                        "the appropriate readgroup header line, for VCF this is taken from the column label in the file-header. " +
                        "After combining all the data from the same group together, an all-against-all comparison is performed. Results are " +
                        "categorized as one of EXPECTED_MATCH, EXPECTED_MISMATCH, UNEXPECTED_MATCH, UNEXPECTED_MISMATCH, or AMBIGUOUS depending " +
                        "on the LOD score and on whether the sample identifiers of the groups agree: LOD scores that are less than LOD_THRESHOLD " +
                        "are considered mismatches, and those greater than -LOD_THRESHOLD are matches (between is ambiguous). " +
                        "If the sample identifiers are equal, the groups are expected to match. They are expected to mismatch otherwise. " +
                        "\n" +
                        "\n" +
                        "The identity check makes use of haplotype blocks defined in the HAPLOTYPE_MAP file to enable it to have higher " +
                        "statistical power for detecting identity or swap by aggregating data from several SNPs in the haplotype block. This " +
                        "enables an identity check of samples with very low coverage (e.g. ~1x mean coverage).\n " +
                        "\n" +
                        "When provided a VCF, the identity check looks at the PL, GL and GT fields (in that order) and uses the first one that " +
                        "it finds. ",
        oneLineSummary = "Checks that all data in the input files appear to have come from the same individual",
        programGroup = Fingerprinting.class
)
@DocumentedFeature
public class CrosscheckFingerprints extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input files (or lists of files) with which to compare fingerprints.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = "SI", optional = true, mutex={"MATRIX_OUTPUT"},
            doc = "A second set of input files (or lists of files) with which to compare fingerprints. If this option is provided " +
                    "the tool compares each sample in INPUT with the sample from SECOND_INPUT that has the same sample ID." +
                    " In addition, data will be grouped by SAMPLE regardless of the value of CROSSCHECK_BY."+
                    " When operating in this mode, each sample in INPUT must also have a corresponding sample in SECOND_INPUT. " +
                    "If this is violated, the tool will proceed to check the matching samples, but report the missing samples" +
                    " and return a non-zero error-code." )
    public List<File> SECOND_INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true,
            doc = "Optional output file to write metrics to. Default is to write to stdout.")
    public File OUTPUT = null;

    @Argument(shortName = "MO", optional = true,
            doc = "Optional output file to write matrix of LOD scores to. This is less informative than the metrics output " +
                    "and only contains Normal-Normal LOD score (i.e. doesn't account for Loss of Heterozygosity). " +
                    "It is however sometimes easier to use visually.", mutex = {"SECOND_INPUT"})
    public File MATRIX_OUTPUT = null;

    @Argument(shortName = "H", doc = "The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Argument(shortName = "LOD",
            doc = "If any two groups (with the same sample name) match with a LOD score lower than the threshold " +
                    "the tool will exit with a non-zero code to indicate error." +
                    " Program will also exit with an error if it finds two groups with different sample name that " +
                    "match with a LOD score greater than -LOD_THRESHOLD." +
                    "\n\n" +
                    "LOD score 0 means equal likelihood " +
                    "that the groups match vs. come from different individuals, negative LOD score -N, mean 10^N time more likely " +
                    "that the groups are from different individuals, and +N means 10^N times more likely that " +
                    "the groups are from the same individual. ")
    public double LOD_THRESHOLD = 0;

    @Argument(doc = "Specificies which data-type should be used as the basic comparison unit. Fingerprints from readgroups can " +
            "be \"rolled-up\" to the LIBRARY, SAMPLE, or FILE level before being compared." +
            " Fingerprints from VCF can be be compared by SAMPLE or FILE.")
    public CrosscheckMetric.DataType CROSSCHECK_BY = CrosscheckMetric.DataType.READGROUP;

    @Argument(doc = "The number of threads to use to process files and generate fingerprints.")
    public int NUM_THREADS = 1;

    @Argument(doc = "Allow the use of duplicate reads in performing the comparison. Can be useful when duplicate " +
            "marking has been overly aggressive and coverage is low.")
    public boolean ALLOW_DUPLICATE_READS = false;

    @Argument(doc = "Assumed genotyping error rate that provides a floor on the probability that a genotype comes from " +
            "the expected sample. Must be greater than zero. " )
    public double GENOTYPING_ERROR_RATE = 0.01;

    @Argument(doc = "If true then only groups that do not relate to each other as expected will have their LODs reported.")
    public boolean OUTPUT_ERRORS_ONLY = false;

    @Argument(doc = "The rate at which a heterozygous genotype in a normal sample turns into a homozygous (via loss of heterozygosity) " +
            "in the tumor (model assumes independent events, so this needs to be larger than reality).", optional = true)
    public double LOSS_OF_HET_RATE = 0.5;

    @Argument(doc = "Expect all groups' fingerprints to match, irrespective of their sample names.  By default (with this value set to " +
            "false), groups (readgroups, libraries, files, or samples) with different sample names are expected to mismatch, and those with the " +
            "same sample name are expected to match. ")
    public boolean EXPECT_ALL_GROUPS_TO_MATCH = false;

    @Argument(doc = "When one or more mismatches between groups is detected, exit with this value instead of 0.")
    public int EXIT_CODE_WHEN_MISMATCH = 1;

    private final Log log = Log.getInstance(CrosscheckFingerprints.class);

    private double[][] crosscheckMatrix = null;
    private final List<String> matrixKeys = new ArrayList<>();

    @Override
    protected String[] customCommandLineValidation() {
        if (GENOTYPING_ERROR_RATE <= 0 ) {
            return new String[]{"GENOTYPING_ERROR_RATE must be greater than zero. Found " + GENOTYPING_ERROR_RATE};
        }
        return super.customCommandLineValidation();
    }

    @Override
    protected int doWork() {
        // Check inputs
        IOUtil.assertFilesAreReadable(INPUT);
        IOUtil.assertFilesAreReadable(SECOND_INPUT);

        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        if (OUTPUT != null) IOUtil.assertFileIsWritable(OUTPUT);

        if (!SECOND_INPUT.isEmpty()) {
            log.info("SECOND_INPUT is not empty. NOT doing cross-check. Will only compare each SAMPLE in INPUT against that sample in SECOND_INPUT.");
            if (CROSSCHECK_BY != CrosscheckMetric.DataType.SAMPLE) {
                log.warn("CROSSCHECK_BY is not SAMPLE, This doesn't make sense in non-crosscheck mode. Setting CROSSCHECK_BY to SAMPLE.");
                CROSSCHECK_BY = CrosscheckMetric.DataType.SAMPLE;
            }
        }

        if (MATRIX_OUTPUT != null) IOUtil.assertFileIsWritable(MATRIX_OUTPUT);

        final HaplotypeMap map = new HaplotypeMap(HAPLOTYPE_MAP);
        final FingerprintChecker checker = new FingerprintChecker(map);

        checker.setAllowDuplicateReads(ALLOW_DUPLICATE_READS);
        checker.setValidationStringency(VALIDATION_STRINGENCY);


        final List<String> extensions = new ArrayList<>();

        extensions.add(BamFileIoUtils.BAM_FILE_EXTENSION);
        extensions.add(IOUtil.SAM_FILE_EXTENSION);
        extensions.addAll(Arrays.asList(IOUtil.VCF_EXTENSIONS));

        final List<File> unrolledFiles = IOUtil.unrollFiles(INPUT, extensions.toArray(new String[extensions.size()]));
        IOUtil.assertFilesAreReadable(unrolledFiles);

        // unroll and check readable here, as it can be annoying to fingerprint INPUT files and only then discover a problem
        // in a file in SECOND_INPUT
        final List<File> unrolledFiles2 = IOUtil.unrollFiles(SECOND_INPUT, extensions.toArray(new String[extensions.size()]));
        IOUtil.assertFilesAreReadable(unrolledFiles2);

        log.info("Fingerprinting " + unrolledFiles.size() + " INPUT files.");
        final Map<FingerprintIdDetails, Fingerprint> fpMap = checker.fingerprintFiles(unrolledFiles, NUM_THREADS, 1, TimeUnit.DAYS);

        final List<CrosscheckMetric> metrics = new ArrayList<>();
        final int numUnexpected;

        if (SECOND_INPUT.isEmpty()) {
            log.info("Cross-checking all " + CROSSCHECK_BY + " against each other");
            numUnexpected = crossCheckGrouped(fpMap, metrics, getFingerprintIdDetailsStringFunction(CROSSCHECK_BY), CROSSCHECK_BY);
        } else {
            log.info("Fingerprinting " + unrolledFiles2.size() + " SECOND_INPUT files.");

            final Map<FingerprintIdDetails, Fingerprint> fpMap2 = checker.fingerprintFiles(unrolledFiles2, NUM_THREADS, 1, TimeUnit.DAYS);
            log.info("Checking each sample in INPUT with the same sample in SECOND_INPUT.");
            numUnexpected = checkFingerprintsBySample(fpMap, fpMap2, metrics);
        }

        final MetricsFile<CrosscheckMetric, ?> metricsFile = getMetricsFile();
        metricsFile.addAllMetrics(metrics);
        if (OUTPUT != null) {
            metricsFile.write(OUTPUT);
        } else {
            metricsFile.write(new OutputStreamWriter(System.out));
        }

        if (MATRIX_OUTPUT != null) {
            writeMatrix();
        }

        if (numUnexpected > 0) {
            log.warn("At least two read groups did not relate as expected.");
            return EXIT_CODE_WHEN_MISMATCH;
        } else {
            log.info("All read groups related as expected.");
            return 0;
        }
    }

    private void writeMatrix() {

        final NumberFormat format = NumberFormat.getInstance();
        format.setMaximumFractionDigits(4);

        try (final OutputStream stream = new FileOutputStream(MATRIX_OUTPUT);
             final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stream))) {

            // write the type by which the roll-up happened in the top left corner of the matrix
            writer.write(CROSSCHECK_BY.name());

            // write the nams of the keys as the first row
            for (int i = 0; i < crosscheckMatrix.length; i++) {
                writer.write('\t' + matrixKeys.get(i));
            }
            writer.newLine();

            for (int i = 0; i < crosscheckMatrix.length; i++) {
                // write the key in the first column
                writer.write(matrixKeys.get(i));
                // and then write all the values
                for (int j = 0; j < crosscheckMatrix.length; j++) {
                    writer.write('\t' + format.format(crosscheckMatrix[i][j]));
                }
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static Function<FingerprintIdDetails, String> getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType CROSSCHECK_BY) {
        final Function<FingerprintIdDetails, String> groupByTemp;
        switch (CROSSCHECK_BY) {
            case READGROUP:
                groupByTemp = details -> details.platformUnit;
                break;
            case LIBRARY:
                groupByTemp = details -> details.sample + "::" + details.library;
                break;
            case FILE:
                groupByTemp = details -> details.file + "::" + details.sample;
                break;
            case SAMPLE:
                groupByTemp = details -> details.sample;
                break;
            default:
                throw new PicardException("unpossible");
        }

        // if the groupBy string is null (e.g. a vcf file has no read group info) then the hashcode is
        // used intending to be unique per object (ignoring possible collisions)
        return key -> {
            final String temp = groupByTemp.apply(key);
            return temp == null ? Integer.toString(key.hashCode()) : temp;
        };
    }

    public static Map<FingerprintIdDetails, Fingerprint> mergeFingerprintsBy(
            final Map<FingerprintIdDetails, Fingerprint> fingerprints,
            final Function<FingerprintIdDetails, String> by){

        // collect the various entries according to the grouping "by"

        final Map<String, List<Map.Entry<FingerprintIdDetails, Fingerprint>>> collection =
                fingerprints.entrySet()
                        .stream()
                        .collect(Collectors.groupingBy(entry -> by.apply(entry.getKey())));

        return collection.entrySet().stream()
                .collect(Collectors.toMap(
                        entry -> {
                            // merge the keys (unequal values are eliminated by merge).

                            final FingerprintIdDetails finalId = new FingerprintIdDetails();
                            entry.getValue().forEach(id -> finalId.merge(id.getKey()));
                            finalId.group = entry.getKey();
                            return finalId;

                        }, entry -> {
                            // merge the values by merging the fingerprints.

                            final FingerprintIdDetails firstDetail = entry.getValue().get(0).getKey();
                            //use the "by" function to determine the "info" part of the fingerprint
                            final Fingerprint sampleFp = new Fingerprint(firstDetail.sample, null, by.apply(firstDetail));
                            entry.getValue().stream().map(Map.Entry::getValue).collect(Collectors.toSet()).forEach(sampleFp::merge);
                            return sampleFp;

                        }));
    }

    /**
     * Method that combines the fingerprint evidence across all the read groups for the same sample
     * and then produces a matrix of LOD scores for comparing every sample with every other sample.
     */
    private int crossCheckGrouped(final Map<FingerprintIdDetails, Fingerprint> fingerprints, final List<CrosscheckMetric> metrics,
                                  final Function<FingerprintIdDetails, String> by,
                                  final CrosscheckMetric.DataType type) {

        final Map<FingerprintIdDetails, Fingerprint> fingerprintsByGroup = mergeFingerprintsBy(fingerprints, by);

        if (MATRIX_OUTPUT != null) {
            crosscheckMatrix = new double[fingerprintsByGroup.size()][];
            for (int i = 0; i < fingerprintsByGroup.size(); i++) {
                crosscheckMatrix[i] = new double[fingerprintsByGroup.size()];
            }
        }
        return crossCheckFingerprints(fingerprintsByGroup, type, metrics, by);
    }

    /**
     * Method that pairwise checks every pair of groups and reports a LOD score for the two groups
     * coming from the same individual.
     */
    private int crossCheckFingerprints(final Map<FingerprintIdDetails, Fingerprint> fingerprints, final CrosscheckMetric.DataType type, final List<CrosscheckMetric> metrics, Function<FingerprintIdDetails, String> by) {
        int unexpectedResults = 0;
        long checksMade = 0;

        // TODO: create a general-use logger in htsjdk (not one that is tied to Chrom/Pos) and use that
        // TODO: to record progress

        final int logEvery = 100_000;

        final List<FingerprintIdDetails> fingerprintIdDetails = new ArrayList<>(fingerprints.keySet());

        // use 1L to promote size() to a long and avoid possible overflow
        final long totalChecks = fingerprintIdDetails.size() * (fingerprintIdDetails.size() + 1L) / 2;

        for (int i = 0; i < fingerprintIdDetails.size(); i++) {
            final FingerprintIdDetails lhsRg = fingerprintIdDetails.get(i);
            if (MATRIX_OUTPUT !=null) {
                matrixKeys.add(by.apply(lhsRg));
            }
            for (int j = i; j < fingerprintIdDetails.size(); j++) {
                final FingerprintIdDetails rhsRg = fingerprintIdDetails.get(j);
                final boolean expectedToMatch = EXPECT_ALL_GROUPS_TO_MATCH || lhsRg.sample.equals(rhsRg.sample);

                final MatchResults results = FingerprintChecker.calculateMatchResults(fingerprints.get(lhsRg), fingerprints.get(rhsRg), GENOTYPING_ERROR_RATE, LOSS_OF_HET_RATE);
                final FingerprintResult result = getMatchResults(expectedToMatch, results);

                if (!OUTPUT_ERRORS_ONLY || result == FingerprintResult.INCONCLUSIVE || !result.isExpected()) {
                    metrics.add(getMatchDetails(result, results, lhsRg, rhsRg, type));
                }
                if (result != FingerprintResult.INCONCLUSIVE && !result.isExpected()) unexpectedResults++;
                if (crosscheckMatrix != null) {
                    crosscheckMatrix[i][j] = results.getLOD();
                    crosscheckMatrix[j][i] = results.getLOD();
                }

                if (++checksMade % logEvery == 0) {
                    log.info("Compared " + checksMade + " of " + totalChecks);
                }
            }
        }
        return unexpectedResults;
    }

    /**
     * Method that checks each sample from fingerprints1 against that sample from fingerprints2 and reports a LOD score for the two groups
     * coming from the same individual.
     */
    private int checkFingerprintsBySample(final Map<FingerprintIdDetails, Fingerprint> fingerprints1, final Map<FingerprintIdDetails, Fingerprint> fingerprints2,
                                          final List<CrosscheckMetric> metrics) {
        int unexpectedResults = 0;

        final Map<FingerprintIdDetails, Fingerprint> fingerprints1BySample = mergeFingerprintsBy(fingerprints1, getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType.SAMPLE));
        final Map<FingerprintIdDetails, Fingerprint> fingerprints2BySample = mergeFingerprintsBy(fingerprints2, getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType.SAMPLE));

        final Map<String, FingerprintIdDetails> sampleToDetail1 = fingerprints1BySample.keySet().stream().collect(Collectors.toMap(id->id.group, id -> id));
        final Map<String, FingerprintIdDetails> sampleToDetail2 = fingerprints2BySample.keySet().stream().collect(Collectors.toMap(id->id.group, id -> id));

        Set<String> samples = new HashSet<>();
        samples.addAll(sampleToDetail1.keySet());
        samples.addAll(sampleToDetail2.keySet());

        for (final String sample : samples) {
            final FingerprintIdDetails lhsID = sampleToDetail1.get(sample);
            final FingerprintIdDetails rhsID = sampleToDetail2.get(sample);

            if (lhsID == null || rhsID == null) {
                log.error(String.format("sample %s is missing from %s group", sample, lhsID == null ? "LEFT" : "RIGHT"));
                unexpectedResults++;
                continue;
            }

            final MatchResults results = FingerprintChecker.calculateMatchResults(fingerprints1BySample.get(lhsID), fingerprints2BySample.get(rhsID), GENOTYPING_ERROR_RATE, LOSS_OF_HET_RATE);
            final CrosscheckMetric.FingerprintResult result = getMatchResults(true, results);

            if (!OUTPUT_ERRORS_ONLY || !result.isExpected()) {
                metrics.add(getMatchDetails(result, results, lhsID, rhsID, CrosscheckMetric.DataType.SAMPLE));
            }
            if (result != FingerprintResult.INCONCLUSIVE && !result.isExpected()) unexpectedResults++;
        }
        return unexpectedResults;
    }

    /**
     * Generates tab-delimited string containing details about a possible match between fingerprints on two different SAMReadGroupRecords
     *
     * @param matchResult    String describing the match type.
     * @param results        MatchResults object
     * @param leftPuDetails  left hand side FingerprintIdDetails
     * @param rightPuDetails right hand side FingerprintIdDetails
     * @return tab delimited string containing details about a possible match
     */
    private CrosscheckMetric getMatchDetails(final FingerprintResult matchResult,
                                             final MatchResults results,
                                             final FingerprintIdDetails leftPuDetails,
                                             final FingerprintIdDetails rightPuDetails,
                                             final CrosscheckMetric.DataType type) {
        final CrosscheckMetric metric = new CrosscheckMetric();

        metric.LEFT_GROUP_VALUE = leftPuDetails.group;
        metric.RIGHT_GROUP_VALUE = rightPuDetails.group;

        metric.RESULT = matchResult;
        metric.LOD_SCORE = results.getLOD();
        metric.LOD_SCORE_TUMOR_NORMAL = results.getLodTN();
        metric.LOD_SCORE_NORMAL_TUMOR = results.getLodNT();
        metric.DATA_TYPE = type;

        metric.LEFT_RUN_BARCODE = leftPuDetails.runBarcode;
        metric.LEFT_LANE = leftPuDetails.runLane;
        metric.LEFT_MOLECULAR_BARCODE_SEQUENCE = leftPuDetails.molecularBarcode;
        metric.LEFT_LIBRARY = leftPuDetails.library;
        metric.LEFT_SAMPLE = leftPuDetails.sample;
        metric.LEFT_FILE = leftPuDetails.file;

        metric.RIGHT_RUN_BARCODE = rightPuDetails.runBarcode;
        metric.RIGHT_LANE = rightPuDetails.runLane;
        metric.RIGHT_MOLECULAR_BARCODE_SEQUENCE = rightPuDetails.molecularBarcode;
        metric.RIGHT_LIBRARY = rightPuDetails.library;
        metric.RIGHT_SAMPLE = rightPuDetails.sample;
        metric.RIGHT_FILE = rightPuDetails.file;

        return metric;
    }

    private FingerprintResult getMatchResults(final boolean expectedToMatch, final MatchResults results) {
        if (expectedToMatch) {
            if (results.getLOD() < LOD_THRESHOLD) {
                return FingerprintResult.UNEXPECTED_MISMATCH;
            } else if (results.getLOD() > -LOD_THRESHOLD) {
                return FingerprintResult.EXPECTED_MATCH;
            } else {
                return FingerprintResult.INCONCLUSIVE;
            }
        } else {
            if (results.getLOD() > -LOD_THRESHOLD) {
                return FingerprintResult.UNEXPECTED_MATCH;
            } else if (results.getLOD() < LOD_THRESHOLD) {
                return FingerprintResult.EXPECTED_MISMATCH;
            } else {
                return FingerprintResult.INCONCLUSIVE;
            }
        }
    }
}
