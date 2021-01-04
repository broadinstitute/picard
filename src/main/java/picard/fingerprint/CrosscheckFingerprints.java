
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

import htsjdk.samtools.SamReader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.fingerprint.CrosscheckMetric.FingerprintResult;
import picard.util.TabbedInputParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.file.Path;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import static picard.fingerprint.Fingerprint.CrosscheckMode.CHECK_SAME_SAMPLE;

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
 * In some cases, the groups collected may not have any observations (calls for a vcf, reads for a bam) at fingerprinting sites, or a sample in INPUT may be missing from the SECOND_INPUT.
 * These cases are handled as follows:  If running in CHECK_SAME_SAMPLES mode with INPUT and SECOND_INPUT, and either INPUT or SECOND_INPUT includes a sample
 * not found in the other, or contains a sample with no observations at any fingerprinting sites, an error will be logged and the tool will return EXIT_CODE_WHEN_MISMATCH.
 * In all other running modes, when any group which is being crosschecked does not have any observations at fingerprinting sites, a warning is logged.  As long as there is at least
 * one comparison where both sides have observations at fingerprinting sites, the tool will return zero.  However, if all comparisons have at least one side with no observations
 * at fingerprinting sites, an error will be logged and the tool will return EXIT_CODE_WHEN_NO_VALID_CHECKS.
 * <h3>Examples</h3>
 * <h4>Check that all the readgroups from a sample match each other:</h4>
 * <pre>
 *     java -jar picard.jar CrosscheckFingerprints \
 *          INPUT=sample.with.many.readgroups.bam \
 *          HAPLOTYPE_MAP=fingerprinting_haplotype_database.txt \
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
 *          HAPLOTYPE_MAP=fingerprinting_haplotype_database.txt \
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
 * @author Tim Fennell
 * @author Yossi Farjoun
 */

@CommandLineProgramProperties(
        summary =
                "Checks that all data in the set of input files appear to come from the same " +
                        "individual. Can be used to cross-check readgroups, libraries, samples, or files. " +
                        "Operates on SAM/BAM/CRAM and VCF (including gVCF and gzipped-VCF). " +
                        "\n" +
                        "<h3>Summary</h3>\n" +
                        "Checks if all the genetic data within a set of files appear to come from the same individual. " +
                        "It quickly determines whether a group's genotype matches that of an input file by selective sampling, " +
                        "and has been designed to work well for low-depth SAM (as well as high depth ones and VCFs.) " +
                        "The tool collects fingerprints (essentially, genotype information from different parts of the genome) " +
                        "at the finest level available in the data (readgroup for read-data files " +
                        "and sample for variant-data files) and then optionally aggregates it by library, sample or file, to increase power and provide " +
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
                        "to link a SAM to its correct sample in a multi-sample VCF. In this case one would not case for the cross-checking " +
                        "of the various samples in the VCF against each other, but only in checking the identity of the SAM against the various " +
                        "samples in the VCF. The SECOND_INPUT is provided for this use-case. With SECOND_INPUT provided, CrosscheckFingerprints " +
                        "does the following:\n" +
                        " - aggregation of data happens independently for the input files in INPUT and SECOND_INPUT. \n" +
                        " - aggregation of data happens at the SAMPLE level \n" +
                        " - each samples from INPUT will only be compared to that same sample in SECOND_INPUT. \n" +
                        " - MATRIX_OUTPUT is disabled. " +
                        "\n" +
                        "In some cases, the groups collected may not have any observations (calls for a VCF, reads for a SAM) at fingerprinting sites, or " +
                        "a sample in INPUT may be missing from the SECOND_INPUT. These cases are handled as follows:  If running in CHECK_SAME_SAMPLES mode " +
                        "with INPUT and SECOND_INPUT, and either INPUT or SECOND_INPUT includes a sample not found in the other, or contains a sample with " +
                        "no observations at any fingerprinting sites, an error will be logged and the tool will return EXIT_CODE_WHEN_MISMATCH. In all other " +
                        "running modes, when any group which is being crosschecked does not have any observations at fingerprinting sites, a warning is " +
                        "logged.  As long as there is at least one comparison where both sides have observations at fingerprinting sites, the tool will " +
                        "return zero.  However, if all comparisons have at least one side with no observations at fingerprinting sites, an error will be " +
                        "logged and the tool will return EXIT_CODE_WHEN_NO_VALID_CHECKS." +
                        "\n" +
                        "<hr/>" +
                        "<h3>Examples</h3>" +
                        "<h4>Check that all the readgroups from a sample match each other:</h4>" +
                        "<pre>" +
                        "    java -jar picard.jar CrosscheckFingerprints \\\n" +
                        "          INPUT=sample.with.many.readgroups.bam \\\n" +
                        "          HAPLOTYPE_MAP=fingerprinting_haplotype_database.txt \\\n" +
                        "          LOD_THRESHOLD=-5 \\\n" +
                        "          OUTPUT=sample.crosscheck_metrics" +
                        " </pre>" +
                        "\n" +
                        " <h4>Check that all the readgroups match as expected when providing reads from two samples from the same individual:</h4>" +
                        " <pre>" +
                        "     java -jar picard.jar CrosscheckFingerprints \\\n" +
                        "           INPUT=sample.one.with.many.readgroups.bam \\\n" +
                        "           INPUT=sample.two.with.many.readgroups.bam \\\n" +
                        "           HAPLOTYPE_MAP=fingerprinting_haplotype_database.txt \\\n" +
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
        programGroup = DiagnosticsAndQCProgramGroup.class

)
@DocumentedFeature
public class CrosscheckFingerprints extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input files (or lists of files) with which to compare fingerprints.", minElements = 1)
    public List<String> INPUT;

    @Argument(doc = "A tsv with two columns representing the sample as it appears in the INPUT data (in column 1) and " +
            "the sample as it should be used for comparisons to SECOND_INPUT (in the second column). " +
            "Need only include the samples that change. " +
            "Values in column 1 should be unique. " +
            "Values in column 2 should be unique even in union with the remaining unmapped samples. " +
            "Should only be used with SECOND_INPUT. ", optional = true, mutex = {"INPUT_SAMPLE_FILE_MAP"})
    public File INPUT_SAMPLE_MAP;

    @Argument(doc = "A tsv with two columns representing " +
            "the sample as it should be used for comparisons to SECOND_INPUT (in the first column) and  " +
            "the source file (in INPUT) for the fingerprint (in the second column). " +
            "Need only to include the samples that change. " +
            "Values in column 1 should be unique even in union with the remaining unmapped samples. " +
            "Values in column 2 should be unique in the file. " +
            "Will error if more than one sample is found in a file (multi-sample VCF) pointed to in column 2. " +
            "Should only be used in the presence of SECOND_INPUT. ", optional = true, mutex = {"INPUT_SAMPLE_MAP"})
    public File INPUT_SAMPLE_FILE_MAP;

    @Argument(shortName = "SI", optional = true, mutex = {"MATRIX_OUTPUT"},
            doc = "A second set of input files (or lists of files) with which to compare fingerprints. If this option is provided " +
                    "the tool compares each sample in INPUT with the sample from SECOND_INPUT that has the same sample ID. " +
                    "In addition, data will be grouped by SAMPLE regardless of the value of CROSSCHECK_BY. " +
                    "When operating in this mode, each sample in INPUT must also have a corresponding sample in SECOND_INPUT. " +
                    "If this is violated, the tool will proceed to check the matching samples, but report the missing samples " +
                    "and return a non-zero error-code.")
    public List<String> SECOND_INPUT;

    @Argument(doc = "A tsv with two columns representing the sample as it appears in the SECOND_INPUT data (in column 1) and " +
            "the sample as it should be used for comparisons to INPUT (in the second column). " +
            "Note that in case of unrolling files (file-of-filenames) one would need to reference the final file, i.e. the file that " +
            "contains the genomic data. Need only include the samples that change. " +
            "Values in column 1 should be unique. " +
            "Values in column 2 should be unique even in union with the remaining unmapped samples. " +
            "Should only be used with SECOND_INPUT. ", optional = true)
    public File SECOND_INPUT_SAMPLE_MAP;

    @Argument(doc = "An argument that controls how crosschecking with both INPUT and SECOND_INPUT should occur. ")
    public Fingerprint.CrosscheckMode CROSSCHECK_MODE = CHECK_SAME_SAMPLE;

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

    @Argument(doc = "specifies whether the Tumor-aware result should be calculated. These are time consuming and can roughly double the " +
            "runtime of the tool. When crosschecking many groups not calculating the tumor-aware  results can result in a significant speedup.")
    public boolean CALCULATE_TUMOR_AWARE_RESULTS = true;

    @Argument(doc = "Allow the use of duplicate reads in performing the comparison. Can be useful when duplicate " +
            "marking has been overly aggressive and coverage is low.")
    public boolean ALLOW_DUPLICATE_READS = false;

    @Argument(doc = "Assumed genotyping error rate that provides a floor on the probability that a genotype comes from " +
            "the expected sample. Must be greater than zero. ")
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

    @Argument(doc = "When all LOD score are zero, exit with this value.")
    public int EXIT_CODE_WHEN_NO_VALID_CHECKS = 1;

    @Argument(doc = "Maximal effect of any single haplotype block on outcome (-log10 of maximal likelihood difference between the different values for the three possible genotypes).", minValue = 0)
    public double MAX_EFFECT_OF_EACH_HAPLOTYPE_BLOCK = 3.0;

    @Hidden
    @Argument(doc = "When true code will check for readability on input files (this can be slow on cloud access)")
    public boolean TEST_INPUT_READABILITY = true;

    private final Log log = Log.getInstance(CrosscheckFingerprints.class);

    private double[][] crosscheckMatrix = null;
    private final List<String> lhsMatrixKeys = new ArrayList<>();
    private final List<String> rhsMatrixKeys = new ArrayList<>();

    @Override
    protected String[] customCommandLineValidation() {
        if (GENOTYPING_ERROR_RATE <= 0) {
            return new String[]{"GENOTYPING_ERROR_RATE must be greater than zero. Found " + GENOTYPING_ERROR_RATE};
        }
        if (GENOTYPING_ERROR_RATE >= 1) {
            return new String[]{"GENOTYPING_ERROR_RATE must be strictly less than 1, found " + GENOTYPING_ERROR_RATE};
        }
        if (SECOND_INPUT == null && INPUT_SAMPLE_MAP != null) {
            return new String[]{"INPUT_SAMPLE_MAP can only be used when also using SECOND_INPUT"};
        }
        if (SECOND_INPUT == null && SECOND_INPUT_SAMPLE_MAP != null) {
            return new String[]{"SECOND_INPUT_SAMPLE_MAP can only be used when also using SECOND_INPUT"};
        }


        //check that reference is provided if using crams as input
        if (REFERENCE_SEQUENCE == null) {
            final List<String> allInputs = new ArrayList<>(INPUT);
            allInputs.addAll(SECOND_INPUT);
            for (final String input : allInputs) {
                if (input.endsWith(SamReader.Type.CRAM_TYPE.fileExtension())) {
                    return new String[]{"REFERENCE must be provided when using CRAM as input."};
                }
            }
        }

        return super.customCommandLineValidation();
    }

    @Override
    protected int doWork() {
        // Check inputs

        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        if (OUTPUT != null) {
            IOUtil.assertFileIsWritable(OUTPUT);
        }

        if (!SECOND_INPUT.isEmpty() && CROSSCHECK_MODE == CHECK_SAME_SAMPLE) {
            log.info("SECOND_INPUT is not empty, and CROSSCHECK_MODE==CHECK_SAME_SAMPLE. NOT doing cross-check. Will only compare each SAMPLE in INPUT against that sample in SECOND_INPUT.");
            if (CROSSCHECK_BY != CrosscheckMetric.DataType.SAMPLE) {
                log.warn("CROSSCHECK_BY is not SAMPLE, This doesn't make sense in non-crosscheck mode. Setting CROSSCHECK_BY to SAMPLE.");
                CROSSCHECK_BY = CrosscheckMetric.DataType.SAMPLE;
            }
        }

        if (!SECOND_INPUT.isEmpty() && CROSSCHECK_MODE == Fingerprint.CrosscheckMode.CHECK_ALL_OTHERS) {
            log.info("SECOND_INPUT is not empty, and CROSSCHECK_MODE==CHECK_ALL_OTHERS. Will compare fingerprints from INPUT against all the fingerprints in SECOND_INPUT.");
        }

        if (MATRIX_OUTPUT != null) {
            IOUtil.assertFileIsWritable(MATRIX_OUTPUT);
        }
        if (INPUT_SAMPLE_MAP != null) {
            IOUtil.assertFileIsReadable(INPUT_SAMPLE_MAP);
        }
        if (INPUT_SAMPLE_FILE_MAP != null) {
            IOUtil.assertFileIsReadable(INPUT_SAMPLE_FILE_MAP);
        }
        if (SECOND_INPUT_SAMPLE_MAP != null) {
            IOUtil.assertFileIsReadable(SECOND_INPUT_SAMPLE_MAP);
        }

        final HaplotypeMap map = new HaplotypeMap(HAPLOTYPE_MAP);
        final FingerprintChecker checker = new FingerprintChecker(map);

        checker.setAllowDuplicateReads(ALLOW_DUPLICATE_READS);
        checker.setValidationStringency(VALIDATION_STRINGENCY);
        checker.setReferenceFasta(REFERENCE_SEQUENCE);

        final List<String> extensions = new ArrayList<>();

        extensions.add(SamReader.Type.BAM_TYPE.fileExtension());
        extensions.add(SamReader.Type.SAM_TYPE.fileExtension());
        extensions.add(SamReader.Type.CRAM_TYPE.fileExtension());
        extensions.addAll(FileExtensions.VCF_LIST);

        final List<Path> inputPaths = IOUtil.getPaths(INPUT);

        final List<Path> unrolledFiles = IOUtil.unrollPaths(inputPaths, extensions.toArray(new String[0]));
        if (TEST_INPUT_READABILITY) {
            IOUtil.assertPathsAreReadable(unrolledFiles);
        }

        final List<Path> secondInputsPaths = IOUtil.getPaths(SECOND_INPUT);

        // unroll and check readable here, as it can be annoying to fingerprint INPUT files and only then discover a problem
        // in a file in SECOND_INPUT
        final List<Path> unrolledFiles2 = IOUtil.unrollPaths(secondInputsPaths, extensions.toArray(new String[0]));
        if (TEST_INPUT_READABILITY) {
            IOUtil.assertPathsAreReadable(unrolledFiles2);
        }

        log.info("Fingerprinting " + unrolledFiles.size() + " INPUT files.");

        final Map<FingerprintIdDetails, Fingerprint> uncappedFpMap = checker.fingerprintFiles(unrolledFiles, NUM_THREADS, 1, TimeUnit.DAYS);
        final Map<FingerprintIdDetails, Fingerprint> fpMap = capFingerprints(uncappedFpMap);

        if (INPUT_SAMPLE_MAP != null) {
            remapFingerprints(fpMap, INPUT_SAMPLE_MAP, "INPUT_SAMPLE_MAP");
        }

        if (INPUT_SAMPLE_FILE_MAP != null) {
            remapFingerprintsFromFiles(fpMap, INPUT_SAMPLE_FILE_MAP);
        }


        final List<CrosscheckMetric> metrics = new ArrayList<>();
        final int numUnexpected;

        if (SECOND_INPUT.isEmpty()) {
            log.info("Cross-checking all " + CROSSCHECK_BY + " against each other");
            numUnexpected = crossCheckGrouped(fpMap, fpMap, metrics, Fingerprint.getFingerprintIdDetailsStringFunction(CROSSCHECK_BY), CROSSCHECK_BY);
        } else {
            log.info("Fingerprinting " + unrolledFiles2.size() + " SECOND_INPUT files.");
            final Map<FingerprintIdDetails, Fingerprint> uncappedFpMap2 = checker.fingerprintFiles(unrolledFiles2, NUM_THREADS, 1, TimeUnit.DAYS);
            final Map<FingerprintIdDetails, Fingerprint> fpMap2 = capFingerprints(uncappedFpMap2);


            if (SECOND_INPUT_SAMPLE_MAP != null) {
                remapFingerprints(fpMap2, SECOND_INPUT_SAMPLE_MAP, "SECOND_INPUT_SAMPLE_MAP");
            }
            switch (CROSSCHECK_MODE) {
                case CHECK_SAME_SAMPLE:
                    log.info("Checking each sample in INPUT with the same sample in SECOND_INPUT.");
                    numUnexpected = checkFingerprintsBySample(fpMap, fpMap2, metrics);
                    break;
                case CHECK_ALL_OTHERS:
                    log.info("Checking each " + CROSSCHECK_BY + " in INPUT with each " + CROSSCHECK_BY + " in SECOND_INPUT.");
                    numUnexpected = crossCheckGrouped(fpMap, fpMap2, metrics, Fingerprint.getFingerprintIdDetailsStringFunction(CROSSCHECK_BY), CROSSCHECK_BY);
                    break;
                default:
                    throw new IllegalArgumentException("Unpossible!");
            }
        }
        //check if all LODs are 0
        if (metrics.stream().noneMatch(m -> m.LOD_SCORE != 0)) {
            log.error("No non-zero results found. This is likely an error. " +
                    "Probable cause: there are no reads or variants at fingerprinting sites ");
            return EXIT_CODE_WHEN_NO_VALID_CHECKS;
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
            log.warn(numUnexpected + " " + CROSSCHECK_BY + "s did not relate as expected.");
            return EXIT_CODE_WHEN_MISMATCH;
        } else {
            log.info("All " + CROSSCHECK_BY + "s are related as expected.");
            return 0;
        }
    }

    private Map<FingerprintIdDetails, Fingerprint> capFingerprints(final Map<FingerprintIdDetails, Fingerprint> fpMap) {
        return fpMap.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, entry -> {
                    final Fingerprint value = entry.getValue();
                    final Fingerprint cappedFp = new Fingerprint(value.getSample(), value.getSource(), value.getInfo());

                    value.values().stream()
                            .map(probabilities -> new CappedHaplotypeProbabilities(probabilities, -MAX_EFFECT_OF_EACH_HAPLOTYPE_BLOCK))
                            .forEach(cappedFp::add);
                    return cappedFp;
                }));
    }

    /**
     * Inspects the contents of sampleMapFile building a map of Sample->Sample.
     * Checks for sanity, and then replaces in the fpMap,
     *
     * @param fpMap         A fingerprint map from a FingerprintIdDetails to a fingerprint
     * @param sampleMapFile a file contining two columns, sample name to be mapped, and mapped sample name.
     */
    private void remapFingerprints(final Map<FingerprintIdDetails, Fingerprint> fpMap, final File sampleMapFile, final String inputFieldName) {
        final Map<String, String> sampleMap = getStringStringMap(sampleMapFile, inputFieldName);

        // check that every key in the sample map is a sample in the fpMap, and warn otherwise
        final Set<String> samplesInFpMap = fpMap.keySet().stream().map(id -> id.sample).collect(Collectors.toSet());
        final Set<String> samplesNotInSampleMap = sampleMap.keySet().stream()
                .filter(((Predicate<String>) samplesInFpMap::contains).negate())
                .collect(Collectors.toSet());
        if (!samplesNotInSampleMap.isEmpty()) {
            log.warn("Some samples in first column in the " + inputFieldName + " were not present as samples in fingerprinted file: [" +
                    String.join(", ", samplesNotInSampleMap) + "].");
        }

        // verify that resulting sample-set is unique
        final List<String> resultingSamples = new ArrayList<>(samplesInFpMap);
        sampleMap.keySet().forEach(s -> {
            if (resultingSamples.remove(s)) {
                resultingSamples.add(sampleMap.get(s));
            }
        });

        if (CollectionUtil.makeSet(resultingSamples.toArray(new String[0])).size() != resultingSamples.size()) {
            final Set<String> duplicates = new HashSet<>();
            final Set<String> unique = new HashSet<>();
            resultingSamples.forEach(s -> {
                if (unique.add(s)) {
                    duplicates.add(s);
                }
            });
            throw new IllegalArgumentException("After applying the mapping found in the " + inputFieldName + " the resulting " +
                    "sample names must be unique when taken together with the remaining unmapped samples. " +
                    "Duplicates are: " + Arrays.toString(duplicates.toArray()));
        }

        // replace samples with their mapped values:
        final Set<FingerprintIdDetails> ids = new HashSet<>(fpMap.keySet());
        ids.forEach(id -> {
            // if sample isn't in sampleMap, leave it alone
            if (!sampleMap.containsKey(id.sample)) {
                return;
            }
            // one needs to replace the item, not simply modify it so that it is placed correctly in the map (since the key is changing)
            final Fingerprint fingerprint = fpMap.remove(id);
            //update the key
            id.sample = sampleMap.get(id.sample);
            //put the fingerprint back in with the updates key
            fpMap.put(id, fingerprint);
        });
    }

    /**
     * Inspects the contents of sampleMapFile building a map of Sample->Sample.
     * Checks for sanity, and then replaces in the fpMap,
     */
    private void remapFingerprintsFromFiles(final Map<FingerprintIdDetails, Fingerprint> fpMap, final File sampleMapFile) {
        final Map<String, String> sampleMap = getStringStringMap(sampleMapFile, "INPUT_SAMPLE_FILE_MAP").entrySet().stream()
                .collect(Collectors.toMap(e -> {
                    try {
                        return IOUtil.getPath(e.getValue()).toUri().toString();
                    } catch (IOException e1) {
                        throw new PicardException("Trouble reading file: " + e.getValue(), e1);
                    }
                }, Map.Entry::getKey));

        // check that every key in the sample map is a sample in the fpMap, and warn otherwise
        final Set<String> filesInFpMap = fpMap.keySet().stream().map(id -> id.file).collect(Collectors.toSet());
        final Set<String> sampleNotInFpMap = sampleMap.keySet().stream()
                .filter(((Predicate<String>) filesInFpMap::contains).negate())
                .collect(Collectors.toSet());
        if (!sampleNotInFpMap.isEmpty()) {
            log.warn("Some samples from the first column in " + "INPUT_SAMPLE_FILE_MAP" + " were not found: " +
                    Arrays.toString(sampleNotInFpMap.toArray()));
        }

        final Map<String, List<FingerprintIdDetails>> fileFpDetailSetMap = fpMap.keySet().stream().collect(Collectors.groupingBy(s -> s.file));

        final Map<String, String> fileSampleMap = new HashMap<>();
        // check that each file in the map points only to one sample
        fileFpDetailSetMap.forEach((key, fingerprintIdDetails) -> {
            final Set<String> samples = fingerprintIdDetails.stream().map(id -> id.sample).collect(Collectors.toSet());
            if (samples.size() > 1) {
                throw new IllegalArgumentException("fingerprinting file (" + key +
                        "in INPUT_SAMPLE_FILE_MAP contains multiple samples: " +
                        String.join("", samples));
            }
            fileSampleMap.put(key, fingerprintIdDetails.get(0).sample);
        });

        // verify that resulting file-set is associated with a unique sample:

        final List<String> resultingSamples = new ArrayList<>(filesInFpMap);
        fileSampleMap.forEach((f, id) -> {
            if (resultingSamples.remove(id)) {
                resultingSamples.add(sampleMap.get(f));
            }
        });

        if (CollectionUtil.makeSet(resultingSamples.toArray(new String[0])).size() != resultingSamples.size()) {
            final Set<String> duplicates = new HashSet<>();
            final Set<String> unique = new HashSet<>();
            resultingSamples.forEach(s -> {
                if (unique.add(s)) {
                    duplicates.add(s);
                }
            });
            throw new IllegalArgumentException("After applying the mapping found in the " + "INPUT_SAMPLE_FILE_MAP" + " the resulting " +
                    "sample names must be unique when taken together with the remaining unmapped samples. " +
                    "Duplicates are: " + Arrays.toString(duplicates.toArray()));
        }

        // replace samples with their mapped values:
        final Set<FingerprintIdDetails> ids = new HashSet<>(fpMap.keySet());
        ids.forEach(id -> {
            // if sample isn't in sampleMap, leave it alone
            if (!sampleMap.containsKey(id.file)) {
                return;
            }
            // one needs to replace the item, not simply modify it so that it is placed correctly in the map (since the key is changing)
            final Fingerprint fingerprint = fpMap.remove(id);
            //update the key
            id.sample = sampleMap.get(id.file);
            //put the fingerprint back in with the updates key
            fpMap.put(id, fingerprint);
        });
    }

    private Map<String, String> getStringStringMap(final File sampleMapFile, final String inputFieldName) {
        final Map<String, String> sampleMap = new HashMap<>();

        final TabbedInputParser parser = new TabbedInputParser(false, sampleMapFile);

        // build the map
        for (final String[] strings : parser) {
            if (strings.length != 2) {
                throw new IllegalArgumentException("Each line of the " + inputFieldName + " must have exactly two strings separated by a tab. " +
                        "Found: [" + Arrays.toString(strings) +
                        "] right before [" + parser.getCurrentLine() + "], in " + sampleMapFile.getAbsolutePath());
            }
            if (sampleMap.containsKey(strings[0])) {
                throw new IllegalArgumentException("Strings in first column of the " + inputFieldName + " must be unique. found [" + strings[0] +
                        "] twice. Right before [" + parser.getCurrentLine() + "] in " + sampleMapFile.getAbsolutePath());
            }
            sampleMap.put(strings[0], strings[1]);
        }
        return sampleMap;
    }

    private void writeMatrix() {

        final NumberFormat format = NumberFormat.getInstance();
        format.setMaximumFractionDigits(4);

        try (final OutputStream stream = new FileOutputStream(MATRIX_OUTPUT);
             final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stream))) {

            // write the type by which the roll-up happened in the top left corner of the matrix
            writer.write(CROSSCHECK_BY.name());

            // write the names of the keys as the first row
            for (String rhsMatrixKey : rhsMatrixKeys) {
                writer.write('\t' + rhsMatrixKey);
            }
            writer.newLine();

            for (int row = 0; row < lhsMatrixKeys.size(); row++) {
                // write the key in the first column
                writer.write(lhsMatrixKeys.get(row));
                // and then write all the values
                for (final double lod : crosscheckMatrix[row]) {
                    writer.write('\t' + format.format(lod));
                }
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Method that crosschecks fingerprints from one list of fingerprints against those in another
     * putting the results in a List of CrosscheckMetics.
     */
    private int crossCheckGrouped(final Map<FingerprintIdDetails, Fingerprint> lhsFingerprints,
                                  final Map<FingerprintIdDetails, Fingerprint> rhsFingerprints,
                                  final List<CrosscheckMetric> metrics,
                                  final Function<FingerprintIdDetails, String> by,
                                  final CrosscheckMetric.DataType type) {

        final Map<FingerprintIdDetails, Fingerprint> lhsFingerprintsByGroup = Fingerprint.mergeFingerprintsBy(lhsFingerprints, by);
        final Map<FingerprintIdDetails, Fingerprint> rhsFingerprintsByGroup = Fingerprint.mergeFingerprintsBy(rhsFingerprints, by);

        for (final Map.Entry<FingerprintIdDetails, Fingerprint> pair : lhsFingerprintsByGroup.entrySet()) {
            if (pair.getValue().size() == 0) {
                log.warn(by.apply(pair.getKey()) + " was not fingerprinted in LEFT group.  It probably has no calls/reads overlapping fingerprinting sites.");
            }
        }

        for (final Map.Entry<FingerprintIdDetails, Fingerprint> pair : rhsFingerprintsByGroup.entrySet()) {
            if (pair.getValue().size() == 0) {
                log.warn(by.apply(pair.getKey()) + " was not fingerprinted in RIGHT group.  It probably has no calls/reads overlapping fingerprinting sites.");
            }
        }

        if (MATRIX_OUTPUT != null) {
            crosscheckMatrix = new double[lhsFingerprintsByGroup.size()][];
            for (int row = 0; row < lhsFingerprintsByGroup.size(); row++) {
                crosscheckMatrix[row] = new double[rhsFingerprintsByGroup.size()];
            }
            lhsFingerprintsByGroup.keySet().forEach(k -> lhsMatrixKeys.add(k.group));
            rhsFingerprintsByGroup.keySet().forEach(k -> rhsMatrixKeys.add(k.group));
        }
        return crossCheckFingerprints(lhsFingerprintsByGroup, rhsFingerprintsByGroup, type, metrics);
    }

    /**
     * Method that pairwise checks every pair of groups and reports a LOD score for the two groups
     * coming from the same individual.
     */
    private int crossCheckFingerprints(final Map<FingerprintIdDetails, Fingerprint> lhsFingerprints, final Map<FingerprintIdDetails, Fingerprint> rhsFingerprints, final CrosscheckMetric.DataType type, final List<CrosscheckMetric> metrics) {
        int unexpectedResults = 0;
        long checksMade = 0;

        final int logEvery = 100_000;

        final List<FingerprintIdDetails> lhsFingerprintIdDetails = new ArrayList<>(lhsFingerprints.keySet());
        final List<FingerprintIdDetails> rhsFingerprintIdDetails = new ArrayList<>(rhsFingerprints.keySet());

        // use 1L to promote size() to a long and avoid possible overflow
        final long totalChecks = lhsFingerprintIdDetails.size() * ((long) rhsFingerprintIdDetails.size());

        for (int row = 0; row < lhsFingerprintIdDetails.size(); row++) {
            final FingerprintIdDetails lhsId = lhsFingerprintIdDetails.get(row);

            for (int col = 0; col < rhsFingerprintIdDetails.size(); col++) {
                final FingerprintIdDetails rhsId = rhsFingerprintIdDetails.get(col);
                final boolean expectedToMatch = EXPECT_ALL_GROUPS_TO_MATCH || lhsId.sample.equals(rhsId.sample);

                final MatchResults results = FingerprintChecker.calculateMatchResults(lhsFingerprints.get(lhsId), rhsFingerprints.get(rhsId),
                        GENOTYPING_ERROR_RATE, LOSS_OF_HET_RATE, false, CALCULATE_TUMOR_AWARE_RESULTS);
                final FingerprintResult result = getMatchResults(expectedToMatch, results);

                if (!OUTPUT_ERRORS_ONLY || result == FingerprintResult.INCONCLUSIVE || !result.isExpected()) {
                    metrics.add(getMatchDetails(result, results, lhsId, rhsId, type));
                }
                if (result != FingerprintResult.INCONCLUSIVE && !result.isExpected()) {
                    unexpectedResults++;
                }
                if (crosscheckMatrix != null) {
                    crosscheckMatrix[row][col] = results.getLOD();
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

        final Map<FingerprintIdDetails, Fingerprint> fingerprints1BySample = Fingerprint.mergeFingerprintsBy(fingerprints1, Fingerprint.getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType.SAMPLE));
        final Map<FingerprintIdDetails, Fingerprint> fingerprints2BySample = Fingerprint.mergeFingerprintsBy(fingerprints2, Fingerprint.getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType.SAMPLE));

        final Map<String, FingerprintIdDetails> sampleToDetail1 = fingerprints1BySample.keySet().stream().collect(Collectors.toMap(id -> id.group, id -> id));
        final Map<String, FingerprintIdDetails> sampleToDetail2 = fingerprints2BySample.keySet().stream().collect(Collectors.toMap(id -> id.group, id -> id));

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

            final Fingerprint lhsFP = fingerprints1BySample.get(lhsID);
            final Fingerprint rhsFP = fingerprints2BySample.get(rhsID);

            if (lhsFP.size() == 0 || rhsFP.size() == 0) {
                log.error(String.format("sample %s from %s group was not fingerprinted.  Probably there are no reads/variants at fingerprinting sites.", sample, lhsFP.size() == 0 ? "LEFT" : "RIGHT"));
                unexpectedResults++;
            }

            final MatchResults results = FingerprintChecker.calculateMatchResults(lhsFP, rhsFP,
                    GENOTYPING_ERROR_RATE, LOSS_OF_HET_RATE, false, CALCULATE_TUMOR_AWARE_RESULTS);
            final CrosscheckMetric.FingerprintResult result = getMatchResults(true, results);

            if (!OUTPUT_ERRORS_ONLY || !result.isExpected()) {
                metrics.add(getMatchDetails(result, results, lhsID, rhsID, CrosscheckMetric.DataType.SAMPLE));
            }
            if (result != FingerprintResult.INCONCLUSIVE && !result.isExpected()) {
                unexpectedResults++;
            }
            if (results.getLOD() == 0) {
                log.error("LOD score of zero found when checking sample fingerprints.  Probably there are no reads/variants at fingerprinting sites for one of the samples");
                unexpectedResults++;
            }
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
