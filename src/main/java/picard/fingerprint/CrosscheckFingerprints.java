
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
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Fingerprinting;

import java.io.*;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.toSet;
import static picard.fingerprint.CrosscheckMetric.FingerprintResult;
import static picard.fingerprint.CrosscheckMetric.FingerprintResult.*;

/**
 * Program to check that all fingerprints within the set of input files appear to come from the same
 * individual. Can be used to cross-check readgroups, libraries, samples, or files.
 * Operates on bams/sams and vcfs.
 *
 * @author Tim Fennell
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = "Checks if all fingerprints within a set of files appear to come from the same individual. " +
                "The fingerprints are calculated initially at the readgroup level (if present) but can be " +
                "\"rolled-up\" by library, sample or file, to increase power and provide results at the " +
                "desired resolution. Regular output is in a \"Moltenized\" format, one row per comparison. " +
                "In this format the output will include the LOD score and also tumor-aware LOD score which can " +
                "help assess identity even in the presence of a severe LOH sample with high purity. " +
                "A matrix output is also available to facilitate visual inspection of crosscheck results." +
                "\n" +
                "A separate CLP, ClusterCrosscheckMetrics, can cluster the results as a connected graph " +
                "according to LOD greater than a threshold. ",
        oneLineSummary = "Checks if all fingerprints appear to come from the same individual.",
        programGroup = Fingerprinting.class
)
public class CrosscheckFingerprints extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "One or more input files (or lists of files) to compare fingerprints for.")
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true,
            doc = "Optional output file to write metrics to. Default is to write to stdout.")
    public File OUTPUT = null;

    @Argument(shortName = "MO", optional = true,
            doc = "Optional output file to write matrix of LOD scores to. This is less informative than the metrics output " +
                    "and only contains Normal-Normal LOD score (i.e. doesn't account for Loss of heterogeneity). " +
                    "It is however sometimes easier to use visually.")
    public File MATRIX_OUTPUT = null;

    @Argument(shortName = "H", doc = "The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
            "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
    public File HAPLOTYPE_MAP;

    @Argument(shortName = "LOD",
            doc = "If any two groups (with the same sample name) match with a LOD score lower than the threshold " +
                    "the program will exit with a non-zero code to indicate error." +
                    " Program will also exit with an error if it finds two groups with different sample name that " +
                    "match with a LOD score greater than -LOD_THRESHOLD." +
                    "\n\n" +
                    "LOD score 0 means equal likelihood" +
                    "that the groups match vs. come from different individuals, negative LOD scores mean N logs more likely " +
                    "that the groups are from different individuals, and positive numbers mean N logs more likely that " +
                    "the groups are from the sample individual. ")
    public double LOD_THRESHOLD = 0;

    @Argument(doc = "Specificies which data-type should be used as the basic comparison unit. Fingerprints from readgroups can " +
            "be \"rolled-up\" to the LIBRARY, SAMPLE, or FILE level before being compared." +
            " Fingerprints from VCF can be be compared by SAMPLE or FILE.")
    public CrosscheckMetric.DataType CROSSCHECK_BY = CrosscheckMetric.DataType.READGROUP;

    @Argument(doc = "The number of threads to use to process files and generate Fingerprints.")
    public int NUM_THREADS = 1;

    @Argument(doc = "Allow the use of duplicate reads in performing the comparison. Can be useful when duplicate " +
            "marking has been overly aggressive and coverage is low.")
    public boolean ALLOW_DUPLICATE_READS = false;

    @Argument(doc = "Assumed genotyping error rate that provides a floor on the probability that a genotype comes from " +
            "the expected sample.")
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
    protected int doWork() {
        // Check inputs
        INPUT.forEach(IOUtil::assertFileIsReadable);
        IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);
        if (OUTPUT != null) IOUtil.assertFileIsWritable(OUTPUT);
        if (MATRIX_OUTPUT != null) IOUtil.assertFileIsWritable(MATRIX_OUTPUT);

        final HaplotypeMap map = new HaplotypeMap(HAPLOTYPE_MAP);
        final FingerprintChecker checker = new FingerprintChecker(map);

        checker.setAllowDuplicateReads(ALLOW_DUPLICATE_READS);
        checker.setValidationStringency(VALIDATION_STRINGENCY);

        log.info("Done checking input files, moving onto fingerprinting files.");

        final List<String> extensions = new ArrayList<>();

        extensions.add(BamFileIoUtils.BAM_FILE_EXTENSION);
        extensions.add(IOUtil.SAM_FILE_EXTENSION);
        extensions.addAll(Arrays.asList(IOUtil.VCF_EXTENSIONS));

        final List<File> unrolledFiles = IOUtil.unrollFiles(INPUT, extensions.toArray(new String[extensions.size()]));

        final Map<FingerprintIdDetails, Fingerprint> fpMap = checker.fingerprintFiles(unrolledFiles, NUM_THREADS, 1, TimeUnit.DAYS);

        log.info("Finished generating fingerprints from files, moving on to cross-checking.");

        final List<CrosscheckMetric> metrics = new ArrayList<>();
        final int numUnexpected;

        numUnexpected = crossCheckGrouped(fpMap, metrics, getFingerprintIdDetailsStringFunction(CROSSCHECK_BY), CROSSCHECK_BY);

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
                groupByTemp = details -> details.file;
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

    static public Map<FingerprintIdDetails, Fingerprint> mergeFingerprintsBy(
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
                            entry.getValue().stream().map(Map.Entry::getValue).collect(toSet()).forEach(sampleFp::merge);
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
     * Method that pairwise checks every pair of read groups and reports a LOD score for the two read groups
     * coming from the same sample.
     */
    private int crossCheckFingerprints(final Map<FingerprintIdDetails, Fingerprint> fingerprints, final CrosscheckMetric.DataType type, final List<CrosscheckMetric> metrics, Function<FingerprintIdDetails, String> by) {
        int unexpectedResults = 0;

        final List<FingerprintIdDetails> fingerprintIdDetails = new ArrayList<>(fingerprints.keySet());

        for (int i = 0; i < fingerprintIdDetails.size(); i++) {
            final FingerprintIdDetails lhsRg = fingerprintIdDetails.get(i);
            matrixKeys.add(by.apply(lhsRg));

            for (int j = i; j < fingerprintIdDetails.size(); j++) {
                final FingerprintIdDetails rhsRg = fingerprintIdDetails.get(j);
                final boolean expectedToMatch = EXPECT_ALL_GROUPS_TO_MATCH || lhsRg.sample.equals(rhsRg.sample);

                final MatchResults results = FingerprintChecker.calculateMatchResults(fingerprints.get(lhsRg), fingerprints.get(rhsRg), GENOTYPING_ERROR_RATE, LOSS_OF_HET_RATE);

                final CrosscheckMetric.FingerprintResult result = getMatchResults(expectedToMatch, results);

                if (!OUTPUT_ERRORS_ONLY || !result.isExpected()) {
                    metrics.add(getMatchDetails(result, results, lhsRg, rhsRg, type));
                }
                if (result != INCONCLUSIVE && !result.isExpected()) unexpectedResults++;
                if (crosscheckMatrix != null) {
                    crosscheckMatrix[i][j] = results.getLOD();
                    crosscheckMatrix[j][i] = results.getLOD();
                }
            }
        }

        return unexpectedResults;
    }

    /**
     * Generates tab delimited string containing details about a possible match between fingerprints on two different SAMReadGroupRecords
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
                return UNEXPECTED_MISMATCH;
            } else if (results.getLOD() > -LOD_THRESHOLD) {
                return EXPECTED_MATCH;
            } else {
                return INCONCLUSIVE;
            }
        } else {
            if (results.getLOD() > -LOD_THRESHOLD) {
                return UNEXPECTED_MATCH;
            } else if (results.getLOD() < LOD_THRESHOLD) {
                return EXPECTED_MISMATCH;
            } else {
                return INCONCLUSIVE;
            }
        }
    }
}
