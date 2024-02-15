package picard.sam.SamErrorMetric;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SamLocusIterator;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static picard.cmdline.CommandLineProgramTest.CHR_M_REFERENCE;

public class CollectSamErrorMetricsTest {
    private static File OUTPUT_DATA_PATH;
    private static final String TEST_DIR = "testdata/picard/sam/BamErrorMetrics";

    @BeforeClass
    public void setup() {
        ReadBaseStratification.setLongHomopolymer(6);
    }

    @BeforeTest()
    public void samMetricsProvider() {
        final File[] files = new File[]{
                simpleSamWithBaseErrors1,
                simpleSamWithBaseErrors2,
                simpleSamWithIndels,
                simpleSingleStrandConsensusSamWithBaseErrors,
                simpleDuplexConsensusSamWithBaseErrors,
                chrMReadsWithClips};

        if (OUTPUT_DATA_PATH != null) {
            throw new PicardException("This should only be called once!");
        }

        OUTPUT_DATA_PATH = IOUtil.createTempDir("CollectSamErrorMetricsTest.tmp").toFile();

        for (final File file : files) {
            final File vcf = new File(TEST_DIR, "NIST.selected.vcf");

            final File outputBaseFileName = new File(OUTPUT_DATA_PATH, file.getName());
            outputBaseFileName.deleteOnExit();

            final String[] args = {
                    "INPUT=" + file,
                    "OUTPUT=" + outputBaseFileName,
                    "REFERENCE_SEQUENCE=" + CHR_M_REFERENCE.getAbsolutePath(),
                    "ERROR_METRICS=" + "ERROR:TWO_BASE_PADDED_CONTEXT", // Not all covariates are included by default, but we still want to test them.
                    "ERROR_METRICS=" + "ERROR:CONSENSUS",
                    "ERROR_METRICS=" + "ERROR:NS_IN_READ",
                    "ERROR_METRICS=" + "ERROR:BINNED_CYCLE",
                    "VCF=" + vcf.getAbsolutePath()
            };

            Assert.assertEquals(new CollectSamErrorMetrics().instanceMain(args), 0);

            errorMetrics.put(file, outputBaseFileName);
        }
    }


    @DataProvider
    public Object[][] parseDirectiveData() {
        return new Object[][]{
                {"ERROR", "error_by_all"},
                {"ERROR:ALL", "error_by_all"},
                {"ERROR:GC_CONTENT", "error_by_gc"},
                {"ERROR:READ_ORDINALITY", "error_by_read_ordinality"},
                {"ERROR:READ_BASE", "error_by_read_base"},
                {"ERROR:REFERENCE_BASE", "error_by_ref_base"},
                {"ERROR:PRE_DINUC", "error_by_pre_dinuc"},
                {"ERROR:POST_DINUC", "error_by_post_dinuc"},
                {"ERROR:HOMOPOLYMER_LENGTH", "error_by_homopolymer_length"},
                {"ERROR:HOMOPOLYMER", "error_by_homopolymer_and_following_ref_base"},
                {"ERROR:BINNED_HOMOPOLYMER", "error_by_binned_length_homopolymer_and_following_ref_base"},
                {"ERROR:FLOWCELL_TILE", "error_by_tile"},
                {"ERROR:READ_DIRECTION", "error_by_read_direction"},
                {"ERROR:CYCLE", "error_by_cycle"},
                {"ERROR:BINNED_CYCLE", "error_by_binned_cycle"},
                {"ERROR:INSERT_LENGTH", "error_by_insert_length"},
                {"ERROR:BASE_QUALITY", "error_by_base_quality"},
                {"ERROR:MAPPING_QUALITY", "error_by_mapping_quality"},
                {"ERROR:READ_GROUP", "error_by_read_group"},
                {"ERROR:MISMATCHES_IN_READ", "error_by_mismatches_in_read"},
                {"ERROR:ONE_BASE_PADDED_CONTEXT", "error_by_one_base_padded_context"},
                {"ERROR:TWO_BASE_PADDED_CONTEXT", "error_by_two_base_padded_context"},
                {"ERROR:CONSENSUS", "error_by_consensus"},
                {"ERROR:NS_IN_READ", "error_by_ns_in_read"},

                {"ERROR:POST_DINUC:BASE_QUALITY", "error_by_post_dinuc_and_base_quality"},
                {"ERROR:POST_DINUC:BASE_QUALITY:GC_CONTENT", "error_by_post_dinuc_and_base_quality_and_gc"},
                {" ERROR : POST_DINUC : BASE_QUALITY : GC_CONTENT ", "error_by_post_dinuc_and_base_quality_and_gc"},

                {"OVERLAPPING_ERROR", "overlapping_error_by_all"},
                {"OVERLAPPING_ERROR:ALL", "overlapping_error_by_all"},
                {"OVERLAPPING_ERROR:GC_CONTENT", "overlapping_error_by_gc"},
                {"OVERLAPPING_ERROR:READ_ORDINALITY", "overlapping_error_by_read_ordinality"},
                {"OVERLAPPING_ERROR:READ_BASE", "overlapping_error_by_read_base"},
                {"OVERLAPPING_ERROR:REFERENCE_BASE", "overlapping_error_by_ref_base"},
                {"OVERLAPPING_ERROR:PRE_DINUC", "overlapping_error_by_pre_dinuc"},
                {"OVERLAPPING_ERROR:POST_DINUC", "overlapping_error_by_post_dinuc"},
                {"OVERLAPPING_ERROR:HOMOPOLYMER_LENGTH", "overlapping_error_by_homopolymer_length"},
                {"OVERLAPPING_ERROR:HOMOPOLYMER", "overlapping_error_by_homopolymer_and_following_ref_base"},
                {"OVERLAPPING_ERROR:FLOWCELL_TILE", "overlapping_error_by_tile"},
                {"OVERLAPPING_ERROR:READ_DIRECTION", "overlapping_error_by_read_direction"},
                {"OVERLAPPING_ERROR:CYCLE", "overlapping_error_by_cycle"},
                {"OVERLAPPING_ERROR:INSERT_LENGTH", "overlapping_error_by_insert_length"},
                {"OVERLAPPING_ERROR:BASE_QUALITY", "overlapping_error_by_base_quality"},
        };
    }

    @Test(dataProvider = "parseDirectiveData")
    public void parseDirectiveGood(final String directive, final String extension) {
        parseDirective0(directive, extension);
    }

    @Test
    void testStratifiersHaveDistinctSuffixes() {
        final Set<String> suffixes = new HashSet<>();

        for (final ReadBaseStratification.Stratifier stratifier : ReadBaseStratification.Stratifier.values()) {
            Assert.assertTrue(suffixes.add(stratifier.makeStratifier().getSuffix()), "found duplicate suffix: " +
                    stratifier.makeStratifier().getSuffix() + " for: " + stratifier.makeStratifier());
        }
    }

    @Test
    void testAggregatorsHaveDistinctSuffixes() {
        final Set<String> suffixes = new HashSet<>();

        for (final ErrorType error : ErrorType.values()) {
            final String suffix = error.getErrorSupplier().get().getSuffix();
            Assert.assertTrue(suffixes.add(suffix), "found duplicate suffix: " + suffix);
        }
    }

    @DataProvider
    public Object[][] parseDirectiveBadData() {
        return new Object[][]{
                {"ERROR:"},
                {"ERRORS:READ_ORDINALITY"},
                {"ERROR;REFERENCE_BASE"},
                {"ERROR:what"},
        };
    }

    @Test(dataProvider = "parseDirectiveBadData", expectedExceptions = IllegalArgumentException.class)
    public void parseDirectiveBad(final String directive) {
        parseDirective0(directive, null);
    }

    private static void parseDirective0(final String directive, final String extension) {
        final BaseErrorAggregation agg = CollectSamErrorMetrics.parseDirective(directive);
        if (extension != null) {
            Assert.assertEquals(agg.getSuffix(), extension);
        }
    }

    private static final Map<File, File> errorMetrics = new HashMap<>();

    private static final File simpleSamWithBaseErrors1 = new File(TEST_DIR, "simpleSamWithBaseErrors1.sam");
    private static final File simpleSamWithBaseErrors2 = new File(TEST_DIR, "simpleSamWithBaseErrors2.sam");
    private static final File simpleSamWithIndels = new File(TEST_DIR, "simpleSamWithIndels.sam");
    private static final File simpleSingleStrandConsensusSamWithBaseErrors = new File(TEST_DIR, "simpleSingleStrandConsensusSamWithBaseErrors.sam");
    private static final File simpleDuplexConsensusSamWithBaseErrors = new File(TEST_DIR, "simpleDuplexConsensusSamWithBaseErrors.sam");
    private static final File chrMReadsWithClips = new File(TEST_DIR, "chrMReadsWithClips.sam");


    @DataProvider
    public Object[][] OneCovariateErrorMetricsDataProvider() {

        final int priorQ = 30;
        //These magic numbers come from a separate implementation of the code in R.
        return new Object[][]{

                // Note that soft clipped bases are not counted.
                {".error_by_all", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("all", 62L, 49L)},
                {".error_by_base_quality", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("32", 52L, 41L)},
                {".error_by_base_quality", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("33", 10L, 8L)},
//                 Note that the homopolymer is counted as the number of bases in the read that match each other before a new base.
//                 This catches mismatches (and matches) of the ends of homopolymers.
                {".error_by_homopolymer_and_following_ref_base", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("9,T,T", 1L, 1L)},
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("2,G,T", 2L, 1L)},

                {".error_by_binned_length_homopolymer_and_following_ref_base", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("SHORT_HOMOPOLYMER,G,T", 6L, 1L)},

                {".error_by_read_ordinality_and_pre_dinuc", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("FIRST,A,A", 3L, 3L)},
                // Using a sam file with a single error it is easy to validate demonstrate these tests should pass
                {".error_by_all", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("all", 72L, 1L)},
                // There are two base qualities in the bam, the error occurs in quality "32"
                {".error_by_base_quality", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("32", 51L, 1L)},
                // There are two base qualities in the bam, the error occurs in quality "32", make sure we detect no errors in quality "33"
                {".error_by_base_quality", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("33", 21L, 0L)},
                // simpleSamWithIndels contains two insertions and one deletion.
                // One of the insertions should not be counted in total bases because it has base quality < 20.
                // There should be zero errors because insertions and deletions aren't included as errors.
                {".error_by_base_quality", simpleSamWithIndels, priorQ,
                        new BaseErrorMetric("32", 71L, 0L)},
                // Total bases should be 71 because low quality bases should not be included
                // Error bases should be zero because a low quality mismatch does not count as an error
                {".error_by_all", simpleSamWithIndels, priorQ,
                        new BaseErrorMetric("all", 71L, 0L)},
                // simpleSamWithBaseErrors2 contains 2 differences from the reference
                // after 2 different homopolymers.
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("2,T,A", 2L, 1L)},
                // Make sure that we can correctly identify an error after a homopolymer
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("3,C,G", 1L, 1L)},
                {".error_by_read_ordinality_and_pre_dinuc", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("FIRST,G,T", 1L, 1L)},
                // The covariate "0.5" was chosen to avoid those that have repeating decimals and could have alternative representations as
                // a string.  GC is calculated over the entire read including clipped bases, while errors are calculated only over unclipped bases.
                {".error_by_gc", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("0.5", 36L, 1L)},
                // Make sure that we can detect errors at a particular cycle
                {".error_by_cycle", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("10", 2L, 1L)},
                // There should be one error in the read with mapping quality 60.
                {".error_by_mapping_quality", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("60", 36L, 1L)},
                // There should be no errors in the read with mapping quality 0.
                {".error_by_mapping_quality", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("50", 36L, 0L)},
                // One base has an error in the read group 62A40.2
                {".error_by_read_group", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("62A40.2", 72L, 1L)},
                // No additional mismatches are found on the read with 1 mismatch.
                {".error_by_mismatches_in_read", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("1", 35L, 0L)},
                // No additional mismatches are found on the read with 1 mismatch. (Just another way to check)
                {".error_by_mismatches_in_read", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("0", 37L, 1L)},
                // There should be no errors in the CAG context because it matches reference
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("CAG", 1L, 0L)},
                // There should be one error in the GTC context
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("GTC", 1L, 1L)},
                // There should be one error in the CTT context
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("CCT", 3L, 0L)},
                // There should be no errors in the ACGGG context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("ACGGG", 1L, 0L)},
                // There should be one error in the GGTCT context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("GGTCT", 1L, 1L)},
                // There should be no errors in the CTTGA context (appears in sam file as TCAaG in second read)
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("CTTGA", 1L, 0L)},
                // There should be one error in the CCGTG context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("CCGTG", 1L, 1L)},
                // Reads that don't have consensus tags should be stratified as UNKNOWN
                {".error_by_consensus", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("UNKNOWN", 72L, 1L)},
                // There should be 2 errors in the one simplex singleton reads.
                {".error_by_consensus", simpleSingleStrandConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("SIMPLEX_SINGLETON", 36L, 2L)},
                // There should be no errors in the simplex consensus read.
                {".error_by_consensus", simpleSingleStrandConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("SIMPLEX_CONSENSUS", 36L, 0L)},
                // There should be one error in duplex singleton read.  Also the N in this read reduces total bases from 36 to 35.
                {".error_by_consensus", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("DUPLEX_SINGLETON", 35L, 1L)},
                // There should be two errors in the duplex consensus read.
                {".error_by_consensus", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("DUPLEX_CONSENSUS", 36L, 2L)},
                // There should be two errors in the read with no Ns.
                {".error_by_ns_in_read", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("0", 36L, 2L)},
                // There should be one errors in the read with one N.
                {".error_by_ns_in_read", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("1", 35L, 1L)},
                // There are two errors, one should show up in QUINTILE_1 and the other in QUINTILE_3
                // QUINTILE_5 has 16 total (2 more than the other bins) bases due to rounding.
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_1", 14L, 1L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_2", 14L, 0L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_3", 14L, 1L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_4", 14L, 0L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_5", 16L, 0L)}
        };
    }

    @Test(dataProvider = "OneCovariateErrorMetricsDataProvider")
    public void testOneCovariateErrorMetrics(final String errorSubscript, final File samFile, final int priorQ, BaseErrorMetric expectedMetric) {

        ErrorMetric.setPriorError(QualityUtil.getErrorProbabilityFromPhredScore(priorQ));
        expectedMetric.calculateDerivedFields();

        // Note that soft clipped bases are not counted
        List<BaseErrorMetric> metrics = MetricsFile.readBeans(new File(errorMetrics.get(samFile).getAbsolutePath() + errorSubscript));
        BaseErrorMetric metric = metrics
                .stream()
                .filter(m -> m.COVARIATE.equals(expectedMetric.COVARIATE))
                .findAny()
                .orElseThrow(() -> new AssertionError("didn't find metric with COVARIATE==" + expectedMetric.COVARIATE));
        Assert.assertEquals(metric, expectedMetric);
    }

    @DataProvider()
    public Object[][] oneCovariateIndelErrorMetricsDataProvider() {
        final int priorQ = 30;

        //These magic numbers come from a separate implementation of the code in R.
        return new Object[][]{

                // Note that soft clipped bases are not counted.
                {".error_by_all", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("all", 62L, 49L)},
                {".error_by_base_quality", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("32", 52L, 41L)},
                {".error_by_base_quality", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("33", 10L, 8L)},
//                 Note that the homopolymer is counted as the number of bases in the read that match each other before a new base.
//                 This catches mismatches (and matches) of the ends of homopolymers.
                {".error_by_homopolymer_and_following_ref_base", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("9,T,T", 1L, 1L)},
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("2,G,T", 2L, 1L)},

                {".error_by_binned_length_homopolymer_and_following_ref_base", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("SHORT_HOMOPOLYMER,G,T", 6L, 1L)},

                {".error_by_read_ordinality_and_pre_dinuc", chrMReadsWithClips, priorQ,
                        new BaseErrorMetric("FIRST,A,A", 3L, 3L)},
                // Using a sam file with a single error it is easy to validate demonstrate these tests should pass
                {".error_by_all", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("all", 72L, 1L)},
                // There are two base qualities in the bam, the error occurs in quality "32"
                {".error_by_base_quality", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("32", 51L, 1L)},
                // There are two base qualities in the bam, the error occurs in quality "32", make sure we detect no errors in quality "33"
                {".error_by_base_quality", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("33", 21L, 0L)},
                // simpleSamWithBaseErrors2 contains 2 differences from the reference
                // after 2 different homopolymers.
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("2,T,A", 2L, 1L)},
                // Make sure that we can correctly identify an error after a homopolymer
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("3,C,G", 1L, 1L)},
                {".error_by_read_ordinality_and_pre_dinuc", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("FIRST,G,T", 1L, 1L)},
                // The covariate "0.5" was chosen to avoid those that have repeating decimals and could have alternative representations as
                // a string.  GC is calculated over the entire read including clipped bases, while errors are calculated only over unclipped bases.
                {".error_by_gc", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("0.5", 36L, 1L)},
                // Make sure that we can detect errors at a particular cycle
                {".error_by_cycle", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("10", 2L, 1L)},
                // There should be one error in the read with mapping quality 60.
                {".error_by_mapping_quality", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("60", 36L, 1L)},
                // There should be no errors in the read with mapping quality 0.
                {".error_by_mapping_quality", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("50", 36L, 0L)},
                // One base has an error in the read group 62A40.2
                {".error_by_read_group", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("62A40.2", 72L, 1L)},
                // No additional mismatches are found on the read with 1 mismatch.
                {".error_by_mismatches_in_read", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("1", 35L, 0L)},
                // No additional mismatches are found on the read with 1 mismatch. (Just another way to check)
                {".error_by_mismatches_in_read", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("0", 37L, 1L)},
                // There should be no errors in the CAG context because it matches reference
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("CAG", 1L, 0L)},
                // There should be one error in the GTC context
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("GTC", 1L, 1L)},
                // There should be one error in the CTT context
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("CCT", 3L, 0L)},
                // There should be no errors in the ACGGG context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("ACGGG", 1L, 0L)},
                // There should be one error in the GGTCT context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("GGTCT", 1L, 1L)},
                // There should be no errors in the CTTGA context (appears in sam file as TCAaG in second read)
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("CTTGA", 1L, 0L)},
                // There should be one error in the CCGTG context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("CCGTG", 1L, 1L)},
                // Reads that don't have consensus tags should be stratified as UNKNOWN
                {".error_by_consensus", simpleSamWithBaseErrors1, priorQ,
                        new BaseErrorMetric("UNKNOWN", 72L, 1L)},
                // There should be 2 errors in the one simplex singleton reads.
                {".error_by_consensus", simpleSingleStrandConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("SIMPLEX_SINGLETON", 36L, 2L)},
                // There should be no errors in the simplex consensus read.
                {".error_by_consensus", simpleSingleStrandConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("SIMPLEX_CONSENSUS", 36L, 0L)},
                // There should be one error in duplex singleton read.  Also the N in this read reduces total bases from 36 to 35.
                {".error_by_consensus", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("DUPLEX_SINGLETON", 35L, 1L)},
                // There should be two errors in the duplex consensus read.
                {".error_by_consensus", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("DUPLEX_CONSENSUS", 36L, 2L)},
                // There should be two errors in the read with no Ns.
                {".error_by_ns_in_read", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("0", 36L, 2L)},
                // There should be one errors in the read with one N.
                {".error_by_ns_in_read", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new BaseErrorMetric("1", 35L, 1L)},
                // There are two errors, one should show up in QUINTILE_1 and the other in QUINTILE_3
                // QUINTILE_5 has 16 total (2 more than the other bins) bases due to rounding.
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_1", 14L, 1L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_2", 14L, 0L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_3", 14L, 1L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_4", 14L, 0L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new BaseErrorMetric("QUINTILE_5", 16L, 0L)}
        };
    }

    @Test(dataProvider = "oneCovariateIndelErrorMetricsDataProvider")
    public void testOneCovariateIndelErrorMetrics(final String errorSubscript, final File samFile, final int priorQ, BaseErrorMetric expectedMetric) {

        final File errorByAll = new File(errorMetrics.get(samFile).getAbsolutePath() + errorSubscript);
        errorByAll.deleteOnExit();

        ErrorMetric.setPriorError(QualityUtil.getErrorProbabilityFromPhredScore(priorQ));
        expectedMetric.calculateDerivedFields();

        // Note that soft clipped bases are not counted
        List<BaseErrorMetric> metrics = MetricsFile.readBeans(errorByAll);

        BaseErrorMetric metric = metrics
                .stream()
                .filter(m -> m.COVARIATE.equals(expectedMetric.COVARIATE))
                .findAny()
                .orElseThrow(() -> new AssertionError("didn't find metric with COVARIATE==" + expectedMetric.COVARIATE));

        Assert.assertEquals(metric, expectedMetric);
    }

    @AfterTest()
    public void cleanup() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider
    public Object[][] readCycleBinData() {
        // Test most edge cases of BaseErrorAggregation.CycleBin
        return new Object[][]{
                {0.0, ReadBaseStratification.CycleBin.QUINTILE_1},
                {0.2, ReadBaseStratification.CycleBin.QUINTILE_1},
                {0.2 + Math.ulp(0.2), ReadBaseStratification.CycleBin.QUINTILE_2},
                {0.4, ReadBaseStratification.CycleBin.QUINTILE_2},
                {0.4 + Math.ulp(0.4), ReadBaseStratification.CycleBin.QUINTILE_3},
                {0.54, ReadBaseStratification.CycleBin.QUINTILE_3},
                {0.6, ReadBaseStratification.CycleBin.QUINTILE_3},
                {0.6 + Math.ulp(0.6), ReadBaseStratification.CycleBin.QUINTILE_4},
                {0.8, ReadBaseStratification.CycleBin.QUINTILE_4},
                {0.8 + Math.ulp(0.8), ReadBaseStratification.CycleBin.QUINTILE_5},
                {1.0, ReadBaseStratification.CycleBin.QUINTILE_5},
        };
    }

    @Test(dataProvider = "readCycleBinData")
    public void testReadCycleBin(final double relativePosition, final ReadBaseStratification.CycleBin cycleBin) {
        Assert.assertEquals(ReadBaseStratification.CycleBin.valueOf(relativePosition), cycleBin);
    }

    @DataProvider
    public Object[][] readCycleBinDataError() {
        // Test cases that should throw an exception with BaseErrorAggregation.CycleBin
        return new Object[][]{
                {1.0 + Math.ulp(1.0)},
                {0.0 - Math.ulp(0.0)},
                {-1.0},
                {-0.1},
                {1.1},
                {100.0}
        };
    }

    @Test(dataProvider = "readCycleBinDataError", expectedExceptions = IllegalArgumentException.class)
    public void testReadCycleBinError(final double relativePosition) {
        ReadBaseStratification.CycleBin.valueOf(relativePosition);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTooManyDirectives() {
        final File input = new File(TEST_DIR, "simpleSamWithBaseErrors1.sam");

        final File referenceFile = CHR_M_REFERENCE;
        final File vcf = new File(TEST_DIR, "NIST.selected.vcf");

        final File outputBaseFileName = new File(OUTPUT_DATA_PATH, "test");
        outputBaseFileName.deleteOnExit();

        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("ERROR");
        for (ReadBaseStratification.Stratifier stratifier : ReadBaseStratification.Stratifier.values()) {
            stringBuilder.append(":");
            stringBuilder.append(stratifier.toString());
        }
        // one more should break it:
        stringBuilder.append(":ALL");
        final String[] args = {
                "INPUT=" + input,
                "OUTPUT=" + outputBaseFileName,
                "REFERENCE_SEQUENCE=" + referenceFile.getAbsolutePath(),
                "ERROR_METRICS=" + stringBuilder.toString(),
                "VCF=" + vcf.getAbsolutePath()
        };
        new CollectSamErrorMetrics().instanceMain(args);
    }

    public SAMRecord createRecordFromCigar(final String cigar) {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
        builder.addFrag("", 0, 100, false, false, cigar, null, 30);
        return builder.getRecords().stream().findFirst().get();
    }

    @DataProvider
    public Object[][] provideForTestHasDeletionBeenProcessed() {

        // Test most edge cases of BaseErrorAggregation.CycleBin
        return new Object[][]{
                {"1D", Arrays.asList(false)},
                {"2D", Arrays.asList(false, true)},
                {"1D1D", Arrays.asList(false, true)},
                {"1D1M1D", Arrays.asList(false, false)},
                {"2D1M1D", Arrays.asList(false, true, false)},
                {"3M2D8M1D", Arrays.asList(false, true, false)},
        };
    }

    @Test(dataProvider = "provideForTestHasDeletionBeenProcessed")
    public void testHasDeletionBeenProcessed(final String cigarString, final List<Boolean> expected) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final SAMRecord deletionRecord = createRecordFromCigar(cigar.toString());

        final CollectSamErrorMetrics collectSamErrorMetrics = new CollectSamErrorMetrics();
        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("chr1", 2000);

        int position = 100;
        int iExpected = 0;

        for(CigarElement cigarElement : cigar.getCigarElements()) {
            for(int iCigarElementPosition = 0; iCigarElementPosition < cigarElement.getLength(); iCigarElementPosition++) {
                final SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, ++position);
                if(cigarElement.getOperator() == CigarOperator.D) {
                    Assert.assertEquals(collectSamErrorMetrics.processDeletionLocus(new SamLocusIterator.RecordAndOffset(deletionRecord, 0), locusInfo), (boolean) expected.get(iExpected++));
                }
            }
        }
    }

    @DataProvider
    public Object[][] provideForTestIndelErrors() {
        return new Object[][] {
                // insertions
                { new String[] { "100M" },          ".indel_error_by_all", new IndelErrorMetric("all", 100, 0, 0, 0, 0) },
                { new String[] { "50M1I50M" },      ".indel_error_by_all", new IndelErrorMetric("all", 101, 1, 1, 0, 0) },
                { new String[] { "2I100M" },        ".indel_error_by_all", new IndelErrorMetric("all", 102, 1, 2, 0, 0) },
                { new String[] { "100M1I" },        ".indel_error_by_all", new IndelErrorMetric("all", 101, 1, 1, 0, 0) },
                { new String[] { "50M2I2M2I50M" },  ".indel_error_by_all", new IndelErrorMetric("all", 106, 2, 4, 0, 0) },
                { new String[] { "50M2I2M2I50M",
                                 "50M2I2M2I50M"},   ".indel_error_by_all", new IndelErrorMetric("all", 212, 4, 8, 0, 0) },
                // deletions
                { new String[] { "50M1D50M" },      ".indel_error_by_all", new IndelErrorMetric("all", 100, 0, 0, 1, 1) },
                { new String[] { "1D100M" },        ".indel_error_by_all", new IndelErrorMetric("all", 100, 0, 0, 1, 1) },
                { new String[] { "100M1D" },        ".indel_error_by_all", new IndelErrorMetric("all", 100, 0, 0, 1, 1) },
                { new String[] { "50M2D2M2D50M" },  ".indel_error_by_all", new IndelErrorMetric("all", 102, 0, 0, 2, 4) },
                { new String[] { "50M2D2M2D50M",
                                 "50M2D2M2D50M"},   ".indel_error_by_all", new IndelErrorMetric("all", 204, 0, 0, 4, 8) },
                // insertions & deletions
                { new String[] { "20M1I20M1D20M" }, ".indel_error_by_all", new IndelErrorMetric("all", 61, 1, 1, 1,1) },
                { new String[] { "20M2I2D20M" },    ".indel_error_by_all", new IndelErrorMetric("all", 42, 1, 2, 1,2) },
                { new String[] { "20M2D2I20M" },    ".indel_error_by_all", new IndelErrorMetric("all", 42, 1, 2, 1,2) },
                { new String[] { "2I2D20M" },       ".indel_error_by_all", new IndelErrorMetric("all", 22, 1, 2, 1,2) },
                { new String[] { "2D2I20M" },       ".indel_error_by_all", new IndelErrorMetric("all", 22, 1, 2, 1,2) },
                { new String[] { "20M2D2I" },       ".indel_error_by_all", new IndelErrorMetric("all", 22, 1, 2, 1,2) },
                { new String[] { "20M2I2D" },       ".indel_error_by_all", new IndelErrorMetric("all", 22, 1, 2, 1,2) },

                // indel_length:
                // insertions
                { new String[] { "50M1I50M" },      ".indel_error_by_indel_length", new IndelErrorMetric("1", 1, 1, 1, 0, 0) },
                { new String[] { "50M2I50M" },      ".indel_error_by_indel_length", new IndelErrorMetric("2", 2, 1, 2, 0, 0) },
                { new String[] { "20M2I20M2I20M" }, ".indel_error_by_indel_length", new IndelErrorMetric("2", 4, 2, 4, 0, 0) },
                { new String[] { "1I10M" },         ".indel_error_by_indel_length", new IndelErrorMetric("1", 1, 1, 1, 0, 0) },
                { new String[] { "10M1I" },         ".indel_error_by_indel_length", new IndelErrorMetric("1", 1, 1, 1, 0, 0) },
                // deletions
                { new String[] { "50M1D50M" },      ".indel_error_by_indel_length", new IndelErrorMetric("1", 0, 0, 0, 1, 1) },
                { new String[] { "50M2D50M" },      ".indel_error_by_indel_length", new IndelErrorMetric("2", 0, 0, 0, 1, 2) },
                { new String[] { "20M2D20M2D20M" }, ".indel_error_by_indel_length", new IndelErrorMetric("2", 0, 0, 0, 2, 4) },
                { new String[] { "1D10M" },         ".indel_error_by_indel_length", new IndelErrorMetric("1", 0, 0, 0, 1, 1) },
                { new String[] { "10M1D" },         ".indel_error_by_indel_length", new IndelErrorMetric("1", 0, 0, 0, 1, 1) },
                // insertions & deletions
                { new String[] { "20M2I20M3D20M" }, ".indel_error_by_indel_length", new IndelErrorMetric("2", 2, 1, 2, 0, 0) },
                { new String[] { "20M2I20M3D20M" }, ".indel_error_by_indel_length", new IndelErrorMetric("3", 0, 0, 0, 1, 3) },
                { new String[] { "2I2D20M" },       ".indel_error_by_indel_length", new IndelErrorMetric("2", 2, 1, 2, 1, 2) },
                { new String[] { "2D2I20M" },       ".indel_error_by_indel_length", new IndelErrorMetric("2", 2, 1, 2, 1, 2) },
                { new String[] { "2M2D2I" },        ".indel_error_by_indel_length", new IndelErrorMetric("2", 2, 1, 2, 1, 2) },
                { new String[] { "20M2I2D" },       ".indel_error_by_indel_length", new IndelErrorMetric("2", 2, 1, 2, 1, 2) },
        };
    }

    @Test(dataProvider = "provideForTestIndelErrors")
    public void testIndelErrors(final String[] readCigars, final String errorSubscript, IndelErrorMetric expectedMetric) throws IOException {

        final File temp = File.createTempFile("Indels", ".bam");
        temp.deleteOnExit();

        final int priorQ = 30;

        try (final ReferenceSequenceFileWalker referenceSequenceFileWalker =
                        new ReferenceSequenceFileWalker(CommandLineProgramTest.CHR_M_REFERENCE)) {

            final SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
            builder.getHeader().setSequenceDictionary(referenceSequenceFileWalker.getSequenceDictionary());

            for (int i = 0; i < readCigars.length; i++) {
                // 10M2I3M4D10M5I
                builder.addFrag("Read" + i, 0, 100, false, false, readCigars[i], null, 30);
            }

            try (final SAMFileWriter writer = new SAMFileWriterFactory()
                    .setCompressionLevel(2)
                    .makeBAMWriter(builder.getHeader(), false, temp)) {
                builder.forEach(writer::addAlignment);
            }
        }

        final File referenceFile = CHR_M_REFERENCE;
        final File vcf = new File(TEST_DIR, "NIST.selected.vcf");

        final File outputBaseFileName = new File(OUTPUT_DATA_PATH, "test");
        final File errorByAll = new File(outputBaseFileName.getAbsolutePath() + errorSubscript);
        errorByAll.deleteOnExit();
        outputBaseFileName.deleteOnExit();

        final String[] args = {
                "INPUT=" + temp,
                "OUTPUT=" + outputBaseFileName,
                "REFERENCE_SEQUENCE=" + referenceFile.getAbsolutePath(),
                "ERROR_METRICS=" + "INDEL_ERROR",
                "ERROR_METRICS=" + "INDEL_ERROR:INDEL_LENGTH",
                "VCF=" + vcf.getAbsolutePath()
        };

        CollectSamErrorMetrics collectSamErrorMetrics = new CollectSamErrorMetrics();
        collectSamErrorMetrics.ERROR_METRICS.clear();

        Assert.assertEquals(collectSamErrorMetrics.instanceMain(args), 0);

        ErrorMetric.setPriorError(QualityUtil.getErrorProbabilityFromPhredScore(priorQ));
        expectedMetric.calculateDerivedFields();

        // Note that soft clipped bases are not counted
        List<IndelErrorMetric> metrics = MetricsFile.readBeans(errorByAll);

        IndelErrorMetric metric = metrics
                .stream()
                .filter(m -> m.COVARIATE.equals(expectedMetric.COVARIATE))
                .findAny()
                .orElseThrow(() -> new AssertionError("didn't find metric with COVARIATE==" + expectedMetric.COVARIATE));

        Assert.assertEquals(metric, expectedMetric);
    }
}

