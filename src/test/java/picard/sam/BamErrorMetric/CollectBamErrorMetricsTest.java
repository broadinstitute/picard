package picard.sam.BamErrorMetric;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.QualityUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class CollectBamErrorMetricsTest {
    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("CollectBamErrorMetricsTest", null);
    private static final String TEST_DIR = "testdata/picard/sam/BamErrorMetrics";


    @DataProvider
    public Object[][] parseDirectiveData() {
        return new Object[][]{
                {"ERROR", "error_by_all"},
                {"ERROR,ALL", "error_by_all"},
                {"ERROR,GC_CONTENT", "error_by_gc"},
                {"ERROR,READ_ORDINALITY", "error_by_read_ordinality"},
                {"ERROR,CONFUSION_MATRIX", "error_by_ref_read_base"},
                {"ERROR,READ_BASE", "error_by_read_base"},
                {"ERROR,REFERENCE_BASE", "error_by_ref_base"},
                {"ERROR,PRE_DINUC", "error_by_pre_dinuc"},
                {"ERROR,POST_DINUC", "error_by_post_dinuc"},
                {"ERROR,HOMOPOLYMER_LENGTH", "error_by_homopolymer_length"},
                {"ERROR,HOMOPOLYMER", "error_by_homopolymer_and_following_ref_base"},
                {"ERROR,BINNED_HOMOPOLYMER", "error_by_binned_length_homopolymer_and_following_ref_base"},
                {"ERROR,FLOWCELL_TILE", "error_by_tile"},
                {"ERROR,READ_DIRECTION", "error_by_read_direction"},
                {"ERROR,CYCLE", "error_by_cycle"},
                {"ERROR,BINNED_CYCLE", "error_by_binned_cycle"},
                {"ERROR,INSERT_LENGTH", "error_by_insert_length"},
                {"ERROR,BASE_QUALITY", "error_by_base_quality"},
                {"ERROR,MAPPING_QUALITY", "error_by_mapping_quality"},
                {"ERROR,READ_GROUP", "error_by_read_group"},
                {"ERROR,MISMATCHES_IN_READ", "error_by_mismatches_in_read"},
                {"ERROR,ONE_BASE_PADDED_CONTEXT", "error_by_one_base_padded_context"},
                {"ERROR,TWO_BASE_PADDED_CONTEXT", "error_by_two_base_padded_context"},
                {"ERROR,CONSENSUS", "error_by_consensus"},
                {"ERROR,NS_IN_READ", "error_by_ns_in_read"},

                {"ERROR,POST_DINUC,BASE_QUALITY", "error_by_post_dinuc_and_base_quality"},
                {"ERROR,POST_DINUC,BASE_QUALITY,GC_CONTENT", "error_by_post_dinuc_and_base_quality_and_gc"},
                {" ERROR , POST_DINUC , BASE_QUALITY , GC_CONTENT ", "error_by_post_dinuc_and_base_quality_and_gc"},

                {"OVERLAPPING_ERROR", "overlapping_bases_error_by_all"},
                {"OVERLAPPING_ERROR,ALL", "overlapping_bases_error_by_all"},
                {"OVERLAPPING_ERROR,GC_CONTENT", "overlapping_bases_error_by_gc"},
                {"OVERLAPPING_ERROR,READ_ORDINALITY", "overlapping_bases_error_by_read_ordinality"},
                {"OVERLAPPING_ERROR,CONFUSION_MATRIX", "overlapping_bases_error_by_ref_read_base"},
                {"OVERLAPPING_ERROR,READ_BASE", "overlapping_bases_error_by_read_base"},
                {"OVERLAPPING_ERROR,REFERENCE_BASE", "overlapping_bases_error_by_ref_base"},
                {"OVERLAPPING_ERROR,PRE_DINUC", "overlapping_bases_error_by_pre_dinuc"},
                {"OVERLAPPING_ERROR,POST_DINUC", "overlapping_bases_error_by_post_dinuc"},
                {"OVERLAPPING_ERROR,HOMOPOLYMER_LENGTH", "overlapping_bases_error_by_homopolymer_length"},
                {"OVERLAPPING_ERROR,HOMOPOLYMER", "overlapping_bases_error_by_homopolymer_and_following_ref_base"},
                {"OVERLAPPING_ERROR,FLOWCELL_TILE", "overlapping_bases_error_by_tile"},
                {"OVERLAPPING_ERROR,READ_DIRECTION", "overlapping_bases_error_by_read_direction"},
                {"OVERLAPPING_ERROR,CYCLE", "overlapping_bases_error_by_cycle"},
                {"OVERLAPPING_ERROR,INSERT_LENGTH", "overlapping_bases_error_by_insert_length"},
                {"OVERLAPPING_ERROR,BASE_QUALITY", "overlapping_bases_error_by_base_quality"},
        };
    }

    @Test(dataProvider = "parseDirectiveData")
    public void parseDirectiveGood(final String directive, final String extension) throws NoSuchMethodException, InstantiationException, IllegalAccessException, InvocationTargetException {
        parseDirective0(directive, extension);

    }

    @Test
    void testStratifiersHaveDistinctSuffixes() {
        Set<String> suffixes = new HashSet<>();

        for (final ReadBaseStratification.Stratifiers stratifier : ReadBaseStratification.Stratifiers.values()) {
            Assert.assertTrue(suffixes.add(stratifier.makeStratifier().getSuffix()), "found duplicate suffix: " +
                    stratifier.makeStratifier().getSuffix() + " for: " + stratifier.makeStratifier());
        }
    }

    @Test
    void testAggregatorsHaveDistinctSuffixes() throws InvocationTargetException, NoSuchMethodException, InstantiationException, IllegalAccessException {
        Set<String> suffixes = new HashSet<>();

        final ReadBaseStratification.Stratifiers stratifier = ReadBaseStratification.Stratifiers.ALL;
        for (final BaseErrorCalculation.Errors error : BaseErrorCalculation.Errors.values()) {
            final String suffix = error.getErrorSupplier().get().getSuffix();
            Assert.assertTrue(suffixes.add(suffix), "found duplicate suffix: " + suffix);
        }
    }

    @DataProvider
    public Object[][] parseDirectiveBadData() {
        return new Object[][]{
                {"ERROR,", ""},
                {"ERRORS,READ_ORDINALITY", ""},
                {"ERROR,CONFUSION_MATRIX,", ""},
                {"ERROR;REFERENCE_BASE", ""},
                {"ERROR,what", ""},
        };
    }

    @Test(dataProvider = "parseDirectiveBadData", expectedExceptions = IllegalArgumentException.class)
    public void parseDirectiveBad(final String directive, final String extension) throws NoSuchMethodException, InstantiationException, IllegalAccessException, InvocationTargetException {
        parseDirective0(directive, extension);
    }

    public void parseDirective0(final String directive, final String extension) throws NoSuchMethodException, InstantiationException, IllegalAccessException, InvocationTargetException {
        final BaseErrorAggregation agg = CollectBamErrorMetrics.parseDirective(directive);
        Assert.assertEquals(agg.getSuffix(), extension);
    }

    @DataProvider(name = "OneCovariateErrorMetricsDataProvider")
    public Object[][] oneCovariateErrorMetricsDataProvider() {
        final File simpleSamWithBaseErrors1 = new File(TEST_DIR, "simpleSamWithBaseErrors1.sam");
        final File simpleSamWithBaseErrors2 = new File(TEST_DIR, "simpleSamWithBaseErrors2.sam");
        final File simpleSingleStrandConsensusSamWithBaseErrors = new File(TEST_DIR, "simpleSingleStrandConsensusSamWithBaseErrors.sam");
        final File simpleDuplexConsensusSamWithBaseErrors = new File(TEST_DIR, "simpleDuplexConsensusSamWithBaseErrors.sam");
        final File chrMReadsWithClips = new File(TEST_DIR, "chrMReadsWithClips.sam");
        final int priorQ = 30;

        //These magic numbers come from a separate implementation of the code in R.
        return new Object[][] {

                // Note that soft clipped bases are not counted.
                {".error_by_all", chrMReadsWithClips, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("all", 62L, 49L)},
                {".error_by_base_quality", chrMReadsWithClips, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("32", 52L, 41L)},
                {".error_by_base_quality", chrMReadsWithClips, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("33", 10L, 8L)},
                // Note that the homopolymer is counted as the number of bases in the read that match each other before a new base.
                // This catches mismatches (and matches) of the ends of homopolymers.
                {".error_by_homopolymer_and_following_ref_base", chrMReadsWithClips, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("T,T,10", 1L, 1L)},
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("G,T,3", 2L, 1L)},
                {".error_by_binned_length_homopolymer_and_following_ref_base", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("G,T,SHORT_HOMOPOLYMER", 6L, 1L)},
                {".error_by_read_ordinality_and_pre_dinuc", chrMReadsWithClips, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("FIRST,A,A", 3L, 3L)},
                // Using a sam file with a single error it is easy to validate demonstrate these tests should pass
                {".error_by_all", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("all", 72L, 1L)},
                // There are two base qualities in the bam, the error occurs in quality "32"
                {".error_by_base_quality", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("32", 51L, 1L)},
                // There are two base qualities in the bam, the error occurs in quality "32", make sure we detect no errors in quality "33"
                {".error_by_base_quality", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("33", 21L, 0L)},
                // After a homopolymer run of 2 Gs, ending with a T there is a single error (curious if this 3 really should be 2, do we have an off by one?)
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("G,T,3", 2L, 1L)},

                // simpleSamWithBaseErrors2 contains 2 differences from the reference
                // after 2 different homopolymers.
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("T,A,3", 2L, 1L)},
                // Make sure that we can correctly identify an error after a homopolymer
                {".error_by_homopolymer_and_following_ref_base", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("C,G,4", 1L, 1L)},
                {".error_by_read_ordinality_and_pre_dinuc", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("FIRST,G,T", 1L, 1L)},
                // The covariate "0.5" was chosen to avoid those that have repeating decimals and could have alternative representations as
                // a string.  GC is calculated over the entire read including clipped bases, while errors are calculated only over unclipped bases.
                {".error_by_gc", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("0.5", 36L, 1L)},
                // Make sure that we can detect errors at a particular cycle
                {".error_by_cycle", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("10", 2L, 1L)},
                // There should be one error in the read with mapping quality 60.
                {".error_by_mapping_quality", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("60", 36L, 1L)},
                // There should be no errors in the read with mapping quality 0.
                {".error_by_mapping_quality", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("50", 36L, 0L)},
                // One base has an error in the read group 62A40.2
                {".error_by_read_group", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("62A40.2", 72L, 1L)},
                // No additional mismatches are found on the read with 1 mismatch.
                {".error_by_mismatches_in_read", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("1", 35L, 0L)},
                // No additional mismatches are found on the read with 1 mismatch. (Just another way to check)
                {".error_by_mismatches_in_read", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("0", 37L, 1L)},
                // There should be no errors in the CAG context because it matches reference
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CAG", 1L, 0L)},
                // There should be one error in the GTC context
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("GTC", 1L, 1L)},
                // There should be one error in the CTT context
                {".error_by_one_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CCT", 3L, 0L)},
                // There should be no errors in the ACGGG context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("ACGGG", 1L, 0L)},
                // There should be one error in the GGTCT context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("GGTCT", 1L, 1L)},
                // There should be no errors in the CTTGA context (appears in sam file as TCAaG in second read)
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CTTGA", 1L, 0L)},
                // There should be one error in the CCGTG context
                {".error_by_two_base_padded_context", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CCGTG", 1L, 1L)},
                // Reads that don't have consensus tags should be stratified as UNKNOWN
                {".error_by_consensus", simpleSamWithBaseErrors1, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("UNKNOWN", 72L, 1L)},
                // There should be 2 errors in the one simplex singleton reads.
                {".error_by_consensus", simpleSingleStrandConsensusSamWithBaseErrors, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("SIMPLEX_SINGLETON", 36L, 2L)},
                // There should be no errors in the simplex consensus read.
                {".error_by_consensus", simpleSingleStrandConsensusSamWithBaseErrors, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("SIMPLEX_CONSENSUS", 36L, 0L)},
                // There should be one error in duplex singleton read.  Also the N in this read reduces total bases from 36 to 35.
                {".error_by_consensus", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("DUPLEX_SINGLETON", 35L, 1L)},
                // There should be two errors in the duplex consensus read.
                {".error_by_consensus", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("DUPLEX_CONSENSUS", 36L, 2L)},
                // There should be two errors in the read with no Ns.
                {".error_by_ns_in_read", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("0", 36L, 2L)},
                // There should be one errors in the read with one N.
                {".error_by_ns_in_read", simpleDuplexConsensusSamWithBaseErrors, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("1", 35L, 1L)},
                // There are two errors, one should show up in CYCLE_QUINTILE_1 and the other in CYCLE_QUINTILE_3
                // CYCLE_QUINTILE_5 has 16 total (2 more than the other bins) bases due to rounding.
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CYCLE_QUINTILE_1", 14L, 1L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CYCLE_QUINTILE_2", 14L, 0L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CYCLE_QUINTILE_3", 14L, 1L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CYCLE_QUINTILE_4", 14L, 0L)},
                {".error_by_binned_cycle", simpleSamWithBaseErrors2, priorQ,
                        new ErrorMetrics.SimpleErrorMetric("CYCLE_QUINTILE_5", 16L, 0L)}
        };
    }

    @Test(dataProvider = "OneCovariateErrorMetricsDataProvider")
    public void testOneCovariateErrorMetrics(final String errorSubscript, final File samFile, final int priorQ, ErrorMetrics.SimpleErrorMetric expectedMetric) {
        final String input = samFile.toString();
        final String referenceFile = TEST_DIR + "chrM.reference.fasta";
        final String vcf = TEST_DIR + "NIST.selected.vcf";

        final File outputBaseFileName = new File(OUTPUT_DATA_PATH, "test");
        final File errorByAll = new File(outputBaseFileName.getAbsolutePath() + errorSubscript);
        errorByAll.deleteOnExit();
        outputBaseFileName.deleteOnExit();

        final String[] args = {
                "INPUT=" + input,
                "OUTPUT=" + outputBaseFileName,
                "REFERENCE_SEQUENCE=" + referenceFile,
                "ERROR_METRICS=" + "ERROR,TWO_BASE_PADDED_CONTEXT", // Not all covariates are included by default, but we still want to test them.
                "ERROR_METRICS=" + "ERROR,CONSENSUS",
                "ERROR_METRICS=" + "ERROR,NS_IN_READ",
                "ERROR_METRICS=" + "ERROR,BINNED_CYCLE",
                "VCF=" + vcf
        };

        Assert.assertEquals(new CollectBamErrorMetrics().instanceMain(args), 0);

        ErrorMetrics.setPriorError(QualityUtil.getErrorProbabilityFromPhredScore(priorQ));
        expectedMetric.calculateDerivedFields();

        // Note that soft clipped bases are not counted
        List<ErrorMetrics.SimpleErrorMetric> metrics = MetricsFile.readBeans(errorByAll);

        ErrorMetrics.SimpleErrorMetric metric = metrics
                .stream()
                .filter(m -> m.COVARIATE.equals(expectedMetric.COVARIATE))
                .findAny()
                .orElseThrow(() -> new AssertionError("didn't find metric with COVARIATE==" + expectedMetric.COVARIATE));

        Assert.assertEquals(metric, expectedMetric);
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
                {1.0 + Math.ulp(1.0)}, {0.0 - Math.ulp(0.0)}, {-1.0}, {-0.1}, {1.1}, {100.0}
        };
    }

    @Test(dataProvider = "readCycleBinDataError", expectedExceptions = IllegalArgumentException.class)
    public void testReadCycleBinError(final double relativePosition) {
        ReadBaseStratification.CycleBin.valueOf(relativePosition);
    }
}

