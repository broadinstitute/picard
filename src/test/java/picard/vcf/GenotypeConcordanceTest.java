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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.vcf.GenotypeConcordanceStates.CallState;
import picard.vcf.GenotypeConcordanceStates.TruthAndCallStates;
import picard.vcf.GenotypeConcordanceStates.TruthState;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class GenotypeConcordanceTest {

    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("GenotypeConcordanceTest", null);
    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");

    // Test VCFs
    private static final File CEU_TRIOS_SNPS_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
    private static final File CEU_TRIOS_INDELS_VCF = new File(TEST_DATA_PATH, "CEUTrio-indels.vcf");

    // Test that missing sites flag for new scheme works for NIST data sets
    private static final File NIST_MISSING_SITES_TRUTH_VCF = new File(TEST_DATA_PATH, "NIST.selected.vcf");

    // Test that we notice a difference on the first line
    private static final File CEU_TRIOS_SNPS_FIRST_LINE_DIFF_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_first_line_diff.vcf");

    // Test that we notice a difference on the last line
    private static final File CEU_TRIOS_SNPS_LAST_LINE_DIFF_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_last_line_diff.vcf");

    // Test that we notice a deleted line
    private static final File CEU_TRIOS_SNPS_DEL_LINE_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_del_line.vcf");

    // Existing/expected base metrics file names
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff";
    private static final String CEU_TRIOS_INDELS_VS_CEU_TRIOS_INDELS_GC = "CEUTrio-indels_vs_CEUTrio-indels_GtConcordanceDiff";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_FIRST_LINE_DIFF_GC = "CEUTrio-snps_CEUTrio-snps_first_line_GtConcordanceDiff";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_LAST_LINE_DIFF_GC = "CEUTrio-snps_CEUTrio-snps_last_line_GtConcordanceDiff";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_DEL_LINE_GC = "CEUTrio-snps_CEUTrio-snps_del_line_GtConcordanceDiff";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_ALL_ROWS = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_AllRows";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_GQ = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_MinGq";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_DP = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_MinDp";
    private static final String NIST_TRUTH_SNPS_VS_CEU_TRIOS_SNPS_GC = "NIST-truth-snps_vs_CEUTrio-snps_GtConcordanceDiff";

    private static final String TRUTH_SAMPLE_NAME = "Foo";
    private static final String CALL_SAMPLE_NAME = "Foo";

    // A [ref] / T at 10
    private static final String snpLoc = "chr1";
    private static final int snpLocStart = 10;
    private static final int snpLocStop = 10;

    private static final Allele Aref = Allele.create("A", true);
    private static final Allele C = Allele.create("C");
    private static final Allele G = Allele.create("G");
    private static final Allele T = Allele.create("T");

    private static final Allele AA = Allele.create("AA");
    private static final Allele AAA = Allele.create("AAA");
    private static final Allele AAAA = Allele.create("AAAA");
    private static final Allele AAAAA = Allele.create("AAAAA");

    private static final String NORMALIZE_ALLELES_TRUTH = "normalize_alleles_truth.vcf";
    private static final String NORMALIZE_ALLELES_CALL = "normalize_alleles_call.vcf";
    private static final String NORMALIZE_NO_CALLS_TRUTH = "normalize_no_calls_truth.vcf";
    private static final String NORMALIZE_NO_CALLS_CALL = "normalize_no_calls_call.vcf";

    @AfterClass
    public void tearDown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider(name = "genotypeConcordanceTestFileData")
    public Object[][] getGenotypeConcordanceTestFileData() {
        return new Object[][]{
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12878", null, null, false, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC},
                {CEU_TRIOS_INDELS_VCF, "NA12878", CEU_TRIOS_INDELS_VCF, "NA12878", null, null, false, false, CEU_TRIOS_INDELS_VS_CEU_TRIOS_INDELS_GC},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_FIRST_LINE_DIFF_VCF, "NA12878", null, null, false, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_FIRST_LINE_DIFF_GC},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_LAST_LINE_DIFF_VCF, "NA12878", null, null, false, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_LAST_LINE_DIFF_GC},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_DEL_LINE_VCF, "NA12878", null, null, false, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_DEL_LINE_GC},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12878", null, null, true, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_ALL_ROWS},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12891", 40, null, false, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_GQ},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12891", null, 40, false, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_DP},
                {NIST_MISSING_SITES_TRUTH_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12878", null, null, false, true, NIST_TRUTH_SNPS_VS_CEU_TRIOS_SNPS_GC}
        };
    }

    @Test(dataProvider = "genotypeConcordanceTestFileData")
    public void testGenotypeConcordance(final File vcf1, final String sample1, final File vcf2, final String sample2,
                                        final Integer minGq, final Integer minDp, final boolean outputAllRows, final boolean missingSitesFlag,
                                        final String expectedOutputFileBaseName) throws Exception {
        final List<Boolean> withVcfs = Arrays.asList(true, false);
        for (final boolean withVcf : withVcfs) {
            final File outputBaseFileName    = new File(OUTPUT_DATA_PATH, "actualGtConc");
            final File outputSummaryFile     = new File(outputBaseFileName.getAbsolutePath() + GenotypeConcordance.SUMMARY_METRICS_FILE_EXTENSION);
            final File outputDetailsFile     = new File(outputBaseFileName.getAbsolutePath() + GenotypeConcordance.DETAILED_METRICS_FILE_EXTENSION);
            final File outputContingencyFile = new File(outputBaseFileName.getAbsolutePath() + GenotypeConcordance.CONTINGENCY_METRICS_FILE_EXTENSION);
            final Path outputVcfFile         = Paths.get(outputBaseFileName.getAbsolutePath() + GenotypeConcordance.OUTPUT_VCF_FILE_EXTENSION);
            outputSummaryFile.deleteOnExit();
            outputDetailsFile.deleteOnExit();
            outputContingencyFile.deleteOnExit();
            outputVcfFile.toFile().deleteOnExit();

            final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
            genotypeConcordance.TRUTH_VCF = vcf1;
            genotypeConcordance.TRUTH_SAMPLE = sample1;
            genotypeConcordance.CALL_VCF = vcf2;
            genotypeConcordance.CALL_SAMPLE = sample2;
            if (minGq != null) genotypeConcordance.MIN_GQ = minGq;
            if (minDp != null) genotypeConcordance.MIN_DP = minDp;
            genotypeConcordance.OUTPUT_ALL_ROWS = outputAllRows;
            genotypeConcordance.OUTPUT = outputBaseFileName;
            genotypeConcordance.MISSING_SITES_HOM_REF = missingSitesFlag;
            if (missingSitesFlag) {
                genotypeConcordance.INTERVALS = Collections.singletonList(new File(TEST_DATA_PATH, "IntervalList1PerChrom.interval_list"));
            }
            genotypeConcordance.OUTPUT_VCF = withVcf;

            Assert.assertEquals(genotypeConcordance.instanceMain(new String[0]), 0);
            assertMetricsFileEqual(outputSummaryFile, new File(TEST_DATA_PATH, expectedOutputFileBaseName + GenotypeConcordance.SUMMARY_METRICS_FILE_EXTENSION));
            assertMetricsFileEqual(outputDetailsFile, new File(TEST_DATA_PATH, expectedOutputFileBaseName + GenotypeConcordance.DETAILED_METRICS_FILE_EXTENSION));
            assertMetricsFileEqual(outputContingencyFile, new File(TEST_DATA_PATH, expectedOutputFileBaseName + GenotypeConcordance.CONTINGENCY_METRICS_FILE_EXTENSION));

            if (withVcf) {
                // An ugly way to compare VCFs
                final Path expectedVcf = Paths.get(TEST_DATA_PATH.getAbsolutePath(), expectedOutputFileBaseName + ".vcf");
                final BufferedLineReader reader = new BufferedLineReader(new GZIPInputStream(new FileInputStream(outputVcfFile.toFile())));
                final Iterator<String> actualLines;
                {
                    final List<String> lines = new ArrayList<>();
                    while (-1 != reader.peek()) {
                        lines.add(reader.readLine());
                    }
                    reader.close();
                    actualLines = lines.iterator();
                }

                final Iterator<String> expectedLines = Files.lines(expectedVcf).iterator();
                while (actualLines.hasNext() && expectedLines.hasNext()) {
                    final String actualLine = actualLines.next();
                    final String expectedLine = expectedLines.next();
                    Assert.assertEquals(actualLine, expectedLine);
                }
                Assert.assertFalse(actualLines.hasNext());
                Assert.assertFalse(expectedLines.hasNext());
            }
        }
    }

    private void assertMetricsFileEqual(final File actualMetricsFile, final File expectedMetricsFile) throws FileNotFoundException {
        // Actual metrics file
        final MetricsFile<GenotypeConcordanceSummaryMetrics, Comparable<?>> actual = new MetricsFile<GenotypeConcordanceSummaryMetrics, Comparable<?>>();
        actual.read(new FileReader(actualMetricsFile));

        // Expected metrics file
        final MetricsFile<GenotypeConcordanceSummaryMetrics, Comparable<?>> expected = new MetricsFile<GenotypeConcordanceSummaryMetrics, Comparable<?>>();
        expected.read(new FileReader(expectedMetricsFile));

        // Note - cannot use .equals as it calls .areHeadersEqual and they are not since the timestamp (at a minimum is different)
        Assert.assertTrue(expected.areMetricsEqual(actual));
        Assert.assertTrue(expected.areHistogramsEqual(actual));
    }

    public static GenotypeConcordanceCounts getGenotypeConcordanceCounts(final File truthVCF, final File callVCF, final String callSample, final boolean missingSitesFlag, List<File> intervalFiles){
        //gets the Genotype Concordance Counts for the detail tests for each scheme
        final File outputBaseFileName = new File(OUTPUT_DATA_PATH, "actualGtConc");
        final File outputSummaryFile = new File(outputBaseFileName.getAbsolutePath() + GenotypeConcordance.SUMMARY_METRICS_FILE_EXTENSION);
        final File outputDetailsFile = new File(outputBaseFileName.getAbsolutePath() + GenotypeConcordance.DETAILED_METRICS_FILE_EXTENSION);
        outputSummaryFile.deleteOnExit();
        outputDetailsFile.deleteOnExit();

        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = truthVCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = callVCF;
        genotypeConcordance.CALL_SAMPLE = callSample;
        genotypeConcordance.MISSING_SITES_HOM_REF = missingSitesFlag;
        genotypeConcordance.INTERVALS = intervalFiles;
        genotypeConcordance.OUTPUT = outputBaseFileName;

        Assert.assertEquals(genotypeConcordance.instanceMain(new String[0]), 0);

        return genotypeConcordance.getSnpCounter();
    }

    public static void assertNonZeroCountsAgree(final GenotypeConcordanceCounts counter, final Map<TruthAndCallStates, Integer> expectedCountMap) {
        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                Integer expectedCount = expectedCountMap.get(new TruthAndCallStates(truthState, callState));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(counter.getCount(truthState, callState), expectedCount.intValue());
            }
        }
    }

    @DataProvider(name = "genotypeConcordanceDetermineStateDataProvider")
    public Object[][] genotypeConcordanceDetermineStateDataProvider() {
        final Object[][] originalUnitTestData = new Object[][]{
                {Aref, Aref, TruthState.HOM_REF, Aref, Aref, CallState.HOM_REF},

                {Aref, Aref, TruthState.HOM_REF, Aref, C, CallState.HET_REF_VAR1},
                {Aref, Aref, TruthState.HOM_REF, Aref, G, CallState.HET_REF_VAR1},
                {Aref, Aref, TruthState.HOM_REF, Aref, T, CallState.HET_REF_VAR1},

                {Aref, Aref, TruthState.HOM_REF, C, G, CallState.HET_VAR1_VAR2},
                {Aref, Aref, TruthState.HOM_REF, C, T, CallState.HET_VAR1_VAR2},
                {Aref, Aref, TruthState.HOM_REF, G, T, CallState.HET_VAR1_VAR2},

                {Aref, Aref, TruthState.HOM_REF, C, C, CallState.HOM_VAR1},
                {Aref, Aref, TruthState.HOM_REF, G, G, CallState.HOM_VAR1},
                {Aref, Aref, TruthState.HOM_REF, T, T, CallState.HOM_VAR1},

                //---
                {Aref, C, TruthState.HET_REF_VAR1, Aref, Aref, CallState.HOM_REF},
                {Aref, G, TruthState.HET_REF_VAR1, Aref, Aref, CallState.HOM_REF},
                {Aref, T, TruthState.HET_REF_VAR1, Aref, Aref, CallState.HOM_REF},

                {Aref, C, TruthState.HET_REF_VAR1, Aref, C, CallState.HET_REF_VAR1},
                {Aref, C, TruthState.HET_REF_VAR1, Aref, G, CallState.HET_REF_VAR2},
                {Aref, C, TruthState.HET_REF_VAR1, Aref, T, CallState.HET_REF_VAR2},
                {Aref, G, TruthState.HET_REF_VAR1, Aref, C, CallState.HET_REF_VAR2},
                {Aref, G, TruthState.HET_REF_VAR1, Aref, G, CallState.HET_REF_VAR1},
                {Aref, G, TruthState.HET_REF_VAR1, Aref, T, CallState.HET_REF_VAR2},
                {Aref, T, TruthState.HET_REF_VAR1, Aref, C, CallState.HET_REF_VAR2},
                {Aref, T, TruthState.HET_REF_VAR1, Aref, G, CallState.HET_REF_VAR2},
                {Aref, T, TruthState.HET_REF_VAR1, Aref, T, CallState.HET_REF_VAR1},

                {Aref, C, TruthState.HET_REF_VAR1, C, G, CallState.HET_VAR1_VAR2},
                {Aref, C, TruthState.HET_REF_VAR1, C, T, CallState.HET_VAR1_VAR2},
                {Aref, C, TruthState.HET_REF_VAR1, G, T, CallState.HET_VAR3_VAR4},  // Why isn't this called HET_VAR2_VAR3???
                {Aref, G, TruthState.HET_REF_VAR1, C, G, CallState.HET_VAR1_VAR2},
                {Aref, G, TruthState.HET_REF_VAR1, C, T, CallState.HET_VAR3_VAR4},
                {Aref, G, TruthState.HET_REF_VAR1, G, T, CallState.HET_VAR1_VAR2},
                {Aref, T, TruthState.HET_REF_VAR1, C, G, CallState.HET_VAR3_VAR4},
                {Aref, T, TruthState.HET_REF_VAR1, C, T, CallState.HET_VAR1_VAR2},
                {Aref, T, TruthState.HET_REF_VAR1, G, T, CallState.HET_VAR1_VAR2},

                {Aref, C, TruthState.HET_REF_VAR1, C, C, CallState.HOM_VAR1},
                {Aref, C, TruthState.HET_REF_VAR1, G, G, CallState.HOM_VAR2},
                {Aref, C, TruthState.HET_REF_VAR1, T, T, CallState.HOM_VAR2},
                {Aref, G, TruthState.HET_REF_VAR1, C, C, CallState.HOM_VAR2},
                {Aref, G, TruthState.HET_REF_VAR1, G, G, CallState.HOM_VAR1},
                {Aref, G, TruthState.HET_REF_VAR1, T, T, CallState.HOM_VAR2},
                {Aref, T, TruthState.HET_REF_VAR1, C, C, CallState.HOM_VAR2},
                {Aref, T, TruthState.HET_REF_VAR1, G, G, CallState.HOM_VAR2},
                {Aref, T, TruthState.HET_REF_VAR1, T, T, CallState.HOM_VAR1},

                //---
                {C, G, TruthState.HET_VAR1_VAR2, Aref, Aref, CallState.HOM_REF},
                {C, T, TruthState.HET_VAR1_VAR2, Aref, Aref, CallState.HOM_REF},
                {G, T, TruthState.HET_VAR1_VAR2, Aref, Aref, CallState.HOM_REF},

                {C, G, TruthState.HET_VAR1_VAR2, Aref, C, CallState.HET_REF_VAR1},
                {C, G, TruthState.HET_VAR1_VAR2, Aref, G, CallState.HET_REF_VAR1},
                {C, G, TruthState.HET_VAR1_VAR2, Aref, T, CallState.HET_REF_VAR3},
                {C, T, TruthState.HET_VAR1_VAR2, Aref, C, CallState.HET_REF_VAR1},
                {C, T, TruthState.HET_VAR1_VAR2, Aref, G, CallState.HET_REF_VAR3},
                {C, T, TruthState.HET_VAR1_VAR2, Aref, T, CallState.HET_REF_VAR1},
                {G, T, TruthState.HET_VAR1_VAR2, Aref, C, CallState.HET_REF_VAR3},
                {G, T, TruthState.HET_VAR1_VAR2, Aref, G, CallState.HET_REF_VAR1},
                {G, T, TruthState.HET_VAR1_VAR2, Aref, T, CallState.HET_REF_VAR1},

                {C, G, TruthState.HET_VAR1_VAR2, C, C, CallState.HOM_VAR1},
                {C, G, TruthState.HET_VAR1_VAR2, G, G, CallState.HOM_VAR1},
                {C, G, TruthState.HET_VAR1_VAR2, T, T, CallState.HOM_VAR3},
                {C, T, TruthState.HET_VAR1_VAR2, C, C, CallState.HOM_VAR1},
                {C, T, TruthState.HET_VAR1_VAR2, G, G, CallState.HOM_VAR3},
                {C, T, TruthState.HET_VAR1_VAR2, T, T, CallState.HOM_VAR1},
                {G, T, TruthState.HET_VAR1_VAR2, C, C, CallState.HOM_VAR3},
                {G, T, TruthState.HET_VAR1_VAR2, G, G, CallState.HOM_VAR1},
                {G, T, TruthState.HET_VAR1_VAR2, T, T, CallState.HOM_VAR1},

                {C, G, TruthState.HET_VAR1_VAR2, C, G, CallState.HET_VAR1_VAR2},
                {C, G, TruthState.HET_VAR1_VAR2, C, T, CallState.HET_VAR1_VAR3},
                {C, G, TruthState.HET_VAR1_VAR2, G, T, CallState.HET_VAR1_VAR3},
                {C, T, TruthState.HET_VAR1_VAR2, C, G, CallState.HET_VAR1_VAR3},
                {C, T, TruthState.HET_VAR1_VAR2, C, T, CallState.HET_VAR1_VAR2},
                {C, T, TruthState.HET_VAR1_VAR2, G, T, CallState.HET_VAR1_VAR3},
                {G, T, TruthState.HET_VAR1_VAR2, C, G, CallState.HET_VAR1_VAR3},
                {G, T, TruthState.HET_VAR1_VAR2, C, T, CallState.HET_VAR1_VAR3},
                {G, T, TruthState.HET_VAR1_VAR2, G, T, CallState.HET_VAR1_VAR2},

                //---
                {C, C, TruthState.HOM_VAR1, Aref, Aref, CallState.HOM_REF},
                {G, G, TruthState.HOM_VAR1, Aref, Aref, CallState.HOM_REF},
                {T, T, TruthState.HOM_VAR1, Aref, Aref, CallState.HOM_REF},

                {C, C, TruthState.HOM_VAR1, Aref, C, CallState.HET_REF_VAR1},
                {C, C, TruthState.HOM_VAR1, Aref, G, CallState.HET_REF_VAR2},
                {C, C, TruthState.HOM_VAR1, Aref, T, CallState.HET_REF_VAR2},
                {G, G, TruthState.HOM_VAR1, Aref, C, CallState.HET_REF_VAR2},
                {G, G, TruthState.HOM_VAR1, Aref, G, CallState.HET_REF_VAR1},
                {G, G, TruthState.HOM_VAR1, Aref, T, CallState.HET_REF_VAR2},
                {T, T, TruthState.HOM_VAR1, Aref, C, CallState.HET_REF_VAR2},
                {T, T, TruthState.HOM_VAR1, Aref, G, CallState.HET_REF_VAR2},
                {T, T, TruthState.HOM_VAR1, Aref, T, CallState.HET_REF_VAR1},

                {C, C, TruthState.HOM_VAR1, C, C, CallState.HOM_VAR1},
                {C, C, TruthState.HOM_VAR1, G, G, CallState.HOM_VAR2},
                {C, C, TruthState.HOM_VAR1, T, T, CallState.HOM_VAR2},
                {G, G, TruthState.HOM_VAR1, C, C, CallState.HOM_VAR2},
                {G, G, TruthState.HOM_VAR1, G, G, CallState.HOM_VAR1},
                {G, G, TruthState.HOM_VAR1, T, T, CallState.HOM_VAR2},
                {T, T, TruthState.HOM_VAR1, C, C, CallState.HOM_VAR2},
                {T, T, TruthState.HOM_VAR1, G, G, CallState.HOM_VAR2},
                {T, T, TruthState.HOM_VAR1, T, T, CallState.HOM_VAR1},

                {C, C, TruthState.HOM_VAR1, C, G, CallState.HET_VAR1_VAR2},
                {C, C, TruthState.HOM_VAR1, C, T, CallState.HET_VAR1_VAR2},
                {C, C, TruthState.HOM_VAR1, G, T, CallState.HET_VAR3_VAR4},
                {G, G, TruthState.HOM_VAR1, C, G, CallState.HET_VAR1_VAR2},
                {G, G, TruthState.HOM_VAR1, C, T, CallState.HET_VAR3_VAR4},
                {G, G, TruthState.HOM_VAR1, G, T, CallState.HET_VAR1_VAR2},
                {T, T, TruthState.HOM_VAR1, C, G, CallState.HET_VAR3_VAR4},
                {T, T, TruthState.HOM_VAR1, C, T, CallState.HET_VAR1_VAR2},
                {T, T, TruthState.HOM_VAR1, G, T, CallState.HET_VAR1_VAR2},
                // Some Indel Cases
                {AA, AA, TruthState.HOM_VAR1, AAAA, AAAAA, CallState.HET_VAR3_VAR4},
                {AA, AAA, TruthState.HET_VAR1_VAR2, AAAA, AAAAA, CallState.HET_VAR3_VAR4},

                // Mixed Cases
                {C, AA, TruthState.IS_MIXED, AAAA, AAAAA, CallState.HET_VAR1_VAR2},
                {AA, C, TruthState.IS_MIXED, AAAA, AAAAA, CallState.HET_VAR1_VAR2},

                {AA, AA, TruthState.HOM_VAR1, C, AAAAA, CallState.IS_MIXED},
                {AA, AAA, TruthState.HET_VAR1_VAR2, AAAA, C, CallState.IS_MIXED},

                // No Call cases
                {Allele.NO_CALL, Aref, TruthState.NO_CALL, Aref, Aref, CallState.HOM_REF},
                {Aref, Allele.NO_CALL, TruthState.NO_CALL, Aref, Aref, CallState.HOM_REF},
                {Allele.NO_CALL, Allele.NO_CALL, TruthState.NO_CALL, Aref, Aref, CallState.HOM_REF},

                {Aref, Aref, TruthState.HOM_REF, Allele.NO_CALL, Aref, CallState.NO_CALL},
                {Aref, Aref, TruthState.HOM_REF, Aref, Allele.NO_CALL, CallState.NO_CALL},
                {Aref, Aref, TruthState.HOM_REF, Allele.NO_CALL, Allele.NO_CALL, CallState.NO_CALL}
        };
        // Rebuild a new set of unit test data with all permutations of alleles.
        final List<Object[]> allPermutationUnitTestDataList = new ArrayList<Object[]>();
        for (final Object[] unitTestData : originalUnitTestData) {
            allPermutationUnitTestDataList.add(unitTestData);
            final Allele truthAllele1 = (Allele) unitTestData[0];
            final Allele truthAllele2 = (Allele) unitTestData[1];
            final TruthState expectedTruthState = (TruthState) unitTestData[2];
            final Allele callAllele1 = (Allele) unitTestData[3];
            final Allele callAllele2 = (Allele) unitTestData[4];
            final CallState expectedCallState = (CallState) unitTestData[5];
            if (!callAllele1.equals(callAllele2)) {
                allPermutationUnitTestDataList.add(new Object[]{truthAllele1, truthAllele2, expectedTruthState, callAllele2, callAllele1, expectedCallState});
            }
            if (!truthAllele1.equals(truthAllele2)) {
                allPermutationUnitTestDataList.add(new Object[]{truthAllele2, truthAllele1, expectedTruthState, callAllele1, callAllele2, expectedCallState});
                if (!callAllele1.equals(callAllele2)) {
                    allPermutationUnitTestDataList.add(new Object[]{truthAllele2, truthAllele1, expectedTruthState, callAllele2, callAllele1, expectedCallState});
                }
            }
        }
        Object[][] allPermutationUnitTestData = new Object[allPermutationUnitTestDataList.size()][];
        allPermutationUnitTestData = allPermutationUnitTestDataList.toArray(allPermutationUnitTestData);
        return allPermutationUnitTestData;
    }


    @Test(dataProvider = "genotypeConcordanceDetermineStateDataProvider")
    public void testGenotypeConcordanceDetermineState(final Allele truthAllele1, final Allele truthAllele2, final TruthState expectedTruthState,
                                                      final Allele callAllele1, final Allele callAllele2, final CallState expectedCallState) throws Exception {
        final List<Allele> truthAlleles = makeUniqueListOfAlleles(truthAllele1, truthAllele2);
        final Genotype truthGt = GenotypeBuilder.create(TRUTH_SAMPLE_NAME, Arrays.asList(truthAllele1, truthAllele2));

        final VariantContext truthVariantContext = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, truthAlleles).genotypes(truthGt).make();

        final List<Allele> callAlleles = makeUniqueListOfAlleles(callAllele1, callAllele2);
        final Genotype callGt = GenotypeBuilder.create(CALL_SAMPLE_NAME, Arrays.asList(callAllele1, callAllele2));
        final VariantContext callVariantContext = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, callAlleles).genotypes(callGt).make();

        testGenotypeConcordanceDetermineState(truthVariantContext, expectedTruthState, callVariantContext, expectedCallState, 0, 0);
    }

    @Test
    public void testGenotypeConcordanceDetermineStateNull() throws Exception {
        final List<Allele> alleles = makeUniqueListOfAlleles(Aref, C);
        final Genotype gt1 = GenotypeBuilder.create(TRUTH_SAMPLE_NAME, Arrays.asList(Aref, C));
        final VariantContext vc1 = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(gt1).make();

        testGenotypeConcordanceDetermineState(null, TruthState.MISSING, null, CallState.MISSING, 0, 0);
        testGenotypeConcordanceDetermineState(vc1, TruthState.HET_REF_VAR1, null, CallState.MISSING, 0, 0);
        testGenotypeConcordanceDetermineState(null, TruthState.MISSING, vc1, CallState.HET_REF_VAR1, 0, 0);
    }

    @Test
    public void testGenotypeConcordanceDetermineStateFilter() throws Exception {
        final Set<String> filters = new HashSet<String>(Arrays.asList("BAD!"));

        // Filtering on the variant context
        final List<Allele> alleles1 = makeUniqueListOfAlleles(Aref, C);
        final Genotype gt1 = GenotypeBuilder.create(TRUTH_SAMPLE_NAME, Arrays.asList(Aref, C));
        final VariantContext vcFiltered = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles1).genotypes(gt1).filters(filters).make();

        final List<Allele> alleles2 = makeUniqueListOfAlleles(Aref, T);
        final Genotype gt2 = GenotypeBuilder.create(TRUTH_SAMPLE_NAME, Arrays.asList(Aref, T));
        final VariantContext vcNotFiltered = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles2).genotypes(gt2).make();

        testGenotypeConcordanceDetermineState(vcFiltered, TruthState.VC_FILTERED, vcNotFiltered, CallState.HET_REF_VAR1, 0, 0);
        testGenotypeConcordanceDetermineState(vcNotFiltered, TruthState.HET_REF_VAR1, vcFiltered, CallState.VC_FILTERED, 0, 0);
        testGenotypeConcordanceDetermineState(vcFiltered, TruthState.VC_FILTERED, vcFiltered, CallState.VC_FILTERED, 0, 0);

        // Filtering on the genotype
        final List<String> gtFilters = new ArrayList<String>(Arrays.asList("WICKED"));
        final List<Allele> alleles3 = makeUniqueListOfAlleles(Aref, C);
        final Genotype gt3 = new GenotypeBuilder(TRUTH_SAMPLE_NAME, Arrays.asList(Aref, C)).filters(gtFilters).make();
        final VariantContext vcGtFiltered = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles3).genotypes(gt3).make();

        testGenotypeConcordanceDetermineState(vcGtFiltered, TruthState.GT_FILTERED, vcNotFiltered, CallState.HET_REF_VAR1, 0, 0);
        testGenotypeConcordanceDetermineState(vcNotFiltered, TruthState.HET_REF_VAR1, vcGtFiltered, CallState.GT_FILTERED, 0, 0);
        testGenotypeConcordanceDetermineState(vcGtFiltered, TruthState.GT_FILTERED, vcGtFiltered, CallState.GT_FILTERED, 0, 0);
    }

    @Test
    public void testGenotypeConcordanceDetermineStateDp() throws Exception {
        final List<Allele> allelesNormal = makeUniqueListOfAlleles(Aref, C);
        final Genotype gtNormal = GenotypeBuilder.create(TRUTH_SAMPLE_NAME, Arrays.asList(Aref, C));
        final VariantContext vcNormal = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, allelesNormal).genotypes(gtNormal).make();

        final List<Allele> allelesLowDp = makeUniqueListOfAlleles(Aref, C);
        final Genotype gtLowDp = new GenotypeBuilder(TRUTH_SAMPLE_NAME, Arrays.asList(Aref, C)).DP(4).make();
        final VariantContext vcLowDp = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, allelesLowDp).genotypes(gtLowDp).make();

        testGenotypeConcordanceDetermineState(vcLowDp, TruthState.LOW_DP, vcNormal, CallState.HET_REF_VAR1, 0, 20);
        testGenotypeConcordanceDetermineState(vcLowDp, TruthState.HET_REF_VAR1, vcLowDp, CallState.HET_REF_VAR1, 0, 2);

        testGenotypeConcordanceDetermineState(vcNormal, TruthState.HET_REF_VAR1, vcLowDp, CallState.LOW_DP, 0, 20);
        testGenotypeConcordanceDetermineState(vcNormal, TruthState.HET_REF_VAR1, vcLowDp, CallState.HET_REF_VAR1, 0, 2);

        testGenotypeConcordanceDetermineState(vcLowDp, TruthState.LOW_DP, vcLowDp, CallState.LOW_DP, 0, 20);
        testGenotypeConcordanceDetermineState(vcLowDp, TruthState.HET_REF_VAR1, vcLowDp, CallState.HET_REF_VAR1, 0, 2);
    }

    @Test
    public void testGenotypeConcordanceDetermineStateGq() throws Exception {
        final List<Allele> allelesNormal = makeUniqueListOfAlleles(Aref, C);
        final Genotype gtNormal = GenotypeBuilder.create(TRUTH_SAMPLE_NAME, Arrays.asList(Aref, C));
        final VariantContext vcNormal = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, allelesNormal).genotypes(gtNormal).make();

        final List<Allele> allelesLowGq = makeUniqueListOfAlleles(Aref, C);
        final Genotype gtLowGq = new GenotypeBuilder(TRUTH_SAMPLE_NAME, Arrays.asList(Aref, C)).GQ(4).make();
        final VariantContext vcLowGq = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, allelesLowGq).genotypes(gtLowGq).make();

        testGenotypeConcordanceDetermineState(vcLowGq, TruthState.LOW_GQ, vcNormal, CallState.HET_REF_VAR1, 20, 0);
        testGenotypeConcordanceDetermineState(vcLowGq, TruthState.HET_REF_VAR1, vcLowGq, CallState.HET_REF_VAR1, 2, 0);

        testGenotypeConcordanceDetermineState(vcNormal, TruthState.HET_REF_VAR1, vcLowGq, CallState.LOW_GQ, 20, 0);
        testGenotypeConcordanceDetermineState(vcNormal, TruthState.HET_REF_VAR1, vcLowGq, CallState.HET_REF_VAR1, 2, 0);

        testGenotypeConcordanceDetermineState(vcLowGq, TruthState.LOW_GQ, vcLowGq, CallState.LOW_GQ, 20, 0);
        testGenotypeConcordanceDetermineState(vcLowGq, TruthState.HET_REF_VAR1, vcLowGq, CallState.HET_REF_VAR1, 2, 0);
    }

    /**
     * Test method to determine that the expected truth and call states are returned for a pair of truth and call variant contexts.
     * @param truthVariantContext
     * @param expectedTruthState
     * @param callVariantContext
     * @param expectedCallState
     * @param minGq
     * @param minDp
     */
    private void testGenotypeConcordanceDetermineState(final VariantContext truthVariantContext, final TruthState expectedTruthState,
                                                       final VariantContext callVariantContext, final CallState expectedCallState,
                                                       final int minGq, final int minDp) {
        final TruthAndCallStates truthAndCallStates = GenotypeConcordance.determineState(truthVariantContext, TRUTH_SAMPLE_NAME,
                callVariantContext, CALL_SAMPLE_NAME, minGq, minDp);
        Assert.assertEquals(truthAndCallStates.truthState, expectedTruthState);
        Assert.assertEquals(truthAndCallStates.callState, expectedCallState);
    }

    /**
     * Simple method to return a list of unique alleles.
     */
    private List<Allele> makeUniqueListOfAlleles(final Allele... alleles) {
        final Set<Allele> uniqueAlleles = new HashSet<Allele>();
        for (final Allele allele : alleles) {
            if (!allele.equals(Allele.NO_CALL)) {
                uniqueAlleles.add(allele);
            }
        }
        if (!uniqueAlleles.contains(Aref)) uniqueAlleles.add(Aref);
        return new ArrayList<Allele>(uniqueAlleles);
    }

    /**
     * Tests that we normalize indels correctly
     */
    @Test
    public void testNormalizeAllelesForIndels() {

        final Path truthVcfPath = Paths.get(TEST_DATA_PATH.getAbsolutePath(), NORMALIZE_ALLELES_TRUTH);
        final Path callVcfPath  = Paths.get(TEST_DATA_PATH.getAbsolutePath(), NORMALIZE_ALLELES_CALL);

        final VCFFileReader truthReader = new VCFFileReader(truthVcfPath.toFile(), false);
        final VCFFileReader callReader  = new VCFFileReader(callVcfPath.toFile(), false);

        final Iterator<VariantContext> truthIterator = truthReader.iterator();
        final Iterator<VariantContext> callIterator = callReader.iterator();

        final String truthSample = truthReader.getFileHeader().getSampleNamesInOrder().get(0);
        final String callSample  = callReader.getFileHeader().getSampleNamesInOrder().get(0);

        while (truthIterator.hasNext()) {
            final VariantContext truthCtx = truthIterator.next();
            final VariantContext callCtx  = callIterator.next();

            {
                final GenotypeConcordance.Alleles alleles = GenotypeConcordance.normalizeAlleles(truthCtx, truthSample, callCtx, callSample);
                Assert.assertEquals(alleles.truthAllele1, alleles.callAllele1);
                Assert.assertEquals(alleles.truthAllele2, alleles.callAllele2);
            }
            {
                final GenotypeConcordance.Alleles alleles = GenotypeConcordance.normalizeAlleles(callCtx, callSample, truthCtx, truthSample);
                Assert.assertEquals(alleles.truthAllele1, alleles.callAllele1);
                Assert.assertEquals(alleles.truthAllele2, alleles.callAllele2);
            }
        }

        truthReader.close();
        callReader.close();
    }

    @Test
    public void testNoCallVariants() {
        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = new File(TEST_DATA_PATH, "mini.vcf");
        genotypeConcordance.TRUTH_SAMPLE = "NA20801";
        genotypeConcordance.CALL_VCF = new File(TEST_DATA_PATH, "mini.vcf");
        genotypeConcordance.CALL_SAMPLE = "NA19920";
        genotypeConcordance.OUTPUT = new File(OUTPUT_DATA_PATH, "TwoNoCalls");
        genotypeConcordance.OUTPUT_VCF = true;

        Assert.assertEquals(genotypeConcordance.instanceMain(new String[0]), 0);
    }

    @Test
    public void testNormalizeAllelesForWritingVCF() throws FileNotFoundException {
        final File truthVcfPath             = new File(TEST_DATA_PATH.getAbsolutePath(), NORMALIZE_NO_CALLS_TRUTH);
        final File callVcfPath              = new File(TEST_DATA_PATH.getAbsolutePath(), NORMALIZE_NO_CALLS_CALL);
        final File outputBaseFileName       = new File(OUTPUT_DATA_PATH, "MultipleRefAlleles");
        final File outputContingencyMetrics = new File(outputBaseFileName.getAbsolutePath() + GenotypeConcordance.CONTINGENCY_METRICS_FILE_EXTENSION);
        outputContingencyMetrics.deleteOnExit();

        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = truthVcfPath;
        genotypeConcordance.TRUTH_SAMPLE = "truth";
        genotypeConcordance.CALL_VCF = callVcfPath;
        genotypeConcordance.CALL_SAMPLE = "truth";
        genotypeConcordance.OUTPUT = new File(OUTPUT_DATA_PATH, "MultipleRefAlleles");
        genotypeConcordance.OUTPUT_VCF = true;

        Assert.assertEquals(genotypeConcordance.instanceMain(new String[0]), 0);

        final MetricsFile<GenotypeConcordanceContingencyMetrics, Comparable<?>> output = new MetricsFile<GenotypeConcordanceContingencyMetrics, Comparable<?>>();
        output.read(new FileReader(outputContingencyMetrics));

        for (final GenotypeConcordanceContingencyMetrics metrics : output.getMetrics()) {
            if(metrics.VARIANT_TYPE == VariantContext.Type.INDEL){
                Assert.assertEquals(metrics.TP_COUNT, 3);
                Assert.assertEquals(metrics.TN_COUNT, 3);
                Assert.assertEquals(metrics.FP_COUNT, 0);
                Assert.assertEquals(metrics.FN_COUNT, 0);
                Assert.assertEquals(metrics.EMPTY_COUNT, 2);
            }
        }
    }
}
