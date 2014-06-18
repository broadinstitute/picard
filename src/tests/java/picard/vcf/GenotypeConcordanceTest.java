package picard.vcf;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GenotypeConcordanceTest {

    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("GenotypeConcordanceTest", null);
    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");

    // Test VCFs
    private static final File CEU_TRIOS_SNPS_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
    private static final File CEU_TRIOS_INDELS_VCF = new File(TEST_DATA_PATH, "CEUTrio-indels.vcf");

    // Test that we notice a difference on the first line
    private static final File CEU_TRIOS_SNPS_FIRST_LINE_DIFF_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_first_line_diff.vcf");

    // Test that we notice a difference on the last line
    private static final File CEU_TRIOS_SNPS_LAST_LINE_DIFF_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_last_line_diff.vcf");

    // Test that we notice a deleted line
    private static final File CEU_TRIOS_SNPS_DEL_LINE_VCF = new File(TEST_DATA_PATH, "CEUTrio-snps_del_line.vcf");

    private static final String SUMMARY_METRICS_EXTENSION = ".summary_metrics.txt";
    private static final String DETAILED_METRICS_EXTENSION = ".detailed_metrics.txt";
    private static final String VARIANT_METRICS_EXTENSION = ".variant_metrics.txt";

    // Existing/expected base metrics file names
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff";
    private static final String CEU_TRIOS_INDELS_VS_CEU_TRIOS_INDELS_GC = "CEUTrio-indels_vs_CEUTrio-indels_GtConcordanceDiff";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_FIRST_LINE_DIFF_GC = "CEUTrio-snps_CEUTrio-snps_first_line_GtConcordanceDiff";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_LAST_LINE_DIFF_GC = "CEUTrio-snps_CEUTrio-snps_last_line_GtConcordanceDiff";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_DEL_LINE_GC = "CEUTrio-snps_CEUTrio-snps_del_line_GtConcordanceDiff";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_ALL_ROWS = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_AllRows";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_GQ = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_MinGq";
    private static final String CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_DP = "CEUTrio-snps_vs_CEUTrio-snps_GtConcordanceDiff_MinDp";

    private static final File INTERVALS_FILE = new File(TEST_DATA_PATH, "IntervalListChr1Small.interval_list");

    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @DataProvider(name = "genotypeConcordanceTestFileData")
    public Object[][] getGenotypeConcordanceTestFileData() {
        return new Object[][] {
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12878", null, null, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC},
                {CEU_TRIOS_INDELS_VCF, "NA12878", CEU_TRIOS_INDELS_VCF, "NA12878", null, null, false, CEU_TRIOS_INDELS_VS_CEU_TRIOS_INDELS_GC},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_FIRST_LINE_DIFF_VCF, "NA12878", null, null, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_FIRST_LINE_DIFF_GC},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_LAST_LINE_DIFF_VCF, "NA12878", null, null, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_LAST_LINE_DIFF_GC},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_DEL_LINE_VCF, "NA12878", null, null, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_DEL_LINE_GC},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12878", null, null, true, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_ALL_ROWS},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12891", 40, null, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_GQ},
                {CEU_TRIOS_SNPS_VCF, "NA12878", CEU_TRIOS_SNPS_VCF, "NA12891", null, 40, false, CEU_TRIOS_SNPS_VS_CEU_TRIOS_SNPS_GC_MIN_DP}
        };
    }

    @Test(dataProvider = "genotypeConcordanceTestFileData")
    public void testGenotypeConcordance(final File vcf1, final String sample1, final File vcf2, final String sample2,
                                        final Integer minGq, final Integer minDp, final boolean outputAllRows,
                                        final String expectedOutputFileBaseName) throws Exception {
        final File outputBaseFileName = new File(OUTPUT_DATA_PATH, "actualGtConc");
        final File outputSummaryFile = new File(outputBaseFileName.getAbsolutePath() + SUMMARY_METRICS_EXTENSION);
        final File outputDetailsFile = new File(outputBaseFileName.getAbsolutePath() + DETAILED_METRICS_EXTENSION);
        final File outputVariantFile = new File(outputBaseFileName.getAbsolutePath() + VARIANT_METRICS_EXTENSION);
        outputSummaryFile.deleteOnExit();
        outputDetailsFile.deleteOnExit();
        outputVariantFile.deleteOnExit();

        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = vcf1;
        genotypeConcordance.TRUTH_SAMPLE = sample1;
        genotypeConcordance.CALL_VCF = vcf2;
        genotypeConcordance.CALL_SAMPLE = sample2;
        if (minGq != null)
            genotypeConcordance.MIN_GQ = minGq;
        if (minDp != null)
            genotypeConcordance.MIN_DP = minDp;
        genotypeConcordance.OUTPUT_ALL_ROWS = outputAllRows;
        genotypeConcordance.OUTPUT = outputBaseFileName;
        genotypeConcordance.DEBUG_VARIANT_METRICS_FILE = outputVariantFile;
        final int returnCode = genotypeConcordance.instanceMain(new String[0]);

        genotypeConcordance.getSnpCounter();
        Assert.assertEquals(returnCode, 0);

        assertFilesEqual(outputSummaryFile, new File(TEST_DATA_PATH, expectedOutputFileBaseName + SUMMARY_METRICS_EXTENSION));
        assertFilesEqual(outputDetailsFile, new File(TEST_DATA_PATH, expectedOutputFileBaseName + DETAILED_METRICS_EXTENSION));
        assertFilesEqual(outputVariantFile, new File(TEST_DATA_PATH, expectedOutputFileBaseName + VARIANT_METRICS_EXTENSION));
    }

    /**
     * Checks that the two files are the same length, and have the same content, otherwise throws a runtime exception.
     */
    public static void assertFilesEqual(final File f1, final File f2) {
        try {
            final BufferedReader reader1 = IOUtil.openFileForBufferedReading(f1);
            final BufferedReader reader2 = IOUtil.openFileForBufferedReading(f2);

            boolean equal;

            while (true) { // Continue while there are equal lines
                final String line1 = reader1.readLine();
                final String line2 = reader2.readLine();

                if (line1 == null) {     // End of file 1
                    equal = (line2 == null);
                    break;
                }
                else if (line2 == null) {   // file 2 ended before file 1
                    equal = false;
                    break;
                }
                else {
                    if (!line1.startsWith("#") || (!line2.startsWith("#"))) {
                        equal = line1.equals(line2);
                        if (!equal) break;
                    }
                }
            }
            CloserUtil.close(reader1);
            CloserUtil.close(reader2);
            if (!equal) {
                throw new PicardException("Files " + f1.getAbsolutePath() + " and " + f2.getAbsolutePath() + " differ");
            }
        } catch (final IOException e) {
            throw new PicardException("Exception comparing files " + f1 + " and " + f2, e);
        }

    }

    @Test
    public void testGenotypeConcordanceDetails() throws Exception {
        final File outputBaseFileName = new File(OUTPUT_DATA_PATH, "actualGtConc");
        final File outputSummaryFile = new File(outputBaseFileName.getAbsolutePath() + SUMMARY_METRICS_EXTENSION);
        final File outputDetailsFile = new File(outputBaseFileName.getAbsolutePath() + DETAILED_METRICS_EXTENSION);
        outputSummaryFile.deleteOnExit();
        outputDetailsFile.deleteOnExit();

        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.CALL_SAMPLE = "NA12878";
        genotypeConcordance.OUTPUT = outputBaseFileName;

        int returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        final Map<ConcordanceState, Integer> nonZeroCounts = new HashMap<ConcordanceState, Integer>();
        nonZeroCounts.put(new ConcordanceState(VariantCallState.Het, VariantCallState.Het, true), 104);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.HomVar, VariantCallState.HomVar, true), 59);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.FilteredVariant, VariantCallState.FilteredVariant, true), 40);

        ConcordanceResults concordanceResults = genotypeConcordance.getSnpCounter();
        for (final VariantCallState state1 : VariantCallState.values()) {
            for (final VariantCallState state2 : VariantCallState.values()) {
                Integer expectedCount = nonZeroCounts.get(new ConcordanceState(state1, state2, true));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(concordanceResults.getCount(state1, state2, true), expectedCount.intValue());
            }
        }

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceResults.hetSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.hetSpecificity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.hetPpv()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSpecificity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarPpv()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varSpecificity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varPpv()), "1");

        // Now run it again with different samples
        genotypeConcordance.TRUTH_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.CALL_SAMPLE = "NA12891";
        returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        nonZeroCounts.clear();
        nonZeroCounts.put(new ConcordanceState(VariantCallState.HomRef, VariantCallState.Het, true), 31);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.Het, VariantCallState.HomRef, true), 30);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.Het, VariantCallState.Het, true), 50);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.Het, VariantCallState.HomVar, true), 24);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.HomVar, VariantCallState.Het, true), 18);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.HomVar, VariantCallState.HomVar, true), 41);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.FilteredVariant, VariantCallState.FilteredVariant, true), 49);

        concordanceResults = genotypeConcordance.getSnpCounter();
        for (final VariantCallState state1 : VariantCallState.values()) {
            for (final VariantCallState state2 : VariantCallState.values()) {
                Integer expectedCount = nonZeroCounts.get(new ConcordanceState(state1, state2, true));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(concordanceResults.getCount(state1, state2, true), expectedCount.intValue());
            }
        }

        Assert.assertEquals(fmt.format(concordanceResults.hetSensitivity()), "0.480769");
        Assert.assertEquals(fmt.format(concordanceResults.hetSpecificity()), "0.647482");
        Assert.assertEquals(fmt.format(concordanceResults.hetPpv()), "0.505051");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSensitivity()), "0.694915");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSpecificity()), "0.869565");
        Assert.assertEquals(fmt.format(concordanceResults.homVarPpv()), "0.630769");
        Assert.assertEquals(fmt.format(concordanceResults.varSensitivity()), "0.558282");
        Assert.assertEquals(fmt.format(concordanceResults.varSpecificity()), "0.6125");
        Assert.assertEquals(fmt.format(concordanceResults.varPpv()), "0.810976");
    }

    @Test
    public void testGenotypeConcordanceDetailsWithIntervals() throws Exception {
        final File outputBaseFileName = new File(OUTPUT_DATA_PATH, "actualGtConc");
        final File outputSummaryFile = new File(outputBaseFileName.getAbsolutePath() + SUMMARY_METRICS_EXTENSION);
        final File outputDetailsFile = new File(outputBaseFileName.getAbsolutePath() + DETAILED_METRICS_EXTENSION);
        outputSummaryFile.deleteOnExit();
        outputDetailsFile.deleteOnExit();

        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.CALL_SAMPLE = "NA12878";
        genotypeConcordance.INTERVALS = Collections.singletonList(INTERVALS_FILE);
        genotypeConcordance.OUTPUT = outputBaseFileName;

        int returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        final Map<ConcordanceState, Integer> nonZeroCounts = new HashMap<ConcordanceState, Integer>();
        nonZeroCounts.put(new ConcordanceState(VariantCallState.Het, VariantCallState.Het, true), 1);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.FilteredVariant, VariantCallState.FilteredVariant, true), 2);

        ConcordanceResults concordanceResults = genotypeConcordance.getSnpCounter();
        for (final VariantCallState state1 : VariantCallState.values()) {
            for (final VariantCallState state2 : VariantCallState.values()) {
                Integer expectedCount = nonZeroCounts.get(new ConcordanceState(state1, state2, true));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(concordanceResults.getCount(state1, state2, true), expectedCount.intValue());
            }
        }

        final FormatUtil fmt = new FormatUtil();

        Assert.assertEquals(fmt.format(concordanceResults.hetSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.hetSpecificity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.hetPpv()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSensitivity()), "?");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSpecificity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarPpv()), "?");
        Assert.assertEquals(fmt.format(concordanceResults.varSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varSpecificity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varPpv()), "1");

        // Now run it again with different samples
        genotypeConcordance.TRUTH_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.TRUTH_SAMPLE = "NA12878";
        genotypeConcordance.CALL_VCF = CEU_TRIOS_SNPS_VCF;
        genotypeConcordance.CALL_SAMPLE = "NA12891";
        genotypeConcordance.INTERVALS = Collections.singletonList(INTERVALS_FILE);
        returnCode = genotypeConcordance.instanceMain(new String[0]);
        Assert.assertEquals(returnCode, 0);

        nonZeroCounts.clear();
        nonZeroCounts.put(new ConcordanceState(VariantCallState.HomRef, VariantCallState.Het, true), 1);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.Het, VariantCallState.Het, true), 1);
        nonZeroCounts.put(new ConcordanceState(VariantCallState.FilteredVariant, VariantCallState.FilteredVariant, true), 2);

        concordanceResults = genotypeConcordance.getSnpCounter();
        for (final VariantCallState state1 : VariantCallState.values()) {
            for (final VariantCallState state2 : VariantCallState.values()) {
                Integer expectedCount = nonZeroCounts.get(new ConcordanceState(state1, state2, true));
                if (expectedCount == null) expectedCount = 0;
                Assert.assertEquals(concordanceResults.getCount(state1, state2, true), expectedCount.intValue());
            }
        }

        Assert.assertEquals(fmt.format(concordanceResults.hetSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.hetSpecificity()), "0.666667");
        Assert.assertEquals(fmt.format(concordanceResults.hetPpv()), "0.5");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSensitivity()), "?");
        Assert.assertEquals(fmt.format(concordanceResults.homVarSpecificity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.homVarPpv()), "?");
        Assert.assertEquals(fmt.format(concordanceResults.varSensitivity()), "1");
        Assert.assertEquals(fmt.format(concordanceResults.varSpecificity()), "0.666667");
        Assert.assertEquals(fmt.format(concordanceResults.varPpv()), "0.5");
    }

    @Test
    public void testGenotypeConcordanceDoAltAllelesAgree() {
        // A [ref] / T at 10
        final String snpLoc = "chr1";
        final int snpLocStart = 10;
        final int snpLocStop = 10;

        final Allele Aref = Allele.create("A", true);
        final Allele T = Allele.create("T");
        final Allele C = Allele.create("C");
        final Allele Tref = Allele.create("T", true);
        final List<Allele> alleles = Arrays.asList(Aref, T, C);

        final String sampleName = "FOO";
        final GenotypeConcordance genotypeConcordance = new GenotypeConcordance();
        genotypeConcordance.TRUTH_SAMPLE = sampleName;          // Need for lookup of genotypes in doAltAllelesAgree
        genotypeConcordance.CALL_SAMPLE  = sampleName;          // Need for lookup of genotypes in doAltAllelesAgree

        final Genotype gtHomRef = GenotypeBuilder.create(sampleName, Arrays.asList(Aref, Aref));      // Homref (A*A*)
        final VariantContext vcHomRef = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(gtHomRef).make();

        final Genotype gtHetAT = GenotypeBuilder.create(sampleName, Arrays.asList(Aref, T));         // Het (A*T)
        final VariantContext vcHetAT = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(gtHetAT).make();

        final Genotype gtHetCT = GenotypeBuilder.create(sampleName, Arrays.asList(C, T));            // Het (CT)
        final VariantContext vcHetCT = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(gtHetCT).make();

        final Genotype gtHomVarCC = GenotypeBuilder.create(sampleName, Arrays.asList(C, C));            // HomVar (CC)
        final VariantContext vcHomVarCC = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(gtHomVarCC).make();

        final Genotype gtHomVarTT = GenotypeBuilder.create(sampleName, Arrays.asList(T, T));            // HomVar (TT)
        final VariantContext vcHomVarTT = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, alleles).genotypes(gtHomVarTT).make();

        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHomRef, VariantCallState.HomRef, vcHomRef, VariantCallState.HomRef), true);
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHomRef, VariantCallState.HomRef, vcHetAT, VariantCallState.Het), true);
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHetCT, VariantCallState.Het, vcHomRef, VariantCallState.HomRef), true);

        // Check for het vs het.
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHetAT, VariantCallState.Het, vcHetAT, VariantCallState.Het), true);
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHetAT, VariantCallState.Het, vcHetCT, VariantCallState.Het), false);
        // Parameter order should not matter, but we shall check anyway...
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHetCT, VariantCallState.Het, vcHetAT, VariantCallState.Het), false);

        // HomVar vs HomVar
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHomVarCC, VariantCallState.HomVar, vcHomVarCC, VariantCallState.HomVar), true);
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHomVarTT, VariantCallState.HomVar, vcHomVarTT, VariantCallState.HomVar), true);
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHomVarCC, VariantCallState.HomVar, vcHomVarTT, VariantCallState.HomVar), false);
        // Parameter order should not matter, but we shall check anyway...
        Assert.assertEquals(genotypeConcordance.doAllelesAgree(vcHomVarTT, VariantCallState.HomVar, vcHomVarCC, VariantCallState.HomVar), false);

    }
}
