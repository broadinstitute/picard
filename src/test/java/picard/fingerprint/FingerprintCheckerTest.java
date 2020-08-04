package picard.fingerprint;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.stream.Collectors;

import static picard.util.TestNGUtil.compareDoubleWithAccuracy;

/**
 * Created by farjoun on 8/27/15.
 */
public class FingerprintCheckerTest {

    private final double maf = 0.4;
    private final Snp snp = new Snp("test", "chr1", 1, (byte) 'A', (byte) 'C', maf, Collections.singletonList("dummy"));
    private final HaplotypeBlock hb = new HaplotypeBlock(maf);

    private static final double DELTA = 1e-6;

    private static final File TEST_DATA_DIR = new File("testdata/picard/fingerprint/");
    private static final File SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");

    @BeforeClass
    public void setup() {
        hb.addSnp(snp);
    }

    @DataProvider(name = "pLoH")
    public Iterator<Object[]> pLohData() {
        final List<Object[]> listOfDoubles = new ArrayList<>();

        for (int i = 1; i < 20; i++) {
            listOfDoubles.add(new Object[]{i / 40D});
        }
        return listOfDoubles.iterator();
    }

    @Test(dataProvider = "pLoH")
    public void testMatchResults(final double pLoH) {

        final Fingerprint fpObserved = new Fingerprint("test", null, "noop");
        final Fingerprint fpExpected = new Fingerprint("test", null, "noop");

        final HaplotypeProbabilities hpHet = new HaplotypeProbabilitiesFromGenotype(snp, hb, 0.0001, 1.0, 0.0001);
        final HaplotypeProbabilities hpHomRef = new HaplotypeProbabilitiesFromGenotype(snp, hb, 1.0, 0.00001, 0.000000001);

        // Expected is a het
        fpExpected.add(hpHet);

        // Observed is a hom, so possible scenario is that observed is tumor, and expected is normal
        fpObserved.add(hpHomRef);

        // get match results using pLOD
        final MatchResults mr = FingerprintChecker.calculateMatchResults(fpObserved, fpExpected, 0.01, pLoH);

        // make sure that it's more likely to be the same sample, if the observed is "tumor" and the expected is "normal"
        Assert.assertTrue(mr.getLodTN() > mr.getLOD());

        // make sure that the regular LOD is negative (we're comparing a HET to a HOM)
        Assert.assertTrue(mr.getLOD() < 0);

        // make sure that it's more likely to be tumor/normal rather than normal/tumor
        // (a hom normal isn't expected to be measured as a het in the tumor)
        Assert.assertTrue(mr.getLodTN() > mr.getLodNT());
    }

    @DataProvider(name = "checkFingerprintsVcfDataProvider")
    public Object[][] testCheckFingerprintsVcfDataProvider() {
        return new Object[][]{
                {new File(TEST_DATA_DIR, "NA12891.vcf"), new File(TEST_DATA_DIR, "NA12891.fp.vcf"), "NA12891", "NA12891", -1.048021, -2.053484, 1.005462},
                {new File(TEST_DATA_DIR, "NA12891.vcf"), new File(TEST_DATA_DIR, "NA12891.g.vcf"), "NA12891", "NA12891", -1.034969, -2.048616, 1.013647},
                {new File(TEST_DATA_DIR, "NA12892.vcf"), new File(TEST_DATA_DIR, "NA12892.fp.vcf"), "NA12892", "NA12892", -1.105025, -2.166160, 1.061135},
                {new File(TEST_DATA_DIR, "NA12892.vcf"), new File(TEST_DATA_DIR, "NA12892.g.vcf"), "NA12892", "NA12892", -1.091976, -2.161961, 1.069985},
                {new File(TEST_DATA_DIR, "NA12891.vcf"), new File(TEST_DATA_DIR, "NA12892.fp.vcf"), "NA12891", "NA12892", -7.024770, -2.109822, -4.914948},
                {new File(TEST_DATA_DIR, "NA12891.vcf"), new File(TEST_DATA_DIR, "NA12892.g.vcf"), "NA12891", "NA12892", -7.981971, -2.105623, -5.876347},
                {new File(TEST_DATA_DIR, "NA12892.vcf"), new File(TEST_DATA_DIR, "NA12891.fp.vcf"), "NA12892", "NA12891", -7.024770, -2.109822, -4.914948},
                {new File(TEST_DATA_DIR, "NA12892.vcf"), new File(TEST_DATA_DIR, "NA12891.g.vcf"), "NA12892", "NA12891", -7.924964, -2.104955, -5.820009},
                {new File(TEST_DATA_DIR, "emptyNA12892.vcf"), new File(TEST_DATA_DIR, "NA12891.g.vcf"), "NA12892", "NA12891", 0, 0, 0},

        };
    }

    @Test(dataProvider = "checkFingerprintsVcfDataProvider")
    public void testCheckFingerprintsVcf(final File vcfFile, final File genotypesFile, final String observedSampleAlias, final String expectedSampleAlias,
                                         final double llExpectedSample, final double llRandomSample, final double lodExpectedSample) throws IOException {
        final Path indexedInputVcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(vcfFile, "fingerprintcheckertest.tmp.").toPath();
        final Path indexedGenotypesVcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(genotypesFile, "fingerprintcheckertest.tmp.").toPath();

        final FingerprintChecker fpChecker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final List<FingerprintResults> results = fpChecker.checkFingerprintsFromPaths(Collections.singletonList(indexedInputVcf),
                Collections.singletonList(indexedGenotypesVcf),
                observedSampleAlias,
                expectedSampleAlias);
        Assert.assertEquals(results.size(), 1);
        final FingerprintResults fpr = results.get(0);
        Assert.assertNull(fpr.getReadGroup());
        Assert.assertEquals(fpr.getSampleAlias(), observedSampleAlias);
        final MatchResults mr = fpr.getMatchResults().first();
        Assert.assertEquals(mr.getSample(), expectedSampleAlias);
        Assert.assertEquals(mr.getSampleLikelihood(), llExpectedSample, DELTA);
        Assert.assertEquals(mr.getPopulationLikelihood(), llRandomSample, DELTA);
        Assert.assertEquals(mr.getLOD(), lodExpectedSample, DELTA);
    }

    @Test(dataProvider = "checkFingerprintsVcfDataProvider")
    public void testFingerprintVcf(final File vcfFile, final File genotypesFile, final String observedSampleAlias, final String expectedSampleAlias,
                                   final double llExpectedSample, final double llRandomSample, final double lodExpectedSample) {
        final FingerprintChecker fpChecker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final Map<FingerprintIdDetails, Fingerprint> fp1 = fpChecker.fingerprintVcf(vcfFile.toPath());

        Assert.assertFalse(fp1.isEmpty());
    }

    @Test(dataProvider = "checkFingerprintsVcfDataProvider")
    public void testFingerprintSwapEqual(final File vcfFile, final File genotypesFile, final String observedSampleAlias, final String expectedSampleAlias,
                                         final double llExpectedSample, final double llRandomSample, final double lodExpectedSample) {
        final FingerprintChecker fpChecker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final Map<FingerprintIdDetails, Fingerprint> fp1Map = fpChecker.fingerprintVcf(vcfFile.toPath());
        final Map<FingerprintIdDetails, Fingerprint> fp2Map = fpChecker.fingerprintVcf(genotypesFile.toPath());

        for (Fingerprint fp1 : fp1Map.values()) {
            for (Fingerprint fp2 : fp2Map.values()) {
                final MatchResults matchResults12 = FingerprintChecker.calculateMatchResults(fp1, fp2);
                final MatchResults matchResults21 = FingerprintChecker.calculateMatchResults(fp2, fp1);
                compareDoubleWithAccuracy(matchResults12.getLOD(), matchResults21.getLOD(), 1e-10);
                compareDoubleWithAccuracy(matchResults12.getLodTN(), matchResults21.getLodNT(), 1e-10);
                compareDoubleWithAccuracy(matchResults12.getLodNT(), matchResults21.getLodTN(), 1e-10);
            }
        }
    }

    @DataProvider
    Object[][] samFilesforFingerprinting() {
        return new Object[][]{
                {"NA12891.over.fingerprints.r1.sam", true},
                {"aligned_queryname_sorted.sam", false},
                {"aligned_unsorted.sam", false},
        };
    }

    @Test(dataProvider = "samFilesforFingerprinting")
    public void testTerminateOnBadFile(final String fileName, final boolean shouldSucceed) throws Throwable {
        final FingerprintChecker fpChecker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final File samFile = new File(TEST_DATA_DIR, fileName);

        final Assert.ThrowingRunnable runnable = () -> fpChecker.fingerprintFiles(Collections.singletonList(samFile.toPath()), 1, 2, TimeUnit.MINUTES);
        if (shouldSucceed) {
            runnable.run();
        } else {
            Assert.assertThrows(PicardException.class, runnable);
        }
    }

    @DataProvider(name = "checkFingerprintsSamDataProvider")
    public Object[][] testCheckFingerprintsSamDataProvider() {
        final File na12891_r1 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
        final File na12891_r2 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r2.sam");
        final File na12892_r1 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
        final File na12892_r2 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r2.sam");

        final File na12891_noRg = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.noRgTag.sam");

        return new Object[][]{
                {na12891_r1, na12891_r2, true, true},
                {na12892_r1, na12892_r2, true, true},
                {na12892_r1, na12891_r2, false, true},
                {na12891_r1, na12891_noRg, true, true},
                {na12892_r1, na12891_noRg, false, true},

                {na12891_r2, na12891_r2, true, false},
                {na12892_r2, na12892_r2, true, false},
                {na12892_r2, na12891_r2, false, false},
                {na12891_r2, na12891_noRg, true, false},
                {na12892_r2, na12891_noRg, false, false},
        };
    }

    @Test(dataProvider = "checkFingerprintsSamDataProvider")
    public void testCheckFingerprintsSam(final File samFile1, final File samFile2, final boolean expectedMatch, final boolean silent) throws IOException {

        final File metricsFile = File.createTempFile("crosscheck", ".crosscheck_metrics");
        metricsFile.deleteOnExit();

        final String[] args = {
                "EXPECT_ALL_GROUPS_TO_MATCH=true",
                "LOD_THRESHOLD=-1",
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING.getAbsolutePath(),
                "I=" + samFile1.getAbsolutePath(),
                "I=" + samFile2.getAbsolutePath(),
                "VALIDATION_STRINGENCY=" + (silent ? "SILENT" : "LENIENT"),
                "CROSSCHECK_BY=FILE",
                "OUTPUT=" + metricsFile.getAbsolutePath()
        };

        Assert.assertEquals(new CrosscheckFingerprints().instanceMain(args), expectedMatch ? 0 : 1);

        // this part checks that the results are symmetric (i.e LOD x,y == LOD y,x)
        final MetricsFile<CrosscheckMetric, Double> metricsFileReader = new MetricsFile<>();
        metricsFileReader.read(new FileReader(metricsFile));
        final List<CrosscheckMetric> metrics = metricsFileReader.getMetrics();

        final Map<Set<String>, Set<String>> collectedLOD = metrics.stream()
                // we sometimes get -0.0, and so the comparison fails...hence the 'String.valueOf(s.LOD_SCORE + 0)'
                .collect(Collectors.groupingBy(s -> CollectionUtil.makeSet(s.LEFT_GROUP_VALUE, s.RIGHT_GROUP_VALUE), Collectors.mapping(s -> String.valueOf(s.LOD_SCORE + 0), Collectors.toSet())));

        for (Map.Entry<Set<String>, Set<String>> entry : collectedLOD.entrySet()) {
            if (entry.getValue().size() > 1) {

                final List<CrosscheckMetric> mismatchingMetrics = metrics.stream()
                        .filter(s -> CollectionUtil.makeSet(s.LEFT_GROUP_VALUE, s.RIGHT_GROUP_VALUE).equals(entry.getKey())).collect(Collectors.toList());

                Assert.fail("Metrics disagree: LOD scores are: \n[" + String.join(", ", entry.getValue()) +
                        "],\n from the following metrics: \n" + mismatchingMetrics.get(0) + mismatchingMetrics.get(1));
            }
        }

        final Map<Set<String>, Set<String>> collectedTN_LOD = metrics.stream()
                // we sometimes get -0.0, and so the comparison fails...hence the 'String.valueOf(s.LOD_SCORE + 0)'
                .collect(Collectors.groupingBy(s -> CollectionUtil.makeSet(s.LEFT_GROUP_VALUE, s.RIGHT_GROUP_VALUE), Collectors.mapping(s -> String.valueOf(s.LOD_SCORE_TUMOR_NORMAL + 0), Collectors.toSet())));

        final Map<Set<String>, Set<String>> collectedNT_LOD = metrics.stream()
                // we sometimes get -0.0, and so the comparison fails...hence the 'String.valueOf(s.LOD_SCORE + 0)'
                .collect(Collectors.groupingBy(s -> CollectionUtil.makeSet(s.LEFT_GROUP_VALUE, s.RIGHT_GROUP_VALUE), Collectors.mapping(s -> String.valueOf(s.LOD_SCORE_NORMAL_TUMOR + 0), Collectors.toSet())));

        for (Set<String> key : collectedLOD.keySet()) {
            if (key.size() == 1) {
                continue;
            }
            Assert.assertEquals(key.size(), 2, "Odd..should be size 2 here.");

            final Set<String> TN_values = collectedTN_LOD.get(key);
            final Set<String> NT_values = collectedNT_LOD.get(key);

            final Set<String> allValues = new HashSet<>();
            allValues.addAll(TN_values);
            allValues.addAll(NT_values);

            if (allValues.size() > 2) {

                final List<CrosscheckMetric> mismatchingMetrics = metrics.stream()
                        .filter(s -> CollectionUtil.makeSet(s.LEFT_GROUP_VALUE, s.RIGHT_GROUP_VALUE).equals(key)).collect(Collectors.toList());

                Assert.fail(
                        "Metrics disagree: TN_LOD scores are: [" + String.join(", ", TN_values) + "]\n" +
                                "Metrics disagree: TN_LOD scores are: [" + String.join(", ", NT_values) + "]\n" +
                                "from the following metrics: \n" + mismatchingMetrics.get(0) + "\n" + mismatchingMetrics.get(1));
            }

        }
    }

    @DataProvider(name = "checkFingerprintsSamDataProviderFail")
    public Object[][] testCheckFingerprintsSamDataProviderFail() {
        final File na12891_r1 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
        final File na12892_r1 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
        final File na12891_noRg = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.noRgTag.sam");

        return new Object[][]{
                {na12892_r1, na12891_noRg, false},
                {na12891_r1, na12891_noRg, true},
        };
    }

    @Test(dataProvider = "checkFingerprintsSamDataProviderFail", expectedExceptions = PicardException.class)
    public void testCheckFingerprintsFail(final File samFile1, final File samFile2, final boolean expectedMatch) {
        final String[] args = {
                "EXPECT_ALL_GROUPS_TO_MATCH=true",
                "LOD_THRESHOLD=-1",
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING.getAbsolutePath(),
                "I=" + samFile1.getAbsolutePath(),
                "I=" + samFile2.getAbsolutePath(),
                "VALIDATION_STRINGENCY=STRICT",
                "CROSSCHECK_BY=FILE",
        };

        new CrosscheckFingerprints().instanceMain(args);
    }

    @DataProvider(name = "queryableData")
    public Iterator<Object[]> queryableData() throws IOException {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{new File(TEST_DATA_DIR, "NA12891.fp.vcf"), false});
        tests.add(new Object[]{new File(TEST_DATA_DIR, "NA12891.vcf"), false});
        tests.add(new Object[]{VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.vcf"), "fingerprintcheckertest.tmp."), true});
        tests.add(new Object[]{VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.vcf.gz"), "fingerprintcheckertest.tmp."), true});

        return tests.iterator();
    }

    @Test(dataProvider = "queryableData")
    public void testQueryable(final File vcf, boolean expectedQueryable) {
        try (VCFFileReader reader = new VCFFileReader(vcf, false)) {
            Assert.assertEquals(reader.isQueryable(), expectedQueryable);
        }
    }

    @DataProvider()
    Object[][] mergeIsDafeProvider() {
        final HaplotypeProbabilitiesFromSequence hp1 = new HaplotypeProbabilitiesFromSequence(hb);
        final HaplotypeProbabilitiesFromSequence hp2 = new HaplotypeProbabilitiesFromSequence(hb);

        addObservation(hp1, hb, 5, hb.getFirstSnp().getAllele1());
        addObservation(hp1, hb, 1, (byte) (hb.getFirstSnp().getAllele1() + 1));
        addObservation(hp2, hb, 3, hb.getFirstSnp().getAllele1());
        addObservation(hp2, hb, 2, hb.getFirstSnp().getAllele2());

        final HaplotypeProbabilitiesFromContaminatorSequence hpcs1 = new HaplotypeProbabilitiesFromContaminatorSequence(hb, .1);
        final HaplotypeProbabilitiesFromContaminatorSequence hpcs2 = new HaplotypeProbabilitiesFromContaminatorSequence(hb, .1);

        addObservation(hpcs1, hb, 5, hb.getFirstSnp().getAllele1());
        addObservation(hpcs1, hb, 1, (byte) (hb.getFirstSnp().getAllele1() + 1));
        addObservation(hpcs2, hb, 3, hb.getFirstSnp().getAllele1());
        addObservation(hpcs2, hb, 1, hb.getFirstSnp().getAllele1());

        return new Object[][]{
                {new HaplotypeProbabilitiesFromGenotype(snp, hb, 5D, 0D, 10D), new HaplotypeProbabilitiesFromGenotype(snp, hb, 0D, 10D, 100D)},
                {
                        new HaplotypeProbabilityOfNormalGivenTumor(
                                new HaplotypeProbabilitiesFromGenotype(snp, hb, 5D, 0D, 10D), .05),
                        new HaplotypeProbabilityOfNormalGivenTumor(
                                new HaplotypeProbabilitiesFromGenotype(snp, hb, 0D, 10D, 100D), 0.05)},
                {hp1, hp2},
                {hpcs1, hpcs2},
        };
    }

    private static void addObservation(final HaplotypeProbabilitiesFromSequence haplotypeProb,
                                       final HaplotypeBlock haplotypeBlock, final int count, final byte allele) {
        for (int i = 0; i < count; i++) {
            haplotypeProb.addToProbs(haplotypeBlock.getFirstSnp(), allele, (byte) 30);
        }
    }

    @Test(dataProvider = "mergeIsDafeProvider")
    public void testMergeHaplotypeProbabilitiesIsSafe(final HaplotypeProbabilities hp1,
                                                      final HaplotypeProbabilities hp2) {

        final HaplotypeProbabilities merged1 = hp1.deepCopy().merge(hp2);
        final HaplotypeProbabilities merged2 = hp1.deepCopy().merge(hp2);

        Assert.assertEquals(merged1.getLikelihoods(), merged2.getLikelihoods());
    }

    @Test(dataProvider = "mergeIsDafeProvider")
    public void testMergeFingerprintIsSafe(final HaplotypeProbabilities hp1, final HaplotypeProbabilities hp2) {

        final Fingerprint fpA = new Fingerprint("test2", null, "none");
        final Fingerprint fpB = new Fingerprint("test2", null, "none");

        final Fingerprint fp1 = new Fingerprint("test1", null, "none");
        fp1.add(hp1);

        final Fingerprint fp2 = new Fingerprint("test1", null, "none");
        fp2.add(hp2);

        fpA.merge(fp1);
        fpB.merge(fp1);

        Assert.assertNotEquals(fpA.keySet().size(), 0);
        for (HaplotypeBlock hb : fpA.keySet()) {
            Assert.assertEquals(fpA.get(hb).getLikelihoods(), fpB.get(hb).getLikelihoods());
        }

        fpA.merge(fp2);
        fpB.merge(fp2);
        fpB.merge(fp2);

        Assert.assertNotEquals(fpA.keySet().size(), 0);
        for (HaplotypeBlock hb : fpA.keySet()) {
            Assert.assertNotEquals(fpA.get(hb), fpB.get(hb));
        }

        fpA.merge(fp2);

        Assert.assertNotEquals(fpA.keySet().size(), 0);
        for (HaplotypeBlock hb : fpA.keySet()) {
            Assert.assertEquals(fpA.get(hb).getLikelihoods(), fpB.get(hb).getLikelihoods());
        }
    }

    @Test
    public void testMergeIsSafeFromSequence() {
        final Path na12891_r1 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam").toPath();
        final Path na12891_r2 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r2.sam").toPath();
        final Path na12892_r1 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam").toPath();
        final Path na12892_r2 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r2.sam").toPath();

        final List<Path> listOfFiles = Arrays.asList(na12891_r1, na12891_r2, na12892_r1, na12892_r2);
        final FingerprintChecker checker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final Map<FingerprintIdDetails, Fingerprint> fingerprintIdDetailsFingerprintMap = checker.fingerprintFiles(listOfFiles, 1, 0, TimeUnit.DAYS);

        final Fingerprint combinedFp = new Fingerprint("test", null, null);
        fingerprintIdDetailsFingerprintMap.values().forEach(combinedFp::merge);

        final Fingerprint combinedFp2 = new Fingerprint("test2", null, null);
        fingerprintIdDetailsFingerprintMap.values().forEach(combinedFp2::merge);

        Assert.assertNotEquals(combinedFp.keySet().size(), 0);
        FingerprintingTestUtils.assertFingerPrintHPsAreEqual(combinedFp, combinedFp2);

        final Function<FingerprintIdDetails, String> bySample = Fingerprint.getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType.SAMPLE);
        final Map<FingerprintIdDetails, Fingerprint> fingerprintIdDetailsFingerprintMap1 =
                Fingerprint.mergeFingerprintsBy(fingerprintIdDetailsFingerprintMap, bySample);

        final Map<FingerprintIdDetails, Fingerprint> fingerprintIdDetailsFingerprintMap2 =
                Fingerprint.mergeFingerprintsBy(fingerprintIdDetailsFingerprintMap, bySample);

        Assert.assertNotEquals(fingerprintIdDetailsFingerprintMap1.keySet().size(), 0);

        for (final FingerprintIdDetails fpd1 : fingerprintIdDetailsFingerprintMap1.keySet()) {
            final Fingerprint fingerprint1 = fingerprintIdDetailsFingerprintMap1.get(fpd1);
            final Fingerprint fingerprint2 = fingerprintIdDetailsFingerprintMap2.get(fpd1);

            Assert.assertNotEquals(fingerprint1.keySet().size(), 0);
            FingerprintingTestUtils.assertFingerPrintHPsAreEqual(fingerprint1, fingerprint2);
        }
    }

    @Test
    public void testMergeIsSafeFromVCF() {
        final Path na12891_fp = TEST_DATA_DIR.toPath().resolve("NA12891.fp.vcf");
        final Path na12891_g = TEST_DATA_DIR.toPath().resolve("NA12891.vcf");
        final Path na12892_fp = TEST_DATA_DIR.toPath().resolve("NA12892.fp.vcf");
        final Path na12892_g = TEST_DATA_DIR.toPath().resolve("NA12892.vcf");

        final List<Path> listOfFiles = Arrays.asList(na12891_fp, na12891_g, na12892_fp, na12892_g);
        final FingerprintChecker checker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);
        final Map<FingerprintIdDetails, Fingerprint> fingerprintIdDetailsFingerprintMap = checker.fingerprintFiles(listOfFiles, 1, 0, TimeUnit.DAYS);

        final Fingerprint combinedFp = new Fingerprint("test", null, null);
        fingerprintIdDetailsFingerprintMap.values().forEach(combinedFp::merge);

        final Fingerprint combinedFp2 = new Fingerprint("test", null, null);
        fingerprintIdDetailsFingerprintMap.values().forEach(combinedFp2::merge);

        Assert.assertNotEquals(combinedFp.keySet().size(), 0);
        FingerprintingTestUtils.assertFingerPrintHPsAreEqual(combinedFp, combinedFp2);

        final Function<FingerprintIdDetails, String> bySample = Fingerprint.getFingerprintIdDetailsStringFunction(CrosscheckMetric.DataType.SAMPLE);

        final Map<FingerprintIdDetails, Fingerprint> fingerprintIdDetailsFingerprintMap1 =
                Fingerprint.mergeFingerprintsBy(fingerprintIdDetailsFingerprintMap, bySample);
        final Map<FingerprintIdDetails, Fingerprint> fingerprintIdDetailsFingerprintMap2 =
                Fingerprint.mergeFingerprintsBy(fingerprintIdDetailsFingerprintMap, bySample);

        Assert.assertNotEquals(fingerprintIdDetailsFingerprintMap1.keySet().size(), 0);
        Assert.assertEquals(fingerprintIdDetailsFingerprintMap1.keySet().size(), fingerprintIdDetailsFingerprintMap2.size());

        for (final FingerprintIdDetails fpd1 : fingerprintIdDetailsFingerprintMap1.keySet()) {
            final Fingerprint fingerprint1 = fingerprintIdDetailsFingerprintMap1.get(fpd1);
            final Fingerprint fingerprint2 = fingerprintIdDetailsFingerprintMap2.get(fpd1);

            Assert.assertNotEquals(fingerprint1.keySet().size(), 0);
            FingerprintingTestUtils.assertFingerPrintHPsAreEqual(fingerprint1, fingerprint2);
        }
    }

    @Test
    public void testWriteFingerprint() throws IOException {
        final File haplotype_db = new File(TEST_DATA_DIR, "haplotypeMap_small.vcf");
        final File vcfInput = new File(TEST_DATA_DIR, "testSample_small.vcf");
        final File fasta = new File(TEST_DATA_DIR, "reference.fasta");
        final File vcfExpected = new File(TEST_DATA_DIR, "expectedFingerprint_small.vcf");
        final FingerprintChecker fpchecker = new FingerprintChecker(haplotype_db);
        final Fingerprint fp = fpchecker.fingerprintVcf(vcfInput.toPath()).values().iterator().next();

        final File vcfOutput = File.createTempFile("fingerprint", ".vcf");
        FingerprintUtils.writeFingerPrint(fp, vcfOutput, fasta, "Dummy", "Testing");

        VcfTestUtils.assertVcfFilesAreEqual(vcfOutput, vcfExpected);
    }
}