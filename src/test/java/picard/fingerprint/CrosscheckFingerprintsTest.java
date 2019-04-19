package picard.fingerprint;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.vcf.SamTestUtils;
import picard.vcf.VcfTestUtils;

import java.io.*;
import java.lang.reflect.Field;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Tests for CrosscheckFingerprints
 */
public class CrosscheckFingerprintsTest extends CommandLineProgramTest {

    private final File TEST_DATA_DIR = new File("testdata/picard/fingerprint/");
    private final File HAPLOTYPE_MAP = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");
    private final File HAPLOTYPE_MAP_FOR_CRAMS = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.shifted.for.crams.txt");
    private final File HAPLOTYPE_MAP_MOVED = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.moved.txt");

    private final File NA12891_r1_sam = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
    private final File NA12891_r2_sam = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r2.sam");
    private final File NA12891_r1_sam_shifted_for_cram = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.shifted.for.crams.r1.sam");
    private final File NA12891_r2_sam_shifted_for_cram = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.shifted.for.crams.r2.sam");

    //this is a copy of a previous one, but with a different sample name
    private final File NA12891_named_NA12892_r1_sam = new File(TEST_DATA_DIR, "NA12891_named_NA12892.over.fingerprints.r1.sam");

    private final File NA12892_r1_sam = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
    private final File NA12892_r2_sam = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r2.sam");
    private final File NA12892_r1_sam_shifted_for_cram = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.shifted.for.crams.r1.sam");
    private final File NA12892_r2_sam_shifted_for_cram = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.shifted.for.crams.r2.sam");

    private final File NA12891_r1_one_rg_no_fingerprint_sam = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.one.rg.no.fingerprint.sam");

    private File NA12891_r1, NA12891_r2, NA12891_named_NA12892_r1, NA12892_r1, NA12892_r2;
    private File NA12891_r1_cram, NA12891_r2_cram, NA12892_r1_cram, NA12892_r2_cram;
    private File NA12891_r1_shifted_bam, NA12891_r2_shifted_bam, NA12892_r1_shifted_bam, NA12892_r2_shifted_bam;
    private final File referenceForCrams = new File(TEST_DATA_DIR, "reference.shifted.for.crams.fasta");

    private final int NA12891_r1_RGs = 27;
    private final int NA12891_r2_RGs = 26;
    private final int NA12892_r1_RGs = 25;
    private final int NA12892_r2_RGs = 26;

    private static File  NA12891_1_vcf;
    private static File  NA12891_2_vcf;
    private static File  NA12892_1_vcf;
    private static File  NA12891_swapped_nonref_g_vcf;
    private static File  NA12892_2_vcf;
    private static File  NA12891_named_NA12892_vcf;
    private static File  NA12891_g_vcf;
    private static File  NA12892_g_vcf;
    private static File  NA12892_and_NA123891_vcf;
    private static File  NA12892_and_NA123891_part1_vcf;
    private static File  NA12892_and_NA123891_part2_vcf;
    private static File  NA12892_and_NA123891_part3_vcf;
    private static File NA12891_no_fp_sites_vcf;
    private static File NA12891_no_fp_sites_and_NA12892_vcf;

    private static final Map<CrosscheckMetric.DataType, List<String>> lookupMap = new HashMap<>(4);
    
    @BeforeClass
    public void setup() throws IOException {
        NA12891_r1 = SamTestUtils.createIndexedBamOrCram(NA12891_r1_sam, NA12891_r1_sam, SamReader.Type.BAM_TYPE);
        NA12891_r2 = SamTestUtils.createIndexedBamOrCram(NA12891_r2_sam, NA12891_r2_sam, SamReader.Type.BAM_TYPE);
        NA12891_named_NA12892_r1 = SamTestUtils.createIndexedBamOrCram(NA12891_named_NA12892_r1_sam, NA12891_named_NA12892_r1_sam, SamReader.Type.BAM_TYPE);
        NA12892_r1 = SamTestUtils.createIndexedBamOrCram(NA12892_r1_sam, NA12892_r1_sam, SamReader.Type.BAM_TYPE);
        NA12892_r2 = SamTestUtils.createIndexedBamOrCram(NA12892_r2_sam, NA12892_r2_sam, SamReader.Type.BAM_TYPE);

        NA12891_r1_cram = SamTestUtils.createIndexedBamOrCram(NA12891_r1_sam_shifted_for_cram, NA12891_r1_sam_shifted_for_cram, SamReader.Type.CRAM_TYPE, referenceForCrams);
        NA12891_r2_cram = SamTestUtils.createIndexedBamOrCram(NA12891_r2_sam_shifted_for_cram, NA12891_r2_sam_shifted_for_cram, SamReader.Type.CRAM_TYPE, referenceForCrams);
        NA12892_r1_cram = SamTestUtils.createIndexedBamOrCram(NA12892_r1_sam_shifted_for_cram, NA12892_r1_sam_shifted_for_cram, SamReader.Type.CRAM_TYPE, referenceForCrams);
        NA12892_r2_cram = SamTestUtils.createIndexedBamOrCram(NA12892_r2_sam_shifted_for_cram, NA12892_r2_sam_shifted_for_cram, SamReader.Type.CRAM_TYPE, referenceForCrams);

        NA12891_r1_shifted_bam = SamTestUtils.createIndexedBamOrCram(NA12891_r1_sam_shifted_for_cram, NA12891_r1_sam_shifted_for_cram, SamReader.Type.BAM_TYPE);
        NA12891_r2_shifted_bam = SamTestUtils.createIndexedBamOrCram(NA12891_r2_sam_shifted_for_cram, NA12891_r2_sam_shifted_for_cram, SamReader.Type.BAM_TYPE);
        NA12892_r1_shifted_bam = SamTestUtils.createIndexedBamOrCram(NA12892_r1_sam_shifted_for_cram, NA12892_r1_sam_shifted_for_cram, SamReader.Type.BAM_TYPE);
        NA12892_r2_shifted_bam = SamTestUtils.createIndexedBamOrCram(NA12892_r2_sam_shifted_for_cram, NA12892_r2_sam_shifted_for_cram, SamReader.Type.BAM_TYPE);

        lookupMap.put(CrosscheckMetric.DataType.FILE, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.FILE).addAll(Arrays.asList("LEFT_FILE", "RIGHT_FILE"));

        lookupMap.put(CrosscheckMetric.DataType.SAMPLE, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.SAMPLE).addAll(Arrays.asList("LEFT_SAMPLE", "RIGHT_SAMPLE"));
        lookupMap.get(CrosscheckMetric.DataType.SAMPLE).addAll(lookupMap.get(CrosscheckMetric.DataType.FILE));

        lookupMap.put(CrosscheckMetric.DataType.LIBRARY, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.LIBRARY).addAll(Arrays.asList("LEFT_LIBRARY", "RIGHT_LIBRARY"));
        lookupMap.get(CrosscheckMetric.DataType.LIBRARY).addAll(lookupMap.get(CrosscheckMetric.DataType.SAMPLE));

        lookupMap.put(CrosscheckMetric.DataType.READGROUP, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.READGROUP).addAll(Arrays.asList("LEFT_RUN_BARCODE", "LEFT_LANE",
                "LEFT_MOLECULAR_BARCODE_SEQUENCE","RIGHT_RUN_BARCODE",
                "RIGHT_LANE", "RIGHT_MOLECULAR_BARCODE_SEQUENCE"));
        lookupMap.get(CrosscheckMetric.DataType.READGROUP).addAll(lookupMap.get(CrosscheckMetric.DataType.LIBRARY));

        NA12891_1_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.vcf"), "fingerprint");
        NA12891_2_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.fp.vcf"), "fingerprint");
        NA12891_named_NA12892_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891_named_NA12892.vcf"), "fingerprint");
        NA12892_1_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12892.vcf"), "fingerprint");
        NA12892_2_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12892.fp.vcf"), "fingerprint");
        NA12891_g_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.g.vcf"), "fingerprint");
        NA12891_swapped_nonref_g_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.with_swapped_NON_REF.g.vcf"), "fingerprint");

        NA12892_g_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12892.g.vcf"), "fingerprint");
        NA12892_and_NA123891_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891andNA12892.vcf"), "fingerprint");
        NA12892_and_NA123891_part1_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891andNA12892_part1.vcf"), "fingerprint");
        NA12892_and_NA123891_part2_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891andNA12892_part2.vcf"), "fingerprint");
        NA12892_and_NA123891_part3_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891andNA12892_part3.vcf"), "fingerprint");
        NA12891_no_fp_sites_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.no.fp.sites.vcf"), "fingerprint");
        NA12891_no_fp_sites_and_NA12892_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.no.fp.sites.and.NA12892.vcf"), "fingerprint");
    }

    @Override
    public String getCommandLineProgramName() {
        return CrosscheckFingerprints.class.getSimpleName();
    }

    @DataProvider(name = "bamFilesRGs")
    public Object[][] bamAndCramFilesRGs() {
        return new Object[][] {
                {NA12891_r1, NA12891_r2, false, 0, (NA12891_r1_RGs + NA12891_r2_RGs) * (NA12891_r1_RGs + NA12891_r2_RGs )},
                {NA12891_r1, NA12892_r1, false, 0, (NA12891_r1_RGs + NA12892_r1_RGs) * (NA12891_r1_RGs + NA12892_r1_RGs )},
                {NA12891_r1, NA12892_r2, false, 0, (NA12891_r1_RGs + NA12892_r2_RGs) * (NA12891_r1_RGs + NA12892_r2_RGs )},
                {NA12892_r1, NA12892_r2, false, 0, (NA12892_r1_RGs + NA12892_r2_RGs) * (NA12892_r1_RGs + NA12892_r2_RGs )},
                {NA12892_r2, NA12891_r2, false, 0, (NA12892_r2_RGs + NA12891_r2_RGs) * (NA12892_r2_RGs + NA12891_r2_RGs )},
                {NA12892_r2, NA12891_r1, false, 0, (NA12892_r2_RGs + NA12891_r1_RGs) * (NA12892_r2_RGs + NA12891_r1_RGs )},
                {NA12891_r1, NA12891_r2, true,  0, (NA12891_r1_RGs + NA12891_r2_RGs) * (NA12891_r1_RGs + NA12891_r2_RGs )},
                {NA12891_r1, NA12892_r1, true,  1, (NA12891_r1_RGs + NA12892_r1_RGs) * (NA12891_r1_RGs + NA12892_r1_RGs )},
                {NA12891_r1, NA12892_r2, true,  1, (NA12891_r1_RGs + NA12892_r2_RGs) * (NA12891_r1_RGs + NA12892_r2_RGs )},
                {NA12892_r1, NA12892_r2, true,  0, (NA12892_r1_RGs + NA12892_r2_RGs) * (NA12892_r1_RGs + NA12892_r2_RGs )},
                {NA12892_r2, NA12891_r2, true,  1, (NA12892_r2_RGs + NA12891_r2_RGs) * (NA12892_r2_RGs + NA12891_r2_RGs )},
                {NA12892_r2, NA12891_r1, true, 1, (NA12892_r2_RGs + NA12891_r1_RGs) * (NA12892_r2_RGs + NA12891_r1_RGs)},
                {NA12891_r1_cram, NA12891_r2_cram, false, 0, (NA12891_r1_RGs + NA12891_r2_RGs) * (NA12891_r1_RGs + NA12891_r2_RGs)},
                {NA12891_r1_cram, NA12892_r1_cram, false, 0, (NA12891_r1_RGs + NA12892_r1_RGs) * (NA12891_r1_RGs + NA12892_r1_RGs)},
                {NA12891_r1_cram, NA12892_r2_cram, false, 0, (NA12891_r1_RGs + NA12892_r2_RGs) * (NA12891_r1_RGs + NA12892_r2_RGs)},
                {NA12892_r1_cram, NA12892_r2_cram, false, 0, (NA12892_r1_RGs + NA12892_r2_RGs) * (NA12892_r1_RGs + NA12892_r2_RGs)},
                {NA12892_r2_cram, NA12891_r2_cram, false, 0, (NA12892_r2_RGs + NA12891_r2_RGs) * (NA12892_r2_RGs + NA12891_r2_RGs)},
                {NA12892_r2_cram, NA12891_r1_cram, false, 0, (NA12892_r2_RGs + NA12891_r1_RGs) * (NA12892_r2_RGs + NA12891_r1_RGs)},
                {NA12891_r1_cram, NA12891_r2_cram, true, 0, (NA12891_r1_RGs + NA12891_r2_RGs) * (NA12891_r1_RGs + NA12891_r2_RGs)},
                {NA12891_r1_cram, NA12892_r1_cram, true, 1, (NA12891_r1_RGs + NA12892_r1_RGs) * (NA12891_r1_RGs + NA12892_r1_RGs)},
                {NA12891_r1_cram, NA12892_r2_cram, true, 1, (NA12891_r1_RGs + NA12892_r2_RGs) * (NA12891_r1_RGs + NA12892_r2_RGs)},
                {NA12892_r1_cram, NA12892_r2_cram, true, 0, (NA12892_r1_RGs + NA12892_r2_RGs) * (NA12892_r1_RGs + NA12892_r2_RGs)},
                {NA12892_r2_cram, NA12891_r2_cram, true, 1, (NA12892_r2_RGs + NA12891_r2_RGs) * (NA12892_r2_RGs + NA12891_r2_RGs)},
                {NA12892_r2_cram, NA12891_r1_cram, true, 1, (NA12892_r2_RGs + NA12891_r1_RGs) * (NA12892_r2_RGs + NA12891_r1_RGs)}
        };
    }

    @Test(dataProvider = "bamFilesRGs")
    public void testCrossCheckRGs(final File file1, final File file2, final boolean expectAllMatch, final int expectedRetVal, final int expectedNMetrics) throws IOException {

        File metrics = File.createTempFile("Fingerprinting", "NA1291.RG.crosscheck_metrics");
        metrics.deleteOnExit();

        final List<String> args = new ArrayList<>(Arrays.asList("INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "LOD_THRESHOLD=" + -2.0,
                "EXPECT_ALL_GROUPS_TO_MATCH=" + expectAllMatch)
        );

        if (file1.getName().endsWith(SamReader.Type.CRAM_TYPE.fileExtension())) {
            args.add("R=" + referenceForCrams);
            args.add("HAPLOTYPE_MAP=" + HAPLOTYPE_MAP_FOR_CRAMS);
        } else {
            args.add("HAPLOTYPE_MAP=" + HAPLOTYPE_MAP);
        }

        doTest(args.toArray(new String[args.size()]), metrics, expectedRetVal, expectedNMetrics, CrosscheckMetric.DataType.READGROUP, expectAllMatch);
    }

    @DataProvider(name = "cramsWithNoReference")
    public Object[][] cramsWithNoReference() {
        return new Object[][]{
                {NA12891_r1_cram, NA12891_r2_cram},
                {NA12891_r1_cram, NA12891_r2},
                {NA12891_r1, NA12891_r2_cram}
        };
    }

    @Test(dataProvider = "cramsWithNoReference")
    public void testCramsWithNoReference(final File file1, final File file2) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.RG.crosscheck_metrics");
        metrics.deleteOnExit();

        final String[] args = {"INPUT=" + file1.getAbsolutePath(),
                "SECOND_INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "LOD_THRESHOLD=" + -2.0,
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP_FOR_CRAMS
        };
        final CrosscheckFingerprints crossChecker = new CrosscheckFingerprints();
        Assert.assertEquals(crossChecker.instanceMain(args), 1);
    }

    @DataProvider(name = "cramBamComparison")
    public Object[][] cramBamComparison() {
        return new Object[][]{
                {NA12891_r1_shifted_bam, NA12891_r2_shifted_bam, NA12891_r1_cram, NA12891_r2_cram},
                {NA12891_r1_shifted_bam, NA12892_r1_shifted_bam, NA12891_r1_cram, NA12892_r1_cram},
                {NA12891_r1_shifted_bam, NA12892_r2_shifted_bam, NA12891_r1_cram, NA12892_r2_cram},
                {NA12892_r1_shifted_bam, NA12892_r2_shifted_bam, NA12892_r1_cram, NA12892_r2_cram},
                {NA12892_r2_shifted_bam, NA12891_r2_shifted_bam, NA12892_r2_cram, NA12891_r2_cram},
                {NA12892_r2_shifted_bam, NA12891_r1_shifted_bam, NA12892_r2_cram, NA12891_r1_cram}
        };
    }

    @Test(dataProvider = "cramBamComparison")
    public void testCramBamComparison(final File bam1, final File bam2, final File cram1, final File cram2) throws IOException {
        File metricsBam = File.createTempFile("Fingerprinting.bam.comparison", "crosscheck_metrics");
        metricsBam.deleteOnExit();
        File metricsCram = File.createTempFile("Fingerprinting.cram.comparison", "crosscheck_metrics");
        metricsCram.deleteOnExit();

        final List<String> argsBam = new ArrayList<>(Arrays.asList("INPUT=" + bam1.getAbsolutePath(),
                "INPUT=" + bam2.getAbsolutePath(),
                "OUTPUT=" + metricsBam.getAbsolutePath(),
                "LOD_THRESHOLD=" + -2.0,
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP_FOR_CRAMS)
        );

        final List<String> argsCram = new ArrayList<>(Arrays.asList("INPUT=" + cram1.getAbsolutePath(),
                "INPUT=" + cram2.getAbsolutePath(),
                "OUTPUT=" + metricsCram.getAbsolutePath(),
                "LOD_THRESHOLD=" + -2.0,
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP_FOR_CRAMS,
                "R=" + referenceForCrams)
        );

        Assert.assertEquals(runPicardCommandLine(argsBam), 0);
        final MetricsFile<CrosscheckMetric, Comparable<?>> metricsOutputBam = new MetricsFile<>();
        metricsOutputBam.read(new FileReader(metricsBam));

        Assert.assertEquals(runPicardCommandLine(argsCram), 0);
        final MetricsFile<CrosscheckMetric, Comparable<?>> metricsOutputCram = new MetricsFile<>();
        metricsOutputCram.read(new FileReader(metricsCram));

        Assert.assertEquals(metricsOutputBam.getMetrics().size(), metricsOutputCram.getMetrics().size());

        final HashMap<String, CrosscheckMetric> metricMapBam = new HashMap<>(); //lines may not be in same order in output file
        for (final CrosscheckMetric metric : metricsOutputBam.getMetrics()) {
            metricMapBam.put(metric.LEFT_GROUP_VALUE + metric.RIGHT_GROUP_VALUE, metric);
        }

        for (final CrosscheckMetric metric : metricsOutputCram.getMetrics()) {
            Assert.assertEquals(metric.LOD_SCORE, metricMapBam.get(metric.LEFT_GROUP_VALUE + metric.RIGHT_GROUP_VALUE).LOD_SCORE);
            Assert.assertEquals(metric.LOD_SCORE_NORMAL_TUMOR, metricMapBam.get(metric.LEFT_GROUP_VALUE + metric.RIGHT_GROUP_VALUE).LOD_SCORE_NORMAL_TUMOR);
            Assert.assertEquals(metric.LOD_SCORE_TUMOR_NORMAL, metricMapBam.get(metric.LEFT_GROUP_VALUE + metric.RIGHT_GROUP_VALUE).LOD_SCORE_TUMOR_NORMAL);
            Assert.assertEquals(metric.RESULT, metricMapBam.get(metric.LEFT_GROUP_VALUE + metric.RIGHT_GROUP_VALUE).RESULT);
        }
    }

    @DataProvider(name = "bamFilesLBs")
    public Object[][] bamFilesLBs() {

        return new Object[][]{
               {NA12891_r1, NA12891_r2,                   0},
               {NA12891_r1, NA12892_r1,                   0},
               {NA12892_r2, NA12891_r2,                   0},
               {NA12892_r2, NA12891_r1,                   0},
               {NA12891_r1, NA12891_named_NA12892_r1_sam, 1}, //error since expected match but found a mismatch
               {NA12892_r1, NA12891_named_NA12892_r1_sam, 1}, //error since expected mismatch but found a match
               {NA12891_r2, NA12891_named_NA12892_r1_sam, 1}, //error since expected match but found a mismatch
               {NA12892_r2, NA12891_named_NA12892_r1_sam, 1}, //error since expected mismatch but found a match
        };
    }

    @Test
    public void testCrossCheckLBsWithClustering() throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.LB.crosscheck_metrics");
        metrics.deleteOnExit();

        {
            final String[] args = new String[]{
                    "INPUT=" + NA12891_r1.getAbsolutePath(),
                    "INPUT=" + NA12891_r2.getAbsolutePath(),
                    "INPUT=" + NA12892_r1.getAbsolutePath(),
                    "INPUT=" + NA12892_r2.getAbsolutePath(),
                    "INPUT=" + NA12891_named_NA12892_r1_sam.getAbsolutePath(),
                    "OUTPUT=" + metrics.getAbsolutePath(),
                    "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                    "LOD_THRESHOLD=" + -1.0,
                    "CROSSCHECK_BY=LIBRARY"
            };
            final int nLibs = 5;
            doTest(args, metrics, 1, nLibs * nLibs, CrosscheckMetric.DataType.LIBRARY);
        }

        File cluster = File.createTempFile("Fingerprinting", "NA1291.LB.crosscheck_cluster_metrics");
        cluster.deleteOnExit();

        {
            final String[] args = new String[]{
                    "INPUT=" + metrics.getAbsolutePath(),
                    "OUTPUT=" + cluster.getAbsolutePath(),
                    "LOD_THRESHOLD=" + 1.0,
            };

            final ClusterCrosscheckMetrics clusterer = new ClusterCrosscheckMetrics();
            Assert.assertEquals(clusterer.instanceMain(args), 0);
        }

        IOUtil.assertFileIsReadable(cluster);

        final MetricsFile<ClusteredCrosscheckMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(cluster));

        final List<Long> collect = metricsOutput
                .getMetrics()
                .stream()
                .map(m -> m.CLUSTER)
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                .values()
                .stream()
                .sorted()
                .collect(Collectors.toList());

        Assert.assertEquals(collect.get(0), (Long) (2 * 2L));// need to disambiguate from assertEquals(Object, Object)
        Assert.assertEquals(collect.get(1), (Long) (3 * 3L));
    }

    @Test(dataProvider = "bamFilesLBs")
    public void testCrossCheckLBs(final File file1, final File file2, final int expectedRetVal) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.LB.crosscheck_metrics");
        metrics.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -1.0,
                "CROSSCHECK_BY=LIBRARY"
        };
        final int numLibs=2;
        doTest(args, metrics, expectedRetVal, numLibs * numLibs, CrosscheckMetric.DataType.LIBRARY);
    }

    @DataProvider(name = "bamFilesSources")
    public Object[][] bamFilesSources() {

        return new Object[][]{
                {NA12891_r1, NA12891_r2,               0},
                {NA12892_r1, NA12892_r2,               0},
                {NA12891_r1, NA12892_r1,               0},
                {NA12892_r2, NA12891_r2,               0},
                {NA12892_r2, NA12891_named_NA12892_r1, 1}, // the two files contain different samples one has wrong name
                {NA12891_r2, NA12891_named_NA12892_r1, 1}, // unexpected match
                {NA12892_r1, NA12891_named_NA12892_r1, 1}, // the two files contain different samples one has wrong name
                {NA12891_r1, NA12891_named_NA12892_r1, 1}, // unexpected match
        };
    }

    @Test(dataProvider = "bamFilesSources")
    public void testCrossCheckSources(final File file1, final File file2, final int expectedRetVal) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.Sources.crosscheck_metrics");
        metrics.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -1.0,
                "CROSSCHECK_BY=FILE"
        };

        doTest(args, metrics, expectedRetVal, 2 * 2 , CrosscheckMetric.DataType.FILE);
    }

    @DataProvider(name = "bamFilesSMs")
    public Object[][] bamFilesSMs() {

        return new Object[][] {
                {NA12891_r1, NA12891_r2,               0, 1},
                {NA12891_r1, NA12892_r1,               0, 2},
                {NA12892_r2, NA12891_r2,               0, 2},
                {NA12892_r2, NA12891_named_NA12892_r1, 0, 1}, // no error since only one sample in aggregate
                {NA12891_r2, NA12891_named_NA12892_r1, 1, 2}, // unexpected match
                {NA12892_r1, NA12891_named_NA12892_r1, 0, 1}, // no error since only one sample in aggregate
                {NA12891_r1, NA12891_named_NA12892_r1, 1, 2}, // unexpected match
        };
    }

    @Test(dataProvider = "bamFilesSMs")
    public void testCrossCheckSMs(final File file1, final File file2, final int expectedRetVal, final int numberOfSamples) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.SM.crosscheck_metrics");
        metrics.deleteOnExit();

        File matrix = File.createTempFile("Fingerprinting", "NA1291.SM.matrix");
        matrix.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "MATRIX_OUTPUT=" + matrix.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -1.0,
                "CROSSCHECK_BY=SAMPLE"
        };
        doTest(args, metrics, expectedRetVal, numberOfSamples * numberOfSamples  , CrosscheckMetric.DataType.SAMPLE);

        TabbedTextFileWithHeaderParser matrixParser = new TabbedTextFileWithHeaderParser(matrix);
        Assert.assertEquals(matrixParser.columnLabelsList().size(), numberOfSamples + 1 );
    }

    @DataProvider(name = "checkSamplesData")
    public Iterator<Object[]> checkSamplesData() {
        List<Object[]> tests = new ArrayList<>();

        // VCF tests
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_part1_vcf, NA12892_and_NA123891_part2_vcf, NA12892_and_NA123891_part3_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_vcf,NA12891_1_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12891_2_vcf, NA12892_1_vcf), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_2_vcf, NA12892_1_vcf), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_g_vcf, NA12891_named_NA12892_vcf), 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_part1_vcf, NA12892_and_NA123891_part2_vcf, NA12892_and_NA123891_part3_vcf), Arrays.asList(NA12891_2_vcf, NA12892_r1), 0, 2, true});

        // SAM vs. VCF
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1, NA12891_r2, NA12892_r2), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_r1, NA12891_r1), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_r1), Arrays.asList(NA12891_r2, NA12892_1_vcf), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12891_named_NA12892_r1), Arrays.asList(NA12891_r2, NA12892_2_vcf), 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12891_2_vcf, NA12891_named_NA12892_r1), 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_2_vcf, NA12891_named_NA12892_r1), 1, 2, false});

        // SAM tests
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), 0, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1), Arrays.asList(NA12891_r2), 1, 1, true});
        tests.add(new Object[]{Arrays.asList(NA12891_r1), Arrays.asList(NA12891_r2, NA12892_r2), 1, 1, true});
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12891_named_NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1), Arrays.asList(NA12891_r1, NA12891_named_NA12892_r1), 1, 2, false});

        // Swapped Non-ref gvcf
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12891_swapped_nonref_g_vcf, NA12892_g_vcf), 1, 2, false});

        return tests.iterator();
    }

    @Test(dataProvider = "checkSamplesData")
    public void testCheckSamples(final List<File> files1, final List<File> files2, final int expectedRetVal, final int numberOfSamples, boolean ExpectAllMatch) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "test.crosscheck_metrics");
        metrics.deleteOnExit();

        final List<String> args = new ArrayList<>();
        files1.forEach(f->args.add("INPUT="+f.getAbsolutePath()));
        files2.forEach(f->args.add("SECOND_INPUT="+f.getAbsolutePath()));

                args.add("OUTPUT=" + metrics.getAbsolutePath());
                args.add("HAPLOTYPE_MAP=" + HAPLOTYPE_MAP);
                args.add("LOD_THRESHOLD=" + -1.0);
                args.add("CROSSCHECK_BY=SAMPLE");

        doTest(args.toArray(new String[args.size()]), metrics, expectedRetVal, numberOfSamples , CrosscheckMetric.DataType.SAMPLE, ExpectAllMatch);
    }

    @DataProvider(name = "checkSamplesCrosscheckAllData")
    public Iterator<Object[]> checkSamplesCrosscheckAllData() {
        List<Object[]> tests = new ArrayList<>();

        // VCF tests
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_part1_vcf, NA12892_and_NA123891_part2_vcf, NA12892_and_NA123891_part3_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_vcf, NA12891_1_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 1, 2, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12891_2_vcf, NA12892_1_vcf), 0, 2, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_2_vcf, NA12892_1_vcf), 0, 2, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_g_vcf, NA12891_named_NA12892_vcf), 1, 2, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_part1_vcf, NA12892_and_NA123891_part2_vcf, NA12892_and_NA123891_part3_vcf), Arrays.asList(NA12891_2_vcf, NA12892_r1), 0, 2, 2, true});

        // SAM vs. VCF
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1, NA12891_r2, NA12892_r2), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 0, 2, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_r1, NA12891_r1), Arrays.asList(NA12891_g_vcf, NA12892_g_vcf), 1, 2, 2,false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_r1), Arrays.asList(NA12891_r2, NA12892_1_vcf), 0, 2, 2,true});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12891_named_NA12892_r1), Arrays.asList(NA12891_r2, NA12892_2_vcf), 1, 2, 2,false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12891_2_vcf, NA12891_named_NA12892_r1), 1, 2,2, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_2_vcf, NA12891_named_NA12892_r1), 1, 2,2, false});

        // SAM tests
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), 0, 2, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1), Arrays.asList(NA12891_r2), 0, 2, 1, true});
        tests.add(new Object[]{Arrays.asList(NA12891_r1), Arrays.asList(NA12891_r2, NA12892_r2), 0, 1, 2, true});
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12891_named_NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), 1, 2, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1), Arrays.asList(NA12891_r1, NA12891_named_NA12892_r1), 1, 2, 2, false});

        // Swapped Non-ref gvcf
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12891_swapped_nonref_g_vcf, NA12892_g_vcf), 1, 2, 2,false});

        // VCF tests
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12892_g_vcf), 0, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12892_g_vcf), 0, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_part1_vcf, NA12892_and_NA123891_part2_vcf, NA12892_and_NA123891_part3_vcf), Arrays.asList(NA12892_g_vcf), 0, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_vcf, NA12891_1_vcf), Arrays.asList(NA12892_g_vcf), 1, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Arrays.asList(NA12892_1_vcf), 0, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12892_1_vcf), 0, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12892_g_vcf), 0, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_named_NA12892_vcf), 1, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_part1_vcf, NA12892_and_NA123891_part2_vcf, NA12892_and_NA123891_part3_vcf), Arrays.asList(NA12892_r1), 0, 2, 1, false});

        // SAM vs. VCF
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1, NA12891_r2, NA12892_r2), Arrays.asList(NA12892_g_vcf), 0, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_r1, NA12891_r1), Arrays.asList(NA12892_g_vcf), 1, 2, 1, false});

        tests.add(new Object[]{Arrays.asList(NA12892_r1), Arrays.asList(NA12891_r2, NA12892_1_vcf), 0, 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_r1), Arrays.asList(NA12891_r2, NA12892_2_vcf), 1, 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12892_1_vcf), Arrays.asList(NA12891_2_vcf, NA12891_named_NA12892_r1), 1, 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf), Arrays.asList(NA12891_named_NA12892_r1), 1, 2, 1, false});

        // SAM tests
        tests.add(new Object[]{Arrays.asList( NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), 0, 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12892_r1), Arrays.asList(NA12891_r2), 0, 1, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12891_r1), Arrays.asList(NA12892_r2), 0, 1, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), 1, 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12892_r1), Arrays.asList(NA12891_r1, NA12891_named_NA12892_r1), 1, 1, 2, false});

        // Swapped Non-ref gvcf
        tests.add(new Object[]{Arrays.asList(NA12892_1_vcf), Arrays.asList(NA12891_swapped_nonref_g_vcf, NA12892_g_vcf), 0, 1, 2, false});

        return tests.iterator();
    }

    @DataProvider
    public Iterator<Object[]> checkSamplesCrosscheckAllWithMappingData() {
        List<Object[]> tests = new ArrayList<>();

        File NA12892_to_NA12891 = new File(TEST_DATA_DIR,"NA12892_to_NA12891.txt");
        File NotThere_to_NA12892 = new File(TEST_DATA_DIR,"NotThere_to_NA12891.txt");

        tests.add(new Object[]{Arrays.asList( NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), NotThere_to_NA12892, null, 0, 1, 2, false});
        tests.add(new Object[]{Arrays.asList( NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), null, NotThere_to_NA12892, 0, 1, 2, false});
        tests.add(new Object[]{Arrays.asList( NA12892_r1), Arrays.asList(NA12891_r2, NA12892_r2), NotThere_to_NA12892, NotThere_to_NA12892, 0, 1, 2, false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Collections.singletonList(NA12892_g_vcf), null, null, 0, 2, 1, false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf), Collections.singletonList(NA12892_g_vcf), null, NA12892_to_NA12891, 1, 2, 1, false});
        tests.add(new Object[]{Collections.singletonList(NA12892_1_vcf), Collections.singletonList(NA12892_g_vcf), NA12892_to_NA12891, NA12892_to_NA12891, 0, 1, 1, false});
        tests.add(new Object[]{Collections.singletonList(NA12891_named_NA12892_r1), Collections.singletonList(NA12891_r2), NA12892_to_NA12891, null, 0, 1, 1, false});
        tests.add(new Object[]{Collections.singletonList(NA12892_and_NA123891_vcf), Collections.singletonList(NA12891_named_NA12892_r1),null, NA12892_to_NA12891, 0, 2, 1, false});

        return tests.iterator();
    }

    @Test(dataProvider = "checkSamplesCrosscheckAllWithMappingData")
    public void testSecondInputCheckAllWithMapping(final List<File> files1, final List<File> files2, final File inputSampleMap, final File secondInputSampleMap, final int expectedRetVal, final int numberOfSamples1, final int numberOfSamples2, boolean ExpectAllMatch) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "test.crosscheck_metrics");
        metrics.deleteOnExit();

        final List<String> args = new ArrayList<>();
        files1.forEach(f->args.add("INPUT="+f.getAbsolutePath()));
        files2.forEach(f->args.add("SECOND_INPUT="+f.getAbsolutePath()));

        args.add("OUTPUT=" + metrics.getAbsolutePath());
        args.add("HAPLOTYPE_MAP=" + HAPLOTYPE_MAP);
        args.add("LOD_THRESHOLD=" + -1.0);
        args.add("CROSSCHECK_BY=SAMPLE");
        args.add("CROSSCHECK_MODE=CHECK_ALL_OTHERS");

        if(inputSampleMap!=null)  args.add("INPUT_SAMPLE_MAP="+inputSampleMap.getAbsolutePath());
        if(secondInputSampleMap!=null)  args.add("SECOND_INPUT_SAMPLE_MAP="+secondInputSampleMap.getAbsolutePath());

        doTest(args.toArray(new String[args.size()]), metrics, expectedRetVal, numberOfSamples1 * numberOfSamples2 , CrosscheckMetric.DataType.SAMPLE, ExpectAllMatch);
    }

    @DataProvider
    public Iterator<Object[]> checkSampleMapFailuresData() {
        List<Object[]> tests = new ArrayList<>();

        File NA12892_to_NA12891 = new File(TEST_DATA_DIR,"NA12892_to_NA12891.txt");
        File NA12891_to_NA12892 = new File(TEST_DATA_DIR,"NA12891_to_NA12892.txt");
        File tooManyFields = new File(TEST_DATA_DIR,"too_many_fields.txt");
        File tooFewFields = new File(TEST_DATA_DIR,"too_few_fields.txt");
        File duplicateKey = new File(TEST_DATA_DIR,"DuplicateKey.txt");

        tests.add(new Object[]{NA12891_to_NA12892, null});// NA12892 will become non-unique in first file
        tests.add(new Object[]{NA12892_to_NA12891, null});// NA12891 will become non-unique in first file
        tests.add(new Object[]{tooManyFields, null});// file too_many_fields.txt has too many fields
        tests.add(new Object[]{tooFewFields, null});// file too_few_fields.txt has too few fields
        tests.add(new Object[]{duplicateKey, null});// file too_few_fields.txt has too few fields

        tests.add(new Object[]{null, tooManyFields});// file too_many_fields.txt has too many fields
        tests.add(new Object[]{null, tooFewFields});// file too_few_fields.txt has too few fields
        tests.add(new Object[]{null, duplicateKey});// file DuplicateKey.txt has a duplicate key

        return tests.iterator();
    }

    @Test(dataProvider = "checkSampleMapFailuresData", expectedExceptions = IllegalArgumentException.class)
    public void checkSampleMapFailuresData( final File inputSampleMap, final File secondInputSampleMap) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "test.crosscheck_metrics");
        metrics.deleteOnExit();

        final List<File> files1 = Arrays.asList(NA12891_1_vcf, NA12892_1_vcf);
        final List<File> files2 = Collections.singletonList(NA12892_g_vcf);

        final List<String> args = new ArrayList<>();
        files1.forEach(f->args.add("INPUT="+f.getAbsolutePath()));
        files2.forEach(f->args.add("SECOND_INPUT="+f.getAbsolutePath()));

        args.add("OUTPUT=" + metrics.getAbsolutePath());
        args.add("HAPLOTYPE_MAP=" + HAPLOTYPE_MAP);
        args.add("LOD_THRESHOLD=" + -1.0);
        args.add("CROSSCHECK_BY=SAMPLE");
        args.add("CROSSCHECK_MODE=CHECK_ALL_OTHERS");

        if(inputSampleMap!=null)  args.add("INPUT_SAMPLE_MAP="+inputSampleMap.getAbsolutePath());
        if(secondInputSampleMap!=null)  args.add("SECOND_INPUT_SAMPLE_MAP="+secondInputSampleMap.getAbsolutePath());

        doTest(args.toArray(new String[args.size()]), metrics, 0, 0 , CrosscheckMetric.DataType.SAMPLE, false);
    }

    @Test(dataProvider = "checkSamplesCrosscheckAllData")
    public void testSecondInputCheckAll(final List<File> files1, final List<File> files2, final int expectedRetVal, final int numberOfSamples1, final int numberOfSamples2, boolean ExpectAllMatch) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "test.crosscheck_metrics");
        metrics.deleteOnExit();

        final List<String> args = new ArrayList<>();
        files1.forEach(f->args.add("INPUT="+f.getAbsolutePath()));
        files2.forEach(f->args.add("SECOND_INPUT="+f.getAbsolutePath()));

                args.add("OUTPUT=" + metrics.getAbsolutePath());
                args.add("HAPLOTYPE_MAP=" + HAPLOTYPE_MAP);
                args.add("LOD_THRESHOLD=" + -1.0);
                args.add("CROSSCHECK_BY=SAMPLE");
                args.add("CROSSCHECK_MODE=CHECK_ALL_OTHERS");

        doTest(args.toArray(new String[args.size()]), metrics, expectedRetVal, numberOfSamples1 * numberOfSamples2 , CrosscheckMetric.DataType.SAMPLE, ExpectAllMatch);
    }

    @DataProvider(name = "checkFilesData")
    public Iterator<Object[]> checkFilesData() {
        List<Object[]> tests = new ArrayList<>();

        // VCF tests
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf, NA12891_g_vcf, NA12892_g_vcf), 0, 4 * 4, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf, NA12891_g_vcf, NA12892_g_vcf), 0, 4 * 4, false});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_vcf, NA12891_1_vcf, NA12891_g_vcf, NA12892_g_vcf), 1, 4 * 4, false});
        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_1_vcf, NA12891_2_vcf, NA12892_1_vcf), 0, 3 * 3, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf, NA12891_2_vcf, NA12892_1_vcf), 0, 4 * 4, false});
        tests.add(new Object[]{Arrays.asList(NA12892_and_NA123891_vcf, NA12891_g_vcf, NA12891_named_NA12892_vcf), 1, 4 * 4, false});

        // SAM vs. VCF
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1, NA12891_r2, NA12892_r2, NA12891_g_vcf, NA12892_g_vcf), 0, 6 * 6, false});
        tests.add(new Object[]{Arrays.asList(NA12891_named_NA12892_r1, NA12891_r1, NA12891_g_vcf, NA12891_1_vcf), 1, 4*4, false});

        tests.add(new Object[]{Arrays.asList(NA12891_1_vcf, NA12892_r1, NA12891_r2, NA12892_1_vcf), 0, 4*4, true});

        // SAM tests
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1, NA12891_r2, NA12892_r2), 0, 4*4, true});
        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12892_r1, NA12891_r2), 0, 3*3, true});

        tests.add(new Object[]{Arrays.asList(NA12891_r1, NA12891_named_NA12892_r1, NA12891_r2, NA12892_r2), 1, 4*4, false});

        return tests.iterator();
    }

    @Test(dataProvider = "checkFilesData")
    public void testCheckFiles(final List<File> files, final int expectedRetVal, final int numberOfSamples, boolean ExpectAllMatch) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "test.crosscheck_metrics");
        metrics.deleteOnExit();

        final List<String> args = new ArrayList<>();
        files.forEach(f->args.add("INPUT="+f.getAbsolutePath()));

        args.add("OUTPUT=" + metrics.getAbsolutePath());
        args.add("HAPLOTYPE_MAP=" + HAPLOTYPE_MAP);
        args.add("LOD_THRESHOLD=" + -1.0);
        args.add("CROSSCHECK_BY=FILE");

        doTest(args.toArray(new String[args.size()]), metrics, expectedRetVal, numberOfSamples , CrosscheckMetric.DataType.FILE, ExpectAllMatch);
    }

    @DataProvider(name = "checkPathsData")
    public Iterator<Object[]> checkPathsData() throws IOException {
        List<Object[]> tests = new ArrayList<>();

        final File fofn = File.createTempFile("crosscheck",".fofn", TEST_DATA_DIR);
        fofn.deleteOnExit();

        PrintWriter writer = new PrintWriter(new FileWriter(fofn));
        writer.println(NA12892_1_vcf.getAbsolutePath());
        writer.println(NA12891_g_vcf.getAbsolutePath());
        writer.println(NA12892_g_vcf.getAbsolutePath());
        final int fofnSize = 3;
        writer.close();

        // VCF tests
        tests.add(new Object[]{Stream.of(NA12891_1_vcf, NA12892_1_vcf, NA12891_g_vcf, NA12892_g_vcf).map(File::getAbsolutePath).collect(Collectors.toList()), HAPLOTYPE_MAP.getAbsolutePath(), 0, 4 * 4, false});
        tests.add(new Object[]{Stream.of(NA12892_1_vcf, NA12891_g_vcf, NA12892_g_vcf).map(File::getAbsolutePath).collect(Collectors.toList()), HAPLOTYPE_MAP.getAbsolutePath(), 0, 3 * 3, false});
        tests.add(new Object[]{Collections.singletonList(fofn.getAbsolutePath()), HAPLOTYPE_MAP.getAbsolutePath(), 0, fofnSize * fofnSize, false});

        return tests.iterator();
    }

    @Test(dataProvider = "checkPathsData")
    public void testCheckPaths(final List<String> paths, final String haploypeMap, final int expectedRetVal, final int numberOfSamples, boolean ExpectAllMatch) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "test.crosscheck_metrics");
        metrics.deleteOnExit();

        final List<String> args = new ArrayList<>();

        paths.forEach(s -> args.add("INPUT=" + s));

        args.add("OUTPUT=" + metrics.getAbsolutePath());
        args.add("HAPLOTYPE_MAP=" + haploypeMap);
        args.add("LOD_THRESHOLD=" + -1.0);
        args.add("CROSSCHECK_BY=FILE");

        doTest(args.toArray(new String[args.size()]), metrics, expectedRetVal, numberOfSamples , CrosscheckMetric.DataType.FILE, ExpectAllMatch);
    }

    @DataProvider(name = "noFingerprintingSitesData")
    public Object[][] noFingerprintingSitesVCFData() {
        return new Object[][]{
                {NA12891_no_fp_sites_vcf, null, HAPLOTYPE_MAP},
                {NA12891_no_fp_sites_vcf, NA12891_1_vcf, HAPLOTYPE_MAP},
                {NA12891_1_vcf, NA12891_no_fp_sites_vcf, HAPLOTYPE_MAP},
                {NA12891_r1, null, HAPLOTYPE_MAP_MOVED}
        };
    }

    @Test(dataProvider = "noFingerprintingSitesData")
    public void testNoFingerprintingSites(final File input, final File second_input, final File hap_map) throws IOException {
        File metrics = File.createTempFile("Fingerprinting.comparison", "crosscheck_metrics");
        metrics.deleteOnExit();
        final List<String> args = new ArrayList<>(Arrays.asList("INPUT=" + input,
                "OUTPUT=" + metrics.getAbsolutePath(),
                "LOD_THRESHOLD=" + -2.0,
                "HAPLOTYPE_MAP=" + hap_map.getAbsolutePath())
        );

        if (second_input != null) {
            args.add("SECOND_INPUT=" + second_input);
        }

        Assert.assertEquals(runPicardCommandLine(args), 1);
        args.add("CROSSCHECK_MODE=CHECK_ALL_OTHERS");

        Assert.assertEquals(runPicardCommandLine(args), 1);

        for (final CrosscheckMetric.DataType by : CrosscheckMetric.DataType.values()) {
            args.add("CROSSCHECK_BY=" + by);
            Assert.assertEquals(runPicardCommandLine(args), 1);
            args.remove(args.size() - 1);
        }
    }

    @DataProvider(name = "someGroupNoFingerprintingSitesData")
    public Object[][] someGroupNoFingerprintingSitesData() {
        return new Object[][]{
                {NA12891_no_fp_sites_and_NA12892_vcf, NA12892_and_NA123891_vcf, HAPLOTYPE_MAP, 0, 2, CrosscheckMetric.DataType.SAMPLE},
                {NA12891_r1_one_rg_no_fingerprint_sam, null, HAPLOTYPE_MAP, 0, NA12891_r1_RGs * NA12891_r1_RGs, CrosscheckMetric.DataType.READGROUP}
        };
    }

    @Test(dataProvider = "someGroupNoFingerprintingSitesData")
    public void testSomeGroupNoFingerprintingSites(final File input, final File second_input, final File hap_map, final int exptectRetVal, final int expectedNMetrics, final CrosscheckMetric.DataType dataType) throws IOException {
        File metrics = File.createTempFile("Fingerprinting.comparison", "crosscheck_metrics");
        metrics.deleteOnExit();
        final List<String> args = new ArrayList<>(Arrays.asList("INPUT=" + input,
                "OUTPUT=" + metrics.getAbsolutePath(),
                "LOD_THRESHOLD=" + -1.0,
                "HAPLOTYPE_MAP=" + hap_map.getAbsolutePath())
        );

        if (second_input != null) {
            args.add("SECOND_INPUT=" + second_input);
        }

        doTest(args.toArray(new String[args.size()]), metrics, exptectRetVal, expectedNMetrics, dataType, true);
    }

    private void doTest(final String[] args, final File metrics, final int expectedRetVal, final int expectedNMetrics, final CrosscheckMetric.DataType expectedType) throws IOException {
        doTest(args, metrics, expectedRetVal, expectedNMetrics, expectedType, false);
    }

    private void doTest(final String[] args, final File metrics, final int expectedRetVal, final int expectedNMetrics, final CrosscheckMetric.DataType expectedType, final boolean expectAllMatch) throws IOException {
       
        final CrosscheckFingerprints crossChecker = new CrosscheckFingerprints();
        Assert.assertEquals(crossChecker.instanceMain(args), expectedRetVal);

        final MetricsFile<CrosscheckMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(metrics));

        Assert.assertFalse(metricsOutput.getMetrics().stream()
                .anyMatch(m -> m.DATA_TYPE != expectedType));

        Assert.assertFalse(metricsOutput.getMetrics().stream()
                .anyMatch(m-> m.LOD_SCORE_NORMAL_TUMOR == null));
        Assert.assertFalse(metricsOutput.getMetrics().stream()
                .anyMatch(m-> m.LOD_SCORE == null));
        Assert.assertFalse(metricsOutput.getMetrics().stream()
                .anyMatch(m-> m.LOD_SCORE_TUMOR_NORMAL == null));

        if (expectAllMatch) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .allMatch(m -> m.RESULT == CrosscheckMetric.FingerprintResult.INCONCLUSIVE ||
                            m.RESULT.isMatch() == m.LEFT_SAMPLE.equals(m.RIGHT_SAMPLE)));
        } else if (expectedRetVal == 0) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .allMatch(m -> m.RESULT == CrosscheckMetric.FingerprintResult.INCONCLUSIVE ||
                            m.RESULT.isExpected()));
        } else {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .anyMatch(m -> !m.RESULT.isExpected()));
        }

        Assert.assertEquals(metricsOutput.getMetrics().size(), expectedNMetrics);

        // at the readgroup level things aren't always conclusive...
        if (!metricsOutput.getMetrics().isEmpty() && expectedType != CrosscheckMetric.DataType.READGROUP) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .anyMatch(m -> m.RESULT != CrosscheckMetric.FingerprintResult.INCONCLUSIVE));
        }

        //check that fields that should have an actual value, indeed do
        for(final String fieldName : lookupMap.get(expectedType)) {
            try {
                final Field field = CrosscheckMetric.class.getField(fieldName);
                Assert.assertTrue(metricsOutput.getMetrics().stream().allMatch(m -> {
                    try {
                        return field.get(m) != FingerprintIdDetails.multipleValuesString && field.get(m) != null;
                    } catch (IllegalAccessException e) {
                        e.printStackTrace();
                        return false;
                    }
                }));

            } catch (NoSuchFieldException e) {
                e.printStackTrace();
                assert false;
            }
        }
    }

    @Test
    public void canWriteToDevNull() throws IOException {
        File f = new File("/dev/null");
        Assert.assertTrue(f.canRead());

        final OutputStream stream = new FileOutputStream(f);
        final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stream));

        writer.write("Just a test");
        writer.close();
    }
}