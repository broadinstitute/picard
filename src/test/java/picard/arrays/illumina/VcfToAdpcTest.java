package picard.arrays.illumina;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.arrays.VcfToAdpc;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class VcfToAdpcTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_VCF = new File(TEST_DATA_DIR, "TestVcfToAdpc.vcf");
    private static final File SINGLE_SAMPLE_VCF = new File(TEST_DATA_DIR, "TestAdpc1.vcf");
    private static final File MULTI_SAMPLE_VCF = new File(TEST_DATA_DIR, "TestAdpc23.vcf");

    private static final File EXPECTED_ADPC_BIN_FILE = new File(TEST_DATA_DIR, "TestVcfToAdpc.adpc.bin");
    private static final File EXPECTED_SAMPLES_FILE = new File(TEST_DATA_DIR, "TestVcfToAdpc.samples.txt");
    private static final File EXPECTED_NUM_MARKERS_FILE = new File(TEST_DATA_DIR, "TestVcfToAdpc.num_markers.txt");

    // Test with a single sample VCF as input.
    private static final File EXPECTED_SINGLE_SAMPLE_ADPC_BIN_FILE = new File(TEST_DATA_DIR, "TestAdpc1.adpc.bin");
    private static final File EXPECTED_SINGLE_SAMPLE_SAMPLES_FILE = new File(TEST_DATA_DIR, "TestAdpc1.samples.txt");
    private static final File EXPECTED_SINGLE_SAMPLE_NUM_MARKERS_FILE = new File(TEST_DATA_DIR, "TestAdpc.num_markers.txt");
    // Test with a multi sample VCF as input
    private static final File EXPECTED_MULTI_SAMPLE_ADPC_BIN_FILE = new File(TEST_DATA_DIR, "TestAdpc23.adpc.bin");
    private static final File EXPECTED_MULTI_SAMPLE_SAMPLES_FILE = new File(TEST_DATA_DIR, "TestAdpc23.samples.txt");
    private static final File EXPECTED_MULTI_SAMPLE_NUM_MARKERS_FILE = new File(TEST_DATA_DIR, "TestAdpc.num_markers.txt");
    // Test with both a single sample VCF and multi sample VCF as input.  In that order
    private static final File EXPECTED_S_M_ADPC_BIN_FILE = new File(TEST_DATA_DIR, "TestAdpc1_23.adpc.bin");
    private static final File EXPECTED_S_M_SAMPLES_FILE = new File(TEST_DATA_DIR, "TestAdpc1_23.samples.txt");
    private static final File EXPECTED_S_M_NUM_MARKERS_FILE = new File(TEST_DATA_DIR, "TestAdpc.num_markers.txt");
    // Test with both a multi sample VCF and single sample VCF as input.  In that order
    private static final File EXPECTED_M_S_ADPC_BIN_FILE = new File(TEST_DATA_DIR, "TestAdpc23_1.adpc.bin");
    private static final File EXPECTED_M_S_SAMPLES_FILE = new File(TEST_DATA_DIR, "TestAdpc23_1.samples.txt");
    private static final File EXPECTED_M_S_NUM_MARKERS_FILE = new File(TEST_DATA_DIR, "TestAdpc.num_markers.txt");

    @DataProvider(name = "vcfToAdpcBinCombinations")
    public Object[][] vcfToAdpcBinCombinations() {
        return new Object[][]{
                {Collections.singletonList(TEST_VCF), EXPECTED_ADPC_BIN_FILE, EXPECTED_SAMPLES_FILE, EXPECTED_NUM_MARKERS_FILE},
                {Collections.singletonList(SINGLE_SAMPLE_VCF), EXPECTED_SINGLE_SAMPLE_ADPC_BIN_FILE, EXPECTED_SINGLE_SAMPLE_SAMPLES_FILE, EXPECTED_SINGLE_SAMPLE_NUM_MARKERS_FILE},
                {Collections.singletonList(MULTI_SAMPLE_VCF), EXPECTED_MULTI_SAMPLE_ADPC_BIN_FILE, EXPECTED_MULTI_SAMPLE_SAMPLES_FILE, EXPECTED_MULTI_SAMPLE_NUM_MARKERS_FILE},
                {Arrays.asList(SINGLE_SAMPLE_VCF, MULTI_SAMPLE_VCF), EXPECTED_S_M_ADPC_BIN_FILE, EXPECTED_S_M_SAMPLES_FILE, EXPECTED_S_M_NUM_MARKERS_FILE},
                {Arrays.asList(MULTI_SAMPLE_VCF, SINGLE_SAMPLE_VCF), EXPECTED_M_S_ADPC_BIN_FILE, EXPECTED_M_S_SAMPLES_FILE, EXPECTED_M_S_NUM_MARKERS_FILE}
        };
    }

    @Test(dataProvider = "vcfToAdpcBinCombinations")
    public void testVcfToAdpc(final List<File> vcfs, final File expectedAdpcBinFile, final File expectedSamplesFile, final File expectedNumMarkersFile) throws IOException {
        final File output = File.createTempFile("testIlluminaAdpcFileWriter.", ".adpc.bin");
        output.deleteOnExit();
        final File samplesFile = File.createTempFile("testIlluminaAdpcFileWriter.", ".samples.txt");
        samplesFile.deleteOnExit();
        final File numMarkersFile = File.createTempFile("testIlluminaAdpcFileWriter.", ".num_markers.txt");
        samplesFile.deleteOnExit();

        final VcfToAdpc vcfToAdpc = new VcfToAdpc();
        vcfToAdpc.VCF = vcfs;
        vcfToAdpc.OUTPUT = output;
        vcfToAdpc.SAMPLES_FILE = samplesFile;
        vcfToAdpc.NUM_MARKERS_FILE = numMarkersFile;

        Assert.assertEquals(vcfToAdpc.instanceMain(new String[0]), 0);

        IOUtil.assertFilesEqual(expectedAdpcBinFile, output);
        IOUtil.assertFilesEqual(expectedSamplesFile, samplesFile);
        IOUtil.assertFilesEqual(expectedNumMarkersFile, numMarkersFile);
    }

    @Test
    public void testVcfToAdpcFailOnDifferingNumberOfLoci() throws IOException {
        final File output = File.createTempFile("testIlluminaAdpcFileWriter.", ".adpc.bin");
        output.deleteOnExit();
        final File samplesFile = File.createTempFile("testIlluminaAdpcFileWriter.", ".samples.txt");
        samplesFile.deleteOnExit();
        final File numMarkersFile = File.createTempFile("testIlluminaAdpcFileWriter.", ".num_markers.txt");
        samplesFile.deleteOnExit();

        final VcfToAdpc vcfToAdpc = new VcfToAdpc();
        vcfToAdpc.VCF = Arrays.asList(TEST_VCF, SINGLE_SAMPLE_VCF);
        vcfToAdpc.OUTPUT = output;
        vcfToAdpc.SAMPLES_FILE = samplesFile;
        vcfToAdpc.NUM_MARKERS_FILE = numMarkersFile;

        Assert.assertEquals(vcfToAdpc.instanceMain(new String[0]), 1);
    }
}
