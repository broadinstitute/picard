package picard.arrays.illumina;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.arrays.VcfToAdpc;

import java.io.File;
import java.io.IOException;

public class VcfToAdpcTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_VCF = new File(TEST_DATA_DIR, "TestVcfToAdpc.vcf");
    private static final File TEST_EXPECTED_ADPC_BIN_FILE = new File(TEST_DATA_DIR, "TestIlluminaAdpcFileWriter.adpc.bin");

    @Test
    public void testVcfToAdpc() throws IOException {
        final File output = File.createTempFile("testIlluminaAdpcFileWriter.", ".adpc.bin");
        output.deleteOnExit();

        final VcfToAdpc vcfToAdpc = new VcfToAdpc();
        vcfToAdpc.VCF = TEST_VCF;
        vcfToAdpc.OUTPUT = output;

        Assert.assertEquals(vcfToAdpc.instanceMain(new String[0]), 0);

        IOUtil.assertFilesEqual(TEST_EXPECTED_ADPC_BIN_FILE, output);
    }
}
