package picard.arrays;

import picard.vcf.VcfTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class MergePedIntoVcfTest {
    private static final Path TEST_DATA_DIR = Paths.get("testdata/picard/arrays/");

    @Test()
    public void testMergePedIntoVcf() throws IOException {
        final MergePedIntoVcf mergePedIntoVcf = new MergePedIntoVcf();
        final File inputVcf = TEST_DATA_DIR.resolve("input.vcf").toFile();
        final File pedFile = TEST_DATA_DIR.resolve("zcall.output.ped").toFile();
        final File mapFile = TEST_DATA_DIR.resolve("zcall.output.map").toFile();
        final File zcallThresholdsFile = TEST_DATA_DIR.resolve("zcall.thresholds.7.txt").toFile();

        final File output = File.createTempFile("output", ".vcf");
        output.deleteOnExit();
        mergePedIntoVcf.ORIGINAL_VCF = inputVcf;
        mergePedIntoVcf.PED_FILE = pedFile;
        mergePedIntoVcf.MAP_FILE = mapFile;
        mergePedIntoVcf.ZCALL_THRESHOLDS_FILE = zcallThresholdsFile;
        mergePedIntoVcf.ZCALL_VERSION = "1.0.0.0";
        mergePedIntoVcf.OUTPUT = output;

        Assert.assertEquals(mergePedIntoVcf.doWork(), 0);

        final File expected = TEST_DATA_DIR.resolve("expected_output.vcf").toFile();
        VcfTestUtils.assertVcfFilesAreEqual(output, expected);
    }
}