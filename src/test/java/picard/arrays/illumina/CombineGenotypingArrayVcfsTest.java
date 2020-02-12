package picard.arrays.illumina;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.arrays.CombineGenotypingArrayVcfs;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CombineGenotypingArrayVcfsTest {
    private static final Path TEST_DATA_DIR = Paths.get("testdata/picard/arrays/");

    @Test()
    public void testCombineGenotypingArrayVcfs() throws IOException {
        final CombineGenotypingArrayVcfs combineGenotypingArrayVcfs = new CombineGenotypingArrayVcfs();
        final File inputVcf1 = TEST_DATA_DIR.resolve("input_for_combine.vcf").toFile();
        final File inputVcf2 = TEST_DATA_DIR.resolve("input2_for_combine.vcf").toFile();
        final List<File> inputs = Arrays.asList(inputVcf1, inputVcf2);

        final File output = File.createTempFile("output", ".vcf");
        output.deleteOnExit();
        combineGenotypingArrayVcfs.INPUT = inputs;
        combineGenotypingArrayVcfs.OUTPUT = output;

        Assert.assertEquals(combineGenotypingArrayVcfs.doWork(), 0);

        final File expected = TEST_DATA_DIR.resolve("combined_output.vcf").toFile();
        VcfTestUtils.assertVcfFilesAreEqual(output, expected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCombineGenotypingArrayVcfsFail() throws IOException {
        final CombineGenotypingArrayVcfs combineGenotypingArrayVcfs = new CombineGenotypingArrayVcfs();
        final File inputVcf1 = TEST_DATA_DIR.resolve("input_for_combine.vcf").toFile();
        final File inputVcf2 = TEST_DATA_DIR.resolve("input_for_combine.vcf").toFile();
        final List<File> inputs = Arrays.asList(inputVcf1, inputVcf2);

        final File output = File.createTempFile("output", ".vcf");
        output.deleteOnExit();
        combineGenotypingArrayVcfs.INPUT = inputs;
        combineGenotypingArrayVcfs.OUTPUT = output;

        Assert.assertEquals(combineGenotypingArrayVcfs.doWork(), 1);
    }
}
