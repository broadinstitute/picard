package picard.vcf;

import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by farjoun on 10/22/17.
 */
public class GatherVcfsTest extends CommandLineProgramTest {


    @Override
    public String getCommandLineProgramName() {
        return GatherVcfs.class.getSimpleName();
    }

    private File vcf,shard1, shard2, shard3, shard2_bad;
    private File vcf_gz,shard1_gz, shard2_gz, shard3_gz, shard2_bad_gz;

    @BeforeClass
    public void setup() throws IOException {
        final File TEST_DIR = new File("testdata/picard/vcf/GatherVcf");
        vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"input.vcf"),"whole.");
        shard1 = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"shard1.vcf"),"shard1.");
        shard2 = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"shard2.vcf"),"shard2.");
        shard2_bad = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"shard2_bad.vcf"),"shard2_bad.");
        shard3 = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"shard3.vcf"),"shard3");

        vcf_gz = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"input.vcf"),"whole.",".vcf.gz");
        shard1_gz = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"shard1.vcf"),"shard1.",".vcf.gz");
        shard2_gz = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"shard2.vcf"),"shard2.",".vcf.gz");
        shard2_bad_gz = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"shard2_bad.vcf"),"shard2_bad.",".vcf.gz");
        shard3_gz = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DIR,"shard3.vcf"),"shard3.",".vcf.gz");

    }

    @DataProvider
    public Object[][] vcfshards() {
        return new Object[][]{
                {Arrays.asList(shard1,shard2,shard3), vcf, 0},
                {Arrays.asList(shard1,shard2_bad,shard3), vcf, 1},
                {Arrays.asList(shard1,shard3,shard2), vcf, 1},
                {Arrays.asList(shard3,shard1,shard2), vcf, 1} ,
                {Arrays.asList(shard1_gz, shard2_gz, shard3), vcf_gz, 0},
                {Arrays.asList(shard1_gz, shard2_bad_gz, shard3), vcf_gz, 1},
                {Arrays.asList(shard1_gz, shard3_gz, shard2), vcf_gz, 1},
                {Arrays.asList(shard3_gz, shard1_gz, shard2), vcf_gz, 1}
        };
    }

    @Test(dataProvider = "vcfshards")
    public void TestGatherFiles(final List<File> inputFiles, final File expectedOutput, final int expectedRetVal) throws IOException {

        final List<String> args = new ArrayList<>();

        final File output = VcfTestUtils.createTemporaryIndexedFile("result", expectedOutput.getAbsolutePath().endsWith(".vcf") ? ".vcf" : ".vcf.gz");

        inputFiles.forEach(f -> args.add("INPUT=" + f.getAbsolutePath()));
        args.add("OUTPUT=" + output.getAbsolutePath());

        Assert.assertEquals(runPicardCommandLine(args.toArray(new String[args.size()])), expectedRetVal, "Program was expected to run successfully, but didn't.");

        if (expectedRetVal == 0) {
            final VCFFileReader expectedReader = new VCFFileReader(expectedOutput, false);
            final VCFFileReader outputReader = new VCFFileReader(output, false);
            Assert.assertEquals(expectedReader.iterator().stream().count(), outputReader.iterator().stream().count(), "The wrong number of variants was found.");
        }
    }
}