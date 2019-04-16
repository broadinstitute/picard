package picard.sam;

import htsjdk.samtools.SAMException;
import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class BuildBamIndexTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/indices");
    private static final File OUTPUT_FILE = new File(TEST_DATA_DIR, "index_test.bam.bai");
    private static final File EXPECTED_FILE = new File(TEST_DATA_DIR, "index_test.bam.bai.exp");

    public String getCommandLineProgramName() { return BuildBamIndex.class.getSimpleName(); }

    @DataProvider
    public Object[][] indexFilesOK() {
        return new Object[][] {
                {"index_test0.bam"}, // BAM file with no index attached
                {"index_test1.bam"}, // BAM file with BAI index
                {"index_test2.bam"}  // BAM file with CSI index
        };
    }

    @DataProvider
    public Object[][] indexFilesFail() {
        return new Object[][] {
                {"index_test_fail.sam"}, // SAM file
                {"index_test_fail.bam"}  // Unsorted BAM file
        };
    }

    // Test that the index file is created when there is no index file
    @Test(dataProvider = "indexFilesOK")
    public void testBuildBamIndexOK(String bamFileName) throws IOException {
        final List<String> args = new ArrayList<String>();
        args.add("INPUT=" + new File(TEST_DATA_DIR, bamFileName).getAbsolutePath());
        args.add("OUTPUT=" + OUTPUT_FILE.getAbsolutePath());
        runPicardCommandLine(args);
        Assert.assertEquals(FileUtils.readFileToByteArray(OUTPUT_FILE),  FileUtils.readFileToByteArray(EXPECTED_FILE));
        FileUtils.forceDelete(OUTPUT_FILE);
    }

    // Test that the index creation fails when presented with a SAM file
    @Test(dataProvider = "indexFilesFail", expectedExceptions = SAMException.class)
    public void testBuildBamIndexFail(String bamFileName) throws IOException {
        final List<String> args = new ArrayList<String>();
        args.add("INPUT=" + new File(TEST_DATA_DIR, bamFileName).getAbsolutePath());
        args.add("OUTPUT=" + OUTPUT_FILE.getAbsolutePath());
        runPicardCommandLine(args);
        FileUtils.forceDelete(OUTPUT_FILE);
    }
}
