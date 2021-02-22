package picard.sam;

import htsjdk.samtools.SAMException;
import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.AfterTest;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class BuildBamIndexTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/indices/");
    private static final File INPUT_FILE = new File(TEST_DATA_DIR, "index_test.sam");
    private static final File OUTPUT_SORTED_FILE = new File(TEST_DATA_DIR, "index_test_sorted.bam");
    private static final File OUTPUT_UNSORTED_FILE = new File(TEST_DATA_DIR, "/index_test_unsorted.bam");
    private static final File OUTPUT_INDEX_FILE = new File(TEST_DATA_DIR, "/index_test.bam.bai");
    private static final File EXPECTED_BAI_FILE = new File(TEST_DATA_DIR, "index_test_b.bam.bai");

    public String getCommandLineProgramName() { return BuildBamIndex.class.getSimpleName(); }


    // Test that the index file for a sorted BAM is created
    @Test
    public void testBuildBamIndexOK() throws IOException {
        final List<String> args = new ArrayList<String>();
        /* First sort, before indexing */
        new SortSam().instanceMain(new String[]{
                "I=" + INPUT_FILE,
                "O=" + OUTPUT_SORTED_FILE,
                "SORT_ORDER=coordinate"});

        args.add("INPUT=" + OUTPUT_SORTED_FILE);
        args.add("OUTPUT=" + OUTPUT_INDEX_FILE);
        runPicardCommandLine(args);
        Assert.assertEquals(FileUtils.readFileToByteArray(OUTPUT_INDEX_FILE),  FileUtils.readFileToByteArray(EXPECTED_BAI_FILE));
    }

    // Test that the index creation fails when presented with a SAM file
    @Test(expectedExceptions = SAMException.class)
    public void testBuildSamIndexFail() {
        final List<String> args = new ArrayList<String>();
        args.add("INPUT=" + INPUT_FILE);
        args.add("OUTPUT=" + OUTPUT_INDEX_FILE);
        runPicardCommandLine(args);
    }

    // Test that the index creation fails when presented with an unsorted BAM file
    @Test(expectedExceptions = SAMException.class)
    public void testBuildBamIndexFail() {
        final List<String> args = new ArrayList<String>();
        new SamFormatConverter().instanceMain(new String[]{
                "INPUT=" + INPUT_FILE,
                "OUTPUT=" + OUTPUT_UNSORTED_FILE});

        args.add("INPUT=" + OUTPUT_UNSORTED_FILE);
        args.add("OUTPUT=" + OUTPUT_INDEX_FILE);
        runPicardCommandLine(args);
    }

    @AfterTest
    public void cleanup() throws IOException {
        FileUtils.forceDeleteOnExit(OUTPUT_INDEX_FILE);
        FileUtils.forceDeleteOnExit(OUTPUT_SORTED_FILE);
        FileUtils.forceDeleteOnExit(OUTPUT_UNSORTED_FILE);
    }
}
