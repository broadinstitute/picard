package picard.sam;

import htsjdk.beta.io.IOPathUtils;
import htsjdk.io.IOPath;
import htsjdk.samtools.SAMException;
import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class BuildBamIndexTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/indices/");
    private static final PicardHtsPath INPUT_UNSORTED_SAM = new PicardHtsPath(new File(TEST_DATA_DIR, "index_test.sam").getPath());
    private static final File EXPECTED_BAI_FILE = new File(TEST_DATA_DIR, "index_test_b.bam.bai");

    public String getCommandLineProgramName() { return BuildBamIndex.class.getSimpleName(); }

    // Test that the index file for a sorted BAM is created
    @Test
    public void testBuildBamIndexOK() throws IOException {
        final IOPath sortedBAM = IOPathUtils.createTempPath("index_test_sorted", ".bam");
        // don't create the output index file, but mark it for deletion
        final Path indexOutput = sortedBAM.toPath().resolveSibling("index_test.bam.bai");
        indexOutput.toFile().deleteOnExit();

        /* First sort, before indexing */
        new SortSam().instanceMain(new String[]{
                "I=" + INPUT_UNSORTED_SAM,
                "O=" + sortedBAM,
                "SORT_ORDER=coordinate"});

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + sortedBAM);
        args.add("OUTPUT=" + indexOutput);
        runPicardCommandLine(args);
        Assert.assertEquals(FileUtils.readFileToByteArray(indexOutput.toFile()), FileUtils.readFileToByteArray(EXPECTED_BAI_FILE));
    }

    // Test that the index file for a sorted BAM is created in the right place if no output file is specified
    @Test
    public void testBuildBamIndexDefaultOutput() throws IOException {
        final IOPath sortedBAM = IOPathUtils.createTempPath("index_test_sorted", ".bam");
        /* First sort, before indexing */
        new SortSam().instanceMain(new String[]{
                "I=" + INPUT_UNSORTED_SAM,
                "O=" + sortedBAM,
                "SORT_ORDER=coordinate"});

        // don't create the output index file, but construct the expected name, and mark it for deletion
        final String expectedIndexFileName = sortedBAM.getURIString().replace(".bam", ".bai");
        final IOPath indexOutput = new PicardHtsPath(expectedIndexFileName);
        indexOutput.toPath().toFile().deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + sortedBAM);
        runPicardCommandLine(args);
        Assert.assertEquals(FileUtils.readFileToByteArray(indexOutput.toPath().toFile()),  FileUtils.readFileToByteArray(EXPECTED_BAI_FILE));
    }

    // Test that the index creation fails when presented with a SAM file
    @Test(expectedExceptions = SAMException.class)
    public void testBuildSamIndexFail() {
        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + INPUT_UNSORTED_SAM);
        runPicardCommandLine(args);
    }

    // Test that the index creation fails when presented with an unsorted BAM file
    @Test(expectedExceptions = SAMException.class)
    public void testBuildBamIndexFail() {
        final IOPath unsortedBAM = IOPathUtils.createTempPath("index_test_sorted", ".bam");
        new SamFormatConverter().instanceMain(new String[]{
                "INPUT=" + INPUT_UNSORTED_SAM,
                "OUTPUT=" + unsortedBAM});

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + unsortedBAM);
        runPicardCommandLine(args);
    }

}
