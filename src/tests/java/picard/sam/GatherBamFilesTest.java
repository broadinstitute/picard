package picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GatherBamFilesTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/GatherBamFiles");
    private static final File ORIG_BAM = new File(TEST_DATA_DIR, "orig.bam");
    private static final List<File> SPLIT_BAMS = Arrays.asList(
            new File(TEST_DATA_DIR, "indUnknownChrom.bam"),
            new File(TEST_DATA_DIR, "indchr1.bam"),
            new File(TEST_DATA_DIR, "indchr2.bam"),
            new File(TEST_DATA_DIR, "indchr3.bam"),
            new File(TEST_DATA_DIR, "indchr4.bam"),
            new File(TEST_DATA_DIR, "indchr5.bam"),
            new File(TEST_DATA_DIR, "indchr6.bam"),
            new File(TEST_DATA_DIR, "indchr7.bam"),
            new File(TEST_DATA_DIR, "indchr8.bam")
    );

    public String getCommandLineProgramName() {
        return GatherBamFiles.class.getSimpleName();
    }

    @Test
    public void testTheGathering() throws Exception {
        final File outputFile = File.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        outputFile.deleteOnExit();
        final List<String> args = new ArrayList<String>();
        for (final File splitBam : SPLIT_BAMS) {
            args.add("INPUT=" + splitBam.getAbsolutePath());
        }
        args.add("OUTPUT=" + outputFile);
        runPicardCommandLine(args);

        // TODO - Should switch over to using invocation via new PicardCommandLine() - BUT the test here is accessing class members directly.
        final CompareSAMs compareSAMs = new CompareSAMs();
        compareSAMs.samFiles = Arrays.asList(ORIG_BAM, outputFile);
        compareSAMs.doWork();
        Assert.assertTrue(compareSAMs.areEqual());
    }

    @Test
    public void sanityCheckTheGathering() throws Exception {
        final File outputFile = File.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        outputFile.deleteOnExit();
        final List<String> args = new ArrayList<String>();
        for (final File splitBam : SPLIT_BAMS) {
            args.add("INPUT=" + splitBam.getAbsolutePath());
        }
        args.add("OUTPUT=" + outputFile);
        runPicardCommandLine(args);

        // TODO - Should switch over to using invocation via new PicardCommandLine() - BUT the test here is accessing class members directly.
        final CompareSAMs compareSAMs = new CompareSAMs();
        compareSAMs.samFiles = Arrays.asList(ORIG_BAM, SPLIT_BAMS.get(0));
        compareSAMs.doWork();
        Assert.assertFalse(compareSAMs.areEqual());
    }
}
