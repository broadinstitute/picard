package picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.FileExtensions;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.analysis.CollectAlignmentSummaryMetricsTest;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class AddCommentsToBamTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam");

    private static final File INPUT_FILE = new File(TEST_DATA_DIR, "aligned_queryname_sorted.bam");
    private static final File SAM_FILE = new File(CollectAlignmentSummaryMetricsTest.TEST_DATA_DIR, "summary_alignment_stats_test2.sam");

    private static final String[] commentList = new String[]{"test1", "test2", "test3"};

    public String getCommandLineProgramName() {
        return AddCommentsToBam.class.getSimpleName();
    }

    @Test
    public void testAddCommentsToBam() throws Exception {
        final File outputFile = File.createTempFile("addCommentsToBamTest.", BamFileIoUtils.BAM_FILE_EXTENSION);
        outputFile.deleteOnExit();
        runIt(INPUT_FILE, outputFile, commentList);

        final SAMFileHeader newHeader = SamReaderFactory.makeDefault().getFileHeader(outputFile);

        // The original comments are massaged when they're added to the header. Perform the same massaging here,
        // and then compare the lists
        final List<String> massagedComments = new LinkedList<>();
        for (final String comment : commentList) {
            massagedComments.add(SAMTextHeaderCodec.COMMENT_PREFIX + comment);
        }

        Assert.assertEquals(newHeader.getComments(), massagedComments);
        outputFile.delete();
    }

    @Test(expectedExceptions = PicardException.class)
    public void testUsingSam() throws Exception {
        final File outputFile = File.createTempFile("addCommentsToBamTest.samFile", FileExtensions.BAM);
        outputFile.deleteOnExit();
        runIt(SAM_FILE, outputFile, commentList);
        outputFile.delete();
        throw new IllegalStateException("We shouldn't be here!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUsingNewlines() throws Exception {
        final File outputFile = File.createTempFile("addCommentsToBamTest.mewLine", FileExtensions.BAM);
        outputFile.deleteOnExit();
        runIt(SAM_FILE, outputFile, new String[]{"this is\n a crazy\n test"});
        outputFile.delete();
        throw new IllegalStateException("We shouldn't be here!");
    }

    private void runIt(final File inputFile, final File outputFile, final String[] commentList) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "INPUT=" + inputFile.getAbsolutePath(),
                "OUTPUT=" + outputFile.getAbsolutePath()));
        for (final String comment : commentList) {
            args.add("COMMENT=" + comment);
        }
        runPicardCommandLine(args);
    }

}
