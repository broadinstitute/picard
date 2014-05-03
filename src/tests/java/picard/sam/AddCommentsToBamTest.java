package picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMTextHeaderCodec;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class AddCommentsToBamTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam");
    private static final File INPUT_FILE = new File(TEST_DATA_DIR, "aligned_queryname_sorted.bam");
    private static final File SAM_FILE = new File(TEST_DATA_DIR, "summary_alignment_stats_test2.sam");

    private static final List<String> commentList = new ArrayList(Arrays.asList("test1", "test2", "test3"));

    @Test
    public void testAddCommentsToBam() throws Exception {
        final AddCommentsToBam addCommentToBam = new AddCommentsToBam();
        addCommentToBam.INPUT = INPUT_FILE;
        addCommentToBam.OUTPUT = File.createTempFile("addCommentsToBamTest.", BamFileIoUtils.BAM_FILE_EXTENSION);
        addCommentToBam.COMMENT = commentList;
        addCommentToBam.doWork();
        final SAMFileHeader newHeader = new SAMFileReader(addCommentToBam.OUTPUT).getFileHeader();

        // The original comments are massaged when they're added to the header. Perform the same massaging here,
        // and then compare the lists
        final List<String> massagedComments = new LinkedList<String>();
        for (final String comment : commentList) {
            massagedComments.add(SAMTextHeaderCodec.COMMENT_PREFIX + comment);
        }

        Assert.assertEquals(newHeader.getComments(), massagedComments);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testUsingSam() throws Exception {
        final AddCommentsToBam addCommentToBam = new AddCommentsToBam();
        addCommentToBam.INPUT = SAM_FILE;
        addCommentToBam.OUTPUT = File.createTempFile("addCommentsToBamTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        addCommentToBam.COMMENT = commentList;
        addCommentToBam.doWork();
        throw new IllegalStateException("We shouldn't be here!");
    }

    @Test(expectedExceptions = PicardException.class)
    public void testUsingNewlines() throws Exception {
        final AddCommentsToBam addCommentToBam = new AddCommentsToBam();
        addCommentToBam.INPUT = INPUT_FILE;
        addCommentToBam.OUTPUT = File.createTempFile("addCommentsToBamTest.newLine.", BamFileIoUtils.BAM_FILE_EXTENSION);
        addCommentToBam.COMMENT = new ArrayList(Arrays.asList("this is\n a crazy\n test"));
        addCommentToBam.doWork();
        throw new IllegalStateException("We shouldn't be here!");
    }

}
