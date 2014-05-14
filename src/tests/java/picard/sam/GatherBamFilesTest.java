package picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class GatherBamFilesTest {
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

    @Test
    public void testTheGathering() throws Exception {
        final GatherBamFiles gatherBamFiles = new GatherBamFiles();
        gatherBamFiles.INPUT = SPLIT_BAMS;
        gatherBamFiles.OUTPUT = File.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        gatherBamFiles.doWork();

        final CompareSAMs compareSAMs = new CompareSAMs();
        compareSAMs.samFiles = Arrays.asList(ORIG_BAM, gatherBamFiles.OUTPUT);
        compareSAMs.doWork();
        Assert.assertTrue(compareSAMs.areEqual());
    }

    @Test
    public void sanityCheckTheGathering() throws Exception {
        final GatherBamFiles gatherBamFiles = new GatherBamFiles();
        gatherBamFiles.INPUT = SPLIT_BAMS;
        gatherBamFiles.OUTPUT = File.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        gatherBamFiles.doWork();

        final CompareSAMs compareSAMs = new CompareSAMs();
        compareSAMs.samFiles = Arrays.asList(ORIG_BAM, SPLIT_BAMS.get(0));
        compareSAMs.doWork();
        Assert.assertFalse(compareSAMs.areEqual());
    }
}
