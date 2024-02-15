package picard.fastq;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

public class BamToBfqTest {
    private final static File TEST_DATA_DIR = new File("testdata/picard/fastq");
    private final static File INPUT_BAM = new File("testdata/picard/sam/aligned_queryname_sorted.bam");

    @DataProvider(name = "inputs")
    public static Object[][] inputs() throws IOException {
        return new Object[][] {
                {INPUT_BAM, false, "bam_to_bfq_test"},
                {INPUT_BAM, true, "bam_to_bfq_paired_test"},
        };
    }

    @Test(dataProvider = "inputs")
    public void testBamToBfq(final File input, final boolean isPairedRun,
                     final String outputFilePrefix) throws IOException {
        final File analysisDir = IOUtil.createTempDir("BamToBfqTestDir").toFile();
        try {
            final String[] args = new String[] {
                    "INPUT=" + input.getAbsolutePath(),
                    "ANALYSIS_DIR=" + analysisDir.getAbsolutePath(),
                    "OUTPUT_FILE_PREFIX=" + outputFilePrefix,
                    "PAIRED_RUN=" + isPairedRun,
                    "READS_TO_ALIGN=8"
            };
            BamToBfq bamToBfq = new BamToBfq();
            Assert.assertEquals(bamToBfq.instanceMain(args), 0, "Can't process " + input.getAbsolutePath() + " correctly");

            final File output = new File(analysisDir, outputFilePrefix + ".0.1.bfq");
            final File expectedBFQ = new File(TEST_DATA_DIR, outputFilePrefix + ".0.1.bfq");

            Assert.assertEquals(Files.readAllBytes(output.toPath()), Files.readAllBytes(expectedBFQ.toPath()));
        } finally {
            IOUtil.recursiveDelete(analysisDir.toPath());
        }
    }
}
