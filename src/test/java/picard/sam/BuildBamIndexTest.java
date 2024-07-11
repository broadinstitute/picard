package picard.sam;

import htsjdk.beta.io.IOPathUtils;
import htsjdk.io.IOPath;
import htsjdk.samtools.CRAMCRAIIndexer;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.cram.CRAIEntry;
import htsjdk.samtools.cram.CRAIIndex;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.nio.PicardBucketUtils;
import picard.nio.PicardHtsPath;
import picard.nio.PicardIOUtils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class BuildBamIndexTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/indices/");
    private static final PicardHtsPath INPUT_UNSORTED_SAM = new PicardHtsPath(new File(TEST_DATA_DIR, "index_test.sam"));
    private static final File EXPECTED_BAI_FILE = new File(TEST_DATA_DIR, "index_test_b.bam.bai");

    private static final String CLOUD_TEST_DATA_DIR = "gs://hellbender/test/resources/picard/BuildBamIndex/";
    private static final String CLOUD_TEST_OUTPUT_DIR = "gs://hellbender-test-logs/staging/picard/test/BuildBamIndex/";
    // tsato: shouldn't we have a constructor for this?
    private static final PicardHtsPath INPUT_UNSORTED_SAM_CLOUD = new PicardHtsPath(CLOUD_TEST_DATA_DIR, "index_test.sam");

    // tsato: replace with variables defined in the other branches once they merge
    private static final PicardHtsPath SORTED_CRAM_CLOUD = new PicardHtsPath("gs://hellbender/test/resources/picard/BuildBamIndex/CEUTrio.HiSeq.WGS.b37.NA12878.20.21_n10000.cram");
    // tsato: this directory missing a crai...
    private static final PicardHtsPath SORTED_CRAM = new PicardHtsPath("/Users/tsato/workspace/picard/testdata/picard/test/CEUTrio.HiSeq.WGS.b37.NA12878.20.21_n10000.cram");
    private static final PicardHtsPath SORTED_CRAM2 = new PicardHtsPath("/Users/tsato/workspace/picard/testdata/picard/test/CEUTrio.HiSeq.WGS.b37.NA12878.20.21_n10000.cram");



    public String getCommandLineProgramName() { return BuildBamIndex.class.getSimpleName(); }

    @DataProvider(name = "buildBamIndexTestData")
    public Object[][] getBuildBamIndexTestData(){
            return new Object[][]{
                    {INPUT_UNSORTED_SAM, true},
                    {INPUT_UNSORTED_SAM, false},
                    {INPUT_UNSORTED_SAM_CLOUD, true},
                    {INPUT_UNSORTED_SAM_CLOUD, false}};
    }

    // tsato: must add cram test...

    // Test that the index file for a sorted BAM is created
    @Test(dataProvider = "buildBamIndexTestData")
    public void testBuildBamIndexOK(final PicardHtsPath inputUnsortedSam, final boolean specifyOutput) throws IOException {
        final boolean cloud = ! inputUnsortedSam.isLocalPath();
        final String prefix = "index_test_sorted";
        final PicardHtsPath sortedBAM = cloud ? PicardBucketUtils.getTempFilePath(CLOUD_TEST_OUTPUT_DIR, prefix,".bam") :
                PicardBucketUtils.getTempFilePath(null, prefix, ".bam");

        /* First sort, before indexing */ // tsato: do we need to do this dynamically?
        new SortSam().instanceMain(new String[]{
                "I=" + inputUnsortedSam,
                "O=" + sortedBAM,
                "SORT_ORDER=coordinate"});

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + sortedBAM);
        final PicardHtsPath indexOutput;

        if (specifyOutput) {
            indexOutput = cloud ? PicardBucketUtils.getTempFilePath(CLOUD_TEST_OUTPUT_DIR, prefix, ".bai") :
                    PicardBucketUtils.getTempFilePath(null, prefix,".bai");
            args.add("OUTPUT=" + indexOutput);
        } else {
            indexOutput = PicardHtsPath.replaceExtension(sortedBAM, ".bai", false);
            PicardIOUtils.deleteOnExit(indexOutput.toPath());
        }

        runPicardCommandLine(args);
        Assert.assertEquals(Files.readAllBytes(indexOutput.toPath()), Files.readAllBytes(EXPECTED_BAI_FILE.toPath()));
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

    @DataProvider(name = "cramTestData")
    public Object[][] getCramTestData(){
        // tsato: replace with variable
        final PicardHtsPath localRef = new PicardHtsPath("/Users/tsato/workspace/picard/testdata/picard/reference/human_g1k_v37.20.21.fasta");
        final PicardHtsPath cloudRef = new PicardHtsPath("gs://hellbender/test/resources/picard/references/human_g1k_v37.20.21.fasta");

        return new Object[][]{
                {SORTED_CRAM, localRef}};
                // {SORTED_CRAM, cloudRef}}; // ,
                // {SORTED_CRAM_CLOUD, localRef}, // tsato: these are too slow, disable for now.
                // {SORTED_CRAM_CLOUD, cloudRef},
    }

    @Test(dataProvider = "cramTestData")
    public void testCram(final PicardHtsPath cram, final PicardHtsPath reference) throws IOException {
        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + cram);
        args.add("REFERENCE_SEQUENCE=" + reference);

        final PicardHtsPath indexOutput;
        final String prefix = "BuildBamIndex_cram_test";

        final boolean specifyOutput = false; // for now
        if (specifyOutput) {
            indexOutput = PicardBucketUtils.getTempFilePath(CLOUD_TEST_OUTPUT_DIR, prefix, ".crai");
            args.add("OUTPUT=" + indexOutput);
        } else {
            indexOutput = PicardHtsPath.replaceExtension(cram, ".crai", false);
            // tsato: temporarily disable while investigating cram index created this way
            // PicardIOUtils.deleteOnExit(indexOutput.toPath());
        }

        runPicardCommandLine(args);

        final CRAIIndex craiIndex = CRAMCRAIIndexer.readIndex(Files.newInputStream(indexOutput.toPath()));
        final List<CRAIEntry> entries = craiIndex.getCRAIEntries();
        // Let's start with this

        // *** CORRECT ONES ***
        final CRAIIndex craiIndex2 = CRAMCRAIIndexer.readIndex(
                Files.newInputStream(new File("/Users/tsato/workspace/picard/testdata/picard/test/CEUTrio.HiSeq.WGS.b37.NA12878.20.21_n10000.cram.crai.samtools_save").toPath()));
        final List<CRAIEntry> entries2 = craiIndex2.getCRAIEntries();
        // *** END CORRECT ONES ***

        Assert.assertEquals(entries, entries2); // better test but uses a premade samtools index
        Assert.assertEquals(entries.size(), 2); // tsato: more crude test;

        Files.delete(indexOutput.toPath()); // tsato: temporary
    }

    // From https://github.com/samtools/htsjdk/issues/1084
    @Test
    public void testStdin() throws IOException {
        final InputStream inputStream = Files.newInputStream(Paths.get("/dev/stdin/"));
        inputStream.available();
    }
}
