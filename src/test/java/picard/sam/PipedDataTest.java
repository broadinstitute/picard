package picard.sam;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.stream.StreamSupport;

import static picard.sam.MergeBamAlignmentTest.fasta;
import static picard.sam.MergeBamAlignmentTest.unmappedBam;

public class PipedDataTest {
    private static final String classPath = "\"" + System.getProperty("java.class.path") + "\" ";
    private static final String picardCommandlinePreamble = "java -classpath " + classPath + "picard.cmdline.PicardCommandLine ";
    private static final String testBam = "testdata/picard/sam/test.bam";

    /**
     * Creates the command line argument to be used for piped (/dev/stdin) input tests.
     *
     * @param inputSAM the path to the input SAM file.
     * @return commandline for ViewSam.
     */
    private String getViewSamPicardCommand(final String inputSAM){
        return picardCommandlinePreamble +
                "ViewSam " +
                "I=" + inputSAM + " " +
                "ALIGNMENT_STATUS=All " +
                "PF_STATUS=All ";
    }

    @DataProvider(name = "pipedInputTestData")
    public Object[][] getPipedInputTestData() throws IOException {
        // Note the trailing space in each string block.
        // TODO: port ArgumentBuilder from GATK so this would not be necessary.
        final File sortSamOutput = File.createTempFile("SortSam_pipe_test", "bam");
        sortSamOutput.deleteOnExit();
        final String[] sortSamCommand = {
                "/bin/bash",
                "-c",
                getViewSamPicardCommand(testBam) +
                "| " +
                picardCommandlinePreamble +
                "SortSam " +
                "SORT_ORDER=queryname " +
                "I=/dev/stdin " +
                "O=" + sortSamOutput.getAbsolutePath()
        };

        final File revertSamOutput = File.createTempFile("RevertSam_pipe_test", "bam");
        revertSamOutput.deleteOnExit();
        final String[] revertSamCommand = {
            "/bin/bash",
            "-c",
            getViewSamPicardCommand(testBam) +
            "| " +
            picardCommandlinePreamble +
            "RevertSam " +
            "I=/dev/stdin " +
            "O=" + revertSamOutput.getAbsolutePath()
        };


        // Only test the case where the "alignedBam" argument is /dev/stdin, which is the standard use case (bwa -> MergeBamAlignment).
        //
        // Note that the input aligned bam must be query-name sorted (as it would be if it's piped from bwa).
        // MergeBamAlignment fails if the coordinate sorted aligned bam is given, after trying to query-name sort dynamically and failing to find the reference header items,
        // seemingly due to the fact that the input is /dev/stdin.
        final File mergeBamAlignmentOutput = File.createTempFile("MergeBamAlignment_pipe_test", "bam");
        mergeBamAlignmentOutput.deleteOnExit();
        final String[] mergeBamAlignmentCommand = {
                "/bin/bash",
                "-c",
                getViewSamPicardCommand(MergeBamAlignmentTest.alignedQuerynameSortedBam.getAbsolutePath()) + // here we have to use "+", not ",".
                "| " +
                picardCommandlinePreamble +
                "MergeBamAlignment " +
                "ALIGNED_BAM=/dev/stdin " +
                "UNMAPPED=" + unmappedBam + " " +
                "REFERENCE_SEQUENCE=" + fasta + " " +
                "OUTPUT=" + mergeBamAlignmentOutput.getAbsolutePath()
        };

        return new Object[][]{
                {sortSamCommand, sortSamOutput},
                {revertSamCommand, revertSamOutput},
                {mergeBamAlignmentCommand, mergeBamAlignmentOutput}
        };
    }

    @Test(dataProvider = "pipedInputTestData")
    public void testPipedInput(final String[] command, final File tmpOutput) {
        try {
            ProcessBuilder processBuilder = new ProcessBuilder(command);
            processBuilder.inheritIO();

            Process process = processBuilder.start();
            Assert.assertEquals(process.waitFor(), 0);

            final SamReader reader = SamReaderFactory.makeDefault().open(tmpOutput);
            final long numOutputReads = StreamSupport.stream(reader.spliterator(), false).count();
            Assert.assertTrue(numOutputReads > 0);
            Assert.assertTrue(reader.getFileHeader().getReadGroups().size() > 0);
        } catch (Exception e) {
            Assert.fail("Failed to pipe data to picard. The error was: ", e);
        }
    }
}
