package picard.sam;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

import static picard.sam.MergeBamAlignmentTest.fasta;
import static picard.sam.MergeBamAlignmentTest.unmappedBam;

public class PipedDataTest {
    private final String classPath = "\"" + System.getProperty("java.class.path") + "\" ";
    private final String picardCommandlinePreamble = "java -classpath " + classPath + "picard.cmdline.PicardCommandLine ";

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
        final String[] sortSamCommand = {
                "/bin/bash",
                "-c",
                getViewSamPicardCommand("testdata/picard/sam/test.bam") +
                "| " +
                picardCommandlinePreamble +
                "SortSam " +
                "I=/dev/stdin " +
                "O=/dev/null " +
                "SORT_ORDER=queryname"
        };

        final File temp = File.createTempFile("test", "sam");
        temp.deleteOnExit();
        final String[] revertSamCommand = {
            "/bin/bash",
            "-c",
            getViewSamPicardCommand("testdata/picard/sam/test.bam") +
            "| " +
            picardCommandlinePreamble +
            "RevertSam " +
            "I=/dev/stdin " +
            "O=/dev/null "
        };

        // Only test the case where the "alignedBam" argument is /dev/stdin, which is the standard use case (bwa -> mergebamalignment)
        // Make sure that the input aligned bam is query-name sorted (as it would be if it's piped from bwa)
        // MergeBamAlignment fails if the coordinate sorted aligned bam is given, after trying to query-name sort dynamically and failing to find the reference header items,
        // due seemingly to the fact that the input is /dev/stdin.
        final String[] mergeBamAlignmentCommand = {
                "/bin/bash",
                "-c",
                getViewSamPicardCommand(MergeBamAlignmentTest.alignedQuerynameSortedBam.getAbsolutePath()) + // here we have to use "+", not ",". Best to extract a method "piped input argument builder" or something
                "| " +
                picardCommandlinePreamble +
                "MergeBamAlignment " +
                "ALIGNED_BAM=/dev/stdin " +
                "UNMAPPED=" + unmappedBam + " " +
                "REFERENCE_SEQUENCE=" + fasta + " " +
                "OUTPUT=/dev/null "
        };

        return new Object[][]{
                {sortSamCommand},
                {revertSamCommand},
                {mergeBamAlignmentCommand}
        };
    }

    @Test(dataProvider = "pipedInputTestData")
    public void testPipedInput(final String[] command) {
        try {
            ProcessBuilder processBuilder = new ProcessBuilder(command);
            processBuilder.inheritIO();

            Process process = processBuilder.start();
            Assert.assertEquals(process.waitFor(), 0);
        } catch (Exception e) {
            Assert.fail("Failed to pipe data to picard. The error was: ", e);
        }
    }
}
