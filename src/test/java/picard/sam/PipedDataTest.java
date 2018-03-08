package picard.sam;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

public class PipedDataTest {

    private String classPath = "\"" + System.getProperty("java.class.path") + "\" ";

    @Test
    public void testSortSam() {
        String sortCommand = "SortSam " +
                "I=/dev/stdin " +
                "O=/dev/null " +
                "SORT_ORDER=queryname";
        testPiping(sortCommand);
    }

    @Test
    public void testRevertSam() {
        String revertCommand = "RevertSam " +
                "I=/dev/stdin " +
                "O=/dev/null " +
                "SORT_ORDER=queryname " +
                "RESTORE_ORIGINAL_QUALITIES=false " +
                "REMOVE_DUPLICATE_INFORMATION=false " +
                "REMOVE_ALIGNMENT_INFORMATION=false " +
                "ATTRIBUTE_TO_CLEAR=[] SANITIZE=true";
        testPiping(revertCommand);
    }

    private void testPiping(String picardCommand) {
        String[] command = {
                "/bin/bash",
                "-c",
                "java -classpath " +
                        classPath +
                        "picard.cmdline.PicardCommandLine " +
                        "ViewSam " +
                        "I=testdata/picard/sam/test.bam " +
                        "ALIGNMENT_STATUS=All " +
                        "PF_STATUS=All " +
                        "| " +
                        "java -classpath " +
                        classPath +
                        "picard.cmdline.PicardCommandLine " +
                        picardCommand
        };
        try {
            ProcessBuilder processBuilder = new ProcessBuilder(command);
            processBuilder.inheritIO();

            Process process = processBuilder.start();
            Assert.assertEquals(process.waitFor(), 0);

        } catch (Exception e) {
            Assert.fail("Failed to pipe data from htsjdk to picard", e);
        }
    }
}
