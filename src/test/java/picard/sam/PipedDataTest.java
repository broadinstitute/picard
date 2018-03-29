package picard.sam;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

public class PipedDataTest {

    private String classPath = "\"" + System.getProperty("java.class.path") + "\" ";

    @Test
    public void testSortSam() {
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
                        "SortSam " +
                        "I=/dev/stdin " +
                        "O=/dev/null " +
                        "SORT_ORDER=queryname"
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
