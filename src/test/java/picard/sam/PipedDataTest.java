package picard.sam;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

public class PipedDataTest extends CommandLineProgramTest {


    @Override
    public String getCommandLineProgramName() {
        return SortSam.class.getSimpleName();
    }

    @Test
    public void testPipedData() {

        String classPath = "\"" + System.getProperty("java.class.path") + "\" ";

        String[] readCommand = {
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
            ProcessBuilder readProcess = new ProcessBuilder(readCommand);
            readProcess.inheritIO();

            Process read = readProcess.start();
            Assert.assertEquals(read.waitFor(), 0);

        } catch (Exception e) {
            Assert.fail("Failed to pipe data from htsjdk to picard", e);
        }

    }
}
