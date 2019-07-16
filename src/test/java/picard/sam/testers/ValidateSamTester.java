package picard.sam.testers;

import org.testng.Assert;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.ValidateSamFile;

import java.io.File;

/**
 * Created by ggrant on 6/12/14.
 */
public class ValidateSamTester extends CommandLineProgramTest {
    public String getCommandLineProgramName() {
        return ValidateSamFile.class.getSimpleName();
    }

    public void assertSamValid(final File samFile) {
        final int validateExitStatus = runPicardCommandLine(new String[]{"I=" + samFile.getAbsolutePath()});
        Assert.assertEquals(validateExitStatus, 0);
    }
}
