package picard.sam;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;

public class ValidateSamFileTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIR = "testdata/picard/sam/ValidateSamFile/";

    @Override
    public String getCommandLineProgramName() {
        return ValidateSamFile.class.getSimpleName();
    }

    @DataProvider
    public Object[][] samFiles() {
        return new Object[][] {
                {"nofile", ValidateSamFile.ReturnTypes.FAILED.value()},
                {"good/sorted-pair.sam", ValidateSamFile.ReturnTypes.SUCCESSFUL.value()},
                {"bad/unpaired-mate.sam", ValidateSamFile.ReturnTypes.ERRORS.value()},
                {"bad/missing-rg-info.sam", ValidateSamFile.ReturnTypes.ERRORS_WARNINGS.value()},
                {"bad/sorted-pair-missing-rg.sam", ValidateSamFile.ReturnTypes.WARNINGS.value()}
        };
    }

    @Test(dataProvider = "samFiles")
    public void test(final String samFileName, final int exitStatus) {
        final int validateExitStatus = runPicardCommandLine(new String[]{"I=" + new File(TEST_DATA_DIR + samFileName).getAbsolutePath()});
        Assert.assertEquals(validateExitStatus, exitStatus);
    }
}
