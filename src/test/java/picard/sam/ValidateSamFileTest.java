package picard.sam;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;

public class ValidateSamFileTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIR = "testdata/picard/sam/ValidateSamFile/";
    private static final File INDEX_DATA_DIR = new File("testdata/picard/indices/");
    private static final File INDEX_INPUT_FILE = new File(INDEX_DATA_DIR, "index_test.sam");

    @Override
    public String getCommandLineProgramName() {
        return ValidateSamFile.class.getSimpleName();
    }

    @DataProvider
    public Object[][] samFiles() {
        return new Object[][] {
                {"nofile", ValidateSamFile.ReturnTypes.FAILED.value()},
                {"good/sorted-pair.sam", ValidateSamFile.ReturnTypes.SUCCESSFUL.value()},
                {"good/sorted-pair-v1.6.sam", ValidateSamFile.ReturnTypes.SUCCESSFUL.value()},
                {"bad/unpaired-mate.sam", ValidateSamFile.ReturnTypes.ERRORS.value()},
                {"bad/missing-rg-info.sam", ValidateSamFile.ReturnTypes.ERRORS_WARNINGS.value()},
                {"bad/missing-rg-info-v1.6.sam", ValidateSamFile.ReturnTypes.ERRORS_WARNINGS.value()},
                {"bad/sorted-pair-missing-rg.sam", ValidateSamFile.ReturnTypes.WARNINGS.value()}
        };
    }

    @DataProvider
    public Object[][] indexFiles() {
        return new Object[][] {
                {"index_test_b.bam", ValidateSamFile.ReturnTypes.SUCCESSFUL.value()}, // BAM file with BAI index
                {"index_test_c.bam", ValidateSamFile.ReturnTypes.SUCCESSFUL.value()}  // BAM file with CSI index
        };
    }

    @Test(dataProvider = "samFiles")
    public void test(final String samFileName, final int exitStatus) {
        final int validateExitStatus = runPicardCommandLine(new String[]{"I=" + new File(TEST_DATA_DIR + samFileName).getAbsolutePath()});
        Assert.assertEquals(validateExitStatus, exitStatus);
    }

    @Test(dataProvider = "indexFiles")
    public void testIndex(final String bamFileName, final int exitStatus) throws IOException {
        File bamFile = new File(INDEX_DATA_DIR, bamFileName);
        new SortSam().instanceMain(new String[]{
                "I=" + INDEX_INPUT_FILE.getAbsolutePath(),
                "O=" + bamFile.getAbsolutePath(),
                "SORT_ORDER=coordinate"});

        final int validateExitStatus = runPicardCommandLine(
                new String[]{"I=" + bamFile.getAbsolutePath(),
                        "IGNORE=MISSING_PLATFORM_VALUE",
                        "IGNORE=MISMATCH_MATE_ALIGNMENT_START",
                        "IGNORE=MISSING_TAG_NM"});
        Assert.assertEquals(validateExitStatus, exitStatus);
        FileUtils.forceDeleteOnExit(bamFile);
    }
}
