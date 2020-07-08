package picard.sam.markduplicates;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Paths;

public class CheckDuplicateMarkingTest {
    private static final String TEST_FILES_DIR = "testdata/picard/sam/CheckDuplicateMarking";

    @DataProvider(name = "checkDuplicateMarkingDataProvider")
    public Object[][] checkDuplicateMarkingDataProvider() {
        return new Object[][]{
                {"pass_queryname.sam", 0},
                {"pass_coordinate.sam", 0},
                {"fail_mate_queryname.sam", 1},
                {"fail_mate_coordinate.sam", 1},
                {"fail_supplementary_queryname_1.sam", 1},
                {"fail_supplementary_coordinate_1.sam", 1},
                {"fail_supplementary_queryname_2.sam", 1},
                {"fail_supplementary_coordinate_2.sam", 1},
                {"fail_secondary_queryname.sam", 1},
                {"fail_secondary_coordinate.sam", 1},
        };
    }

    @Test(dataProvider = "checkDuplicateMarkingDataProvider")
    public void testCheckDuplicateMarking(final String input, int expectedReturn) {
        final CheckDuplicateMarking cmdLine = new CheckDuplicateMarking();
        final String[] args = {"I=" + Paths.get(TEST_FILES_DIR, input)};
        Assert.assertEquals(cmdLine.instanceMain(args), expectedReturn);
    }
}