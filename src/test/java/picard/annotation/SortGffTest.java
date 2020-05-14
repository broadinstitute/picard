package picard.annotation;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;

public class SortGffTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/annotation/SortGff");
    public String getCommandLineProgramName() {
        return SortGff.class.getSimpleName();
    }

    @DataProvider(name = "testSortGffDataProvider")
    public Object[][] testSortGffDataProvider() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "basic.unsorted.gff3"), new File(TEST_DATA_DIR, "basic.sorted.gff3")},
                {new File(TEST_DATA_DIR, "child.before.parent.unsorted.gff3"), new File(TEST_DATA_DIR, "child.before.parent.sorted.gff3")}
        };
    }

    @Test(dataProvider = "testSortGffDataProvider")
    public void testBasicGff(final File inputGff, final File expectedOutputGff) throws IOException {
        final File outGff = File.createTempFile("testBasicGff", ".gff3");
        outGff.deleteOnExit();

        final String[] args = {
                "I=" + inputGff.getAbsolutePath(),
                "O=" + outGff.getAbsolutePath()
        };

        new SortGff().instanceMain(args);

        Assert.assertTrue(FileUtils.contentEquals(expectedOutputGff, outGff), "output file is different from expected file" + expectedOutputGff);
    }
}