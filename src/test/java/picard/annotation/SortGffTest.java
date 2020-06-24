package picard.annotation;


import htsjdk.samtools.util.IOUtil;
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
                {new File(TEST_DATA_DIR, "basic.unsorted.with.comments.and.directives.gff3"), new File(TEST_DATA_DIR, "basic.sorted.with.comments.and.directives.gff3")},
                {new File(TEST_DATA_DIR, "child.before.parent.belongs.after.parent.unsorted.gff3"), new File(TEST_DATA_DIR, "child.belongs.after.parent.sorted.gff3")},
                {new File(TEST_DATA_DIR, "child.after.parent.belongs.after.parent.unsorted.gff3"), new File(TEST_DATA_DIR, "child.belongs.after.parent.sorted.gff3")},
                {new File(TEST_DATA_DIR, "child.before.parent.belongs.before.parent.unsorted.gff3"), new File(TEST_DATA_DIR, "child.belongs.before.parent.sorted.gff3")},
                {new File(TEST_DATA_DIR, "child.after.parent.belongs.before.parent.unsorted.gff3"), new File(TEST_DATA_DIR, "child.belongs.before.parent.sorted.gff3")}
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

        IOUtil.assertFilesEqual(expectedOutputGff, outGff);
    }

    @DataProvider(name = "testSortByDictDataProvider")
    public Object[][] testSortByDictDataProvider() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "basic.unsorted.gff3"), new File(TEST_DATA_DIR, "standard.dict"), new File(TEST_DATA_DIR, "basic.sorted.gff3")},
                {new File(TEST_DATA_DIR, "basic.unsorted.gff3"), new File(TEST_DATA_DIR, "reverse.dict"), new File(TEST_DATA_DIR, "basic.sorted.reverse.dict.gff3")}
        };
    }

    @Test(dataProvider = "testSortByDictDataProvider")
    public void testSortByDictDataProvider(final File inputGff, final File dict, final File expectedOutputGff) throws IOException {
        final File outGff = File.createTempFile("testBasicGff", ".gff3");
        outGff.deleteOnExit();

        final String[] args = {
                "I=" + inputGff.getAbsolutePath(),
                "O=" + outGff.getAbsolutePath(),
                "SD=" + dict.getAbsolutePath()
        };

        new SortGff().instanceMain(args);

        IOUtil.assertFilesEqual(expectedOutputGff, outGff);
    }
}