package picard.util;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * Tests for IntervalListToBed
 */
public class IntervalListToBedTest {
    private static final String TEST_DATA_DIR = "testdata/picard/util/";
    private final File INTERVAL_LIST = new File(TEST_DATA_DIR, "interval_list_to_bed_test.interval_list");
    private final File UNSORTED_INTERVAL_LIST = new File(TEST_DATA_DIR, "unsorted_interval_list_to_bed_test.interval_list");
    private final File BED_FILE = new File(TEST_DATA_DIR, "interval_list_to_bed_test.bed");
    private final File NO_SORT_BED_FILE = new File(TEST_DATA_DIR, "no_sort_interval_list_to_bed_test.bed");

    @DataProvider()
    Object[][] testConvertILToBedData() {
        return new Object[][]{
                {INTERVAL_LIST, BED_FILE, true},
                {INTERVAL_LIST, NO_SORT_BED_FILE, false},
                {UNSORTED_INTERVAL_LIST, BED_FILE, true},
        };
    }

    @Test(dataProvider = "testConvertILToBedData")
    public void testConvertToBed(final File intervalList, final File bedFile, final boolean sort) throws Exception {
        final IntervalListToBed program = new IntervalListToBed();
        final File tmp = File.createTempFile("interval_list_to_bed_test_output", ".bed");
        tmp.deleteOnExit();

        final String[] args = {
                "INPUT=" + intervalList.getAbsolutePath(),
                "OUTPUT=" + tmp.getAbsolutePath(),
                "SCORE=333",
                "SORT=" + sort
        };
        Assert.assertEquals(program.instanceMain(args), 0);

        final List<String> expected = IOUtil.slurpLines(bedFile);
        final List<String> actual = IOUtil.slurpLines(tmp);

        // Make sure we got the same number of entries!
        Assert.assertEquals(actual.size(), expected.size());

        // Then make sure the entries are the same.
        for (int i = 0; i < expected.size(); ++i) {
            Assert.assertEquals(actual.get(i), expected.get(i));
        }
    }
}
