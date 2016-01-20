package picard.util;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * Tests for IntervalListToBed
 */
public class IntervalListToBedTest {
    private final String TEST_DATA_DIR = "testdata/picard/util/";
    private final File INTERVAL_LIST = new File(TEST_DATA_DIR, "interval_list_to_bed_test.interval_list");
    private final File BED_FILE      = new File(TEST_DATA_DIR, "interval_list_to_bed_test.bed");

    @Test
    public void testConvertToBed() throws Exception {
        final IntervalListToBed program = new IntervalListToBed();
        final File tmp = File.createTempFile("interval_list_to_bed_test_output", ".bed");
        tmp.deleteOnExit();

        final String[] args = {
                "INPUT=" + INTERVAL_LIST.getAbsolutePath(),
                "OUTPUT=" + tmp.getAbsolutePath(),
                "SCORE=333"
        };
        program.instanceMain(args);

        final List<String> expected = IOUtil.slurpLines(BED_FILE);
        final List<String> actual   = IOUtil.slurpLines(tmp);

        // Make sure we got the same number of entries!
        Assert.assertEquals(actual.size(), expected.size());

        // Then make sure the entries are the same.
        for (int i=0; i<expected.size(); ++i) {
            Assert.assertEquals(actual.get(i), expected.get(i));
        }
    }

}
