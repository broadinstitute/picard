package picard.vcf;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by lichtens on 7/11/17.
 */
public class VcfToIntervalListTest extends CommandLineProgramTest {
    private static final File TEST_RESOURCE_DIR = new File("testdata/picard/vcf");

    public String getCommandLineProgramName() {
        return VcfToIntervalList.class.getSimpleName();
    }

    @DataProvider(name="VcfToIntervalListData")
    public Object[][] getVcfToIntervalListData() {
        return new Object[][] {
                // 11 total, 10 defined unique intervals (two variants are adjacent), one is filtered in INFO and
                // one is filtered in FORMAT FT, but only INFO counts
                { new File(TEST_RESOURCE_DIR, "small_m2_more_variants.vcf"), false, 11 - 1 - 1 },

                // 11 total, 10 defined unique intervals (two variants are adjacent)
                { new File(TEST_RESOURCE_DIR, "small_m2_more_variants.vcf"), true, 11 - 1 }
        };
    }

    @Test(dataProvider="VcfToIntervalListData")
    public void testExcludingFiltered(
            final File inputFile,
            final boolean includeFiltered,
            final int expectedIntervalsSize) throws IOException
    {
        final File outputFile = File.createTempFile("vcftointervallist_", ".interval_list");
        outputFile.deleteOnExit();
        final List<String> arguments = new ArrayList<>();
        arguments.add("I=" + inputFile.getAbsolutePath());
        arguments.add("O=" + outputFile.getAbsolutePath());
        if (includeFiltered) {
            arguments.add(VcfToIntervalList.INCLUDE_FILTERED_SHORT_NAME + "=true");
        }
        runPicardCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<Interval> intervals = IntervalList.fromFile(outputFile).getIntervals();

        Assert.assertEquals(intervals.size(), expectedIntervalsSize);
    }
}
