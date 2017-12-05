package picard.sam.markduplicates;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.util.*;

public class MarkDuplicatesExceptionTest extends CommandLineProgramTest {

    protected static String TEST_BASE_NAME = null;
    protected static File TEST_DATA_DIR = null;

    @BeforeClass
    public void setUp() {
        TEST_BASE_NAME = "MD_IT";
        TEST_DATA_DIR = new File("testdata/picard/sam/MarkDuplicates/IntTest");
    }


    @Test(expectedExceptions = PicardException.class,
            expectedExceptionsMessageRegExp = "This program requires input that are either coordinate or query sorted. Found unsorted")
    public void testExceptionCoordinateOrQuerySortRequired() {
        final String input = "sort_exception_test.bam";
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final ArrayList<String> args = new ArrayList<>();
        args.add("INPUT=" + new File(TEST_DATA_DIR, input).getAbsolutePath());
        final File output = new File(outputDir, TEST_BASE_NAME + ".sam");
        args.add("OUTPUT=" + output.getAbsolutePath());
        final File metrics = new File(outputDir, TEST_BASE_NAME + ".integration_metrics");
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());

        Assert.assertEquals(runPicardCommandLine(args), 0);
    }


    @Override
    public String getCommandLineProgramName() { return MarkDuplicates.class.getSimpleName();
    }
}


