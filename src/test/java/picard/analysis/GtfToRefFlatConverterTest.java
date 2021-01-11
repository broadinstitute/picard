package picard.analysis;
import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.GtfToRefFlatConverter;
import picard.cmdline.CommandLineProgramTest;
import java.io.File;
import java.io.IOException;

/**
 * Tests for methods in GtfToRefflat
 *
 * @author Sophie Crennan
 */

public class GtfToRefFlatConverterTest extends CommandLineProgramTest {

    private static final File TEST_DIR = new File("testdata/picard/sam/GtfToRefflat/");

    public String getCommandLineProgramName() {
        return GtfToRefFlatConverter.class.getSimpleName();
    }

    @Test
    public void testBasic() throws IOException {
        final File gtf = new File(TEST_DIR + "/Gtfs/", "basic.gtf");
        final File expectedRefflat = new File(TEST_DIR + "/Refflats/", "basic.refflat");

        File refflat = new GtfToRefFlatConverter(gtf.getAbsoluteFile()).getRefFlat();

        Assert.assertEquals(FileUtils.contentEquals(expectedRefflat, refflat), true);
    }

    @Test
    public void testDups() throws IOException {
        final File gtf = new File(TEST_DIR + "/Gtfs/", "dupsGtf.gtf");
        final File expectedRefflat = new File(TEST_DIR + "/Refflats/", "dupsRefflat.refflat");

        File refflat = new GtfToRefFlatConverter(gtf.getAbsoluteFile()).getRefFlat();

        Assert.assertEquals(FileUtils.contentEquals(expectedRefflat, refflat), true);
    }

    @Test
    public void testEnsembl() throws IOException {
        final File gtf = new File(TEST_DIR + "/Gtfs/", "ensembl.gtf");
        final File expectedRefflat = new File(TEST_DIR + "/Refflats/", "ensembl.refflat");

        File refflat = new GtfToRefFlatConverter(gtf.getAbsoluteFile()).getRefFlat();

        Assert.assertEquals(FileUtils.contentEquals(expectedRefflat, refflat), true);
    }

    @Test
    public void testEnsemblSplicedStops() throws IOException {
        final File gtf = new File(TEST_DIR + "/Gtfs/", "ensemblSplicedStops.gtf");
        final File expectedRefflat = new File(TEST_DIR + "/Refflats/", "ensemblSplicedStops.refflat");

        File refflat = new GtfToRefFlatConverter(gtf.getAbsoluteFile()).getRefFlat();

        Assert.assertEquals(FileUtils.contentEquals(expectedRefflat, refflat), true);
    }

    @Test
    public void testSplitStop() throws IOException {
        final File gtf = new File(TEST_DIR + "/Gtfs/", "splitStop.gtf");
        final File expectedRefflat = new File(TEST_DIR + "/Refflats/", "splitStop.refflat");

        File refflat = new GtfToRefFlatConverter(gtf.getAbsoluteFile()).getRefFlat();

        Assert.assertEquals(FileUtils.contentEquals(expectedRefflat, refflat), true);
    }

}


