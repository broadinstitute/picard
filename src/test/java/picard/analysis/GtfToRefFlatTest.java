package picard.analysis;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;

/**
 * Tests for GtfToRefFlat
 *
 * @author Sophie Crennan
 */

public class GtfToRefFlatTest extends CommandLineProgramTest {

    private static final File GTF_TEST_DIR = new File("testdata/picard/util/GtfToRefflat/Gtfs/");
    private static final File REFFLAT_TEST_DIR = new File("testdata/picard/util/GtfToRefflat/Refflats/");

    public String getCommandLineProgramName() {
        return GtfToRefFlat.class.getSimpleName();
    }

    @Test
    public void testBasic() throws IOException {
        final File gtf = new File(GTF_TEST_DIR, "basic.gtf");
        final File expectedRefFlat = new File(REFFLAT_TEST_DIR, "basic.refflat");

        GtfToRefFlat gtfToRefFlat = new GtfToRefFlat();
        gtfToRefFlat.GTF = gtf;
        gtfToRefFlat.doWork();
        Assert.assertEquals(FileUtils.contentEquals(expectedRefFlat, gtfToRefFlat.getRefFlat()), true);
    }

    @Test
    public void testDups() throws IOException {
        final File gtf = new File(GTF_TEST_DIR, "dupsGtf.gtf");
        final File expectedRefFlat = new File(REFFLAT_TEST_DIR, "dupsRefflat.refflat");

        GtfToRefFlat gtfToRefFlat = new GtfToRefFlat();
        gtfToRefFlat.GTF = gtf;
        gtfToRefFlat.doWork();
        Assert.assertEquals(FileUtils.contentEquals(expectedRefFlat, gtfToRefFlat.getRefFlat()), true);
    }

    @Test
    public void testEnsembl() throws IOException {
        final File gtf = new File(GTF_TEST_DIR, "ensembl.gtf");
        final File expectedRefFlat = new File(REFFLAT_TEST_DIR, "ensembl.refflat");

        GtfToRefFlat gtfToRefFlat = new GtfToRefFlat();
        gtfToRefFlat.GTF = gtf;
        gtfToRefFlat.doWork();
        Assert.assertEquals(FileUtils.contentEquals(expectedRefFlat, gtfToRefFlat.getRefFlat()), true);
    }

    @Test
    public void testEnsemblSplicedStops() throws IOException {
        final File gtf = new File(GTF_TEST_DIR, "ensemblSplicedStops.gtf");
        final File expectedRefFlat = new File(REFFLAT_TEST_DIR, "ensemblSplicedStops.refflat");

        GtfToRefFlat gtfToRefFlat = new GtfToRefFlat();
        gtfToRefFlat.GTF = gtf;
        gtfToRefFlat.doWork();
        Assert.assertEquals(FileUtils.contentEquals(expectedRefFlat, gtfToRefFlat.getRefFlat()), true);
    }

    @Test
    public void testSplitStop() throws IOException {
        final File gtf = new File(GTF_TEST_DIR, "splitStop.gtf");
        final File expectedRefFlat = new File(REFFLAT_TEST_DIR, "splitStop.refflat");

        GtfToRefFlat gtfToRefFlat = new GtfToRefFlat();
        gtfToRefFlat.GTF = gtf;
        gtfToRefFlat.doWork();
        Assert.assertEquals(FileUtils.contentEquals(expectedRefFlat, gtfToRefFlat.getRefFlat()), true);
    }

}