package picard.arrays.illumina;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

public class IlluminaManifestTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_ILLUMINA_MANIFEST_FILE = new File(TEST_DATA_DIR, "HumanExome-12v1-1_A.csv");
    private static final File TEST_BAD_ILLUMINA_MANIFEST_FILE = new File(TEST_DATA_DIR, "HumanExome-12v1-1_A.extended.csv");

    @Test
    public void tesIlluminaManifest() throws IOException {
        final IlluminaManifest manifest = new IlluminaManifest(TEST_ILLUMINA_MANIFEST_FILE);
        Assert.assertEquals(manifest.getDescriptorFileName(), "HumanExome-12v1-1_A.bpm");
        Assert.assertEquals(manifest.getDateManufactured(), "4/23/2012");
        Assert.assertEquals(manifest.getAssayFormat(), "Infinium HD Ultra");
        Assert.assertEquals(manifest.getLociCount(), 4);

        final Iterator<IlluminaManifestRecord> iterator = manifest.iterator();
        int count = 0;
        while (iterator.hasNext()) {
            count++;
            final IlluminaManifestRecord record = iterator.next();
            Assert.assertNotNull(record.getChr());
            Assert.assertTrue(record.getPosition() >= 0);
            Assert.assertNotNull(record.getAlleleAProbeSeq());
            Assert.assertNotNull(record.getSourceSeq());
        }
        Assert.assertEquals(count, 4);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testInvalidIlluminaManifest() throws IOException {
        new IlluminaManifest(TEST_BAD_ILLUMINA_MANIFEST_FILE);
    }
}
