package picard.arrays.illumina;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

public class Build37ExtendedIlluminaManifestTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_EXTENDED_ILLUMINA_MANIFEST_FILE = new File(TEST_DATA_DIR, "HumanExome-12v1-1_A.extended.csv");

    @Test
    public void testBuild37ExtendedIlluminaManifest() throws IOException {
        final Build37ExtendedIlluminaManifest manifest = new Build37ExtendedIlluminaManifest(TEST_EXTENDED_ILLUMINA_MANIFEST_FILE);
        Assert.assertEquals(manifest.getDescriptorFileName(), "HumanExome-12v1-1_A.bpm");
        Assert.assertEquals(manifest.getDateManufactured(), "4/23/2012");
        Assert.assertEquals(manifest.getAssayFormat(), "Infinium HD Ultra");
        Assert.assertEquals(manifest.getLociCount(), 4);
        Assert.assertEquals(manifest.getExtendedManifestVersion(), "1.3");

        final Iterator<Build37ExtendedIlluminaManifestRecord> iterator = manifest.extendedIterator();
        int count = 0;
        int goodCount = 0;
        while (iterator.hasNext()) {
            count++;
            final Build37ExtendedIlluminaManifestRecord record = iterator.next();
            Assert.assertNotNull(record.getB37Chr());
            Assert.assertNotNull(record.getB37Pos());
            Assert.assertNotNull(record.getAlleleA());
            Assert.assertNotNull(record.getAlleleB());
            Assert.assertNotNull(record.getRefAllele());

            if (!record.isBad()) {
                goodCount++;
            }
        }
        Assert.assertEquals(count, 4);
        Assert.assertEquals(goodCount, 4);
    }
}
