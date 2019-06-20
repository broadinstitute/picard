package picard.arrays.illumina;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

public class Build37ExtendedIlluminaManifestTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");

    @DataProvider(name = "testBuild37ExtendedIlluminaManifestDataProvider")
    public Object[][] testBuild37ExtendedIlluminaManifestDataProvider() {
        return new Object[][]{
                { "HumanExome-12v1-1_A.1.3.extended.csv", "HumanExome-12v1-1_A.bpm", "4/23/2012", "Infinium HD Ultra", 8, "1.3", 6 },
                { "MEG_AllofUs_20002558X351448_A2.1.4.extended.csv", "MEG_AllofUs_20002558X351448_A2.bpm", "2/13/2019", "Infinium LCG", 5, "1.4", 1 }
        };
    }

    @Test(dataProvider = "testBuild37ExtendedIlluminaManifestDataProvider")
    public void testBuild37ExtendedIlluminaManifest(final String extendedManifestFilename, final String expectedDescriptorFilename,
                                                    final String expectedDateManufactured, final String expectedAssayFormat,
                                                    final int expectedLociCount, final String expectedExtendedManifestVersion,
                                                    final int expectedNumPass) throws IOException {
        final Build37ExtendedIlluminaManifest manifest = new Build37ExtendedIlluminaManifest(new File(TEST_DATA_DIR, extendedManifestFilename));
        Assert.assertEquals(manifest.getDescriptorFileName(), expectedDescriptorFilename);
        Assert.assertEquals(manifest.getDateManufactured(), expectedDateManufactured);
        Assert.assertEquals(manifest.getAssayFormat(), expectedAssayFormat);
        Assert.assertEquals(manifest.getLociCount(), expectedLociCount);
        Assert.assertEquals(manifest.getExtendedManifestVersion(), expectedExtendedManifestVersion);

        final Iterator<Build37ExtendedIlluminaManifestRecord> iterator = manifest.extendedIterator();
        int count = 0;
        int goodCount = 0;
        while (iterator.hasNext()) {
            count++;
            final Build37ExtendedIlluminaManifestRecord record = iterator.next();
            Assert.assertNotNull(record.getName());
            Assert.assertNotNull(record.getIlmnId());
            Assert.assertNotNull(record.getSnp());
            if (!record.isBad()) {
                Assert.assertNotNull(record.getB37Chr());
                Assert.assertNotNull(record.getB37Pos());
                Assert.assertNotNull(record.getAlleleA());
                Assert.assertNotNull(record.getAlleleB());
                Assert.assertNotNull(record.getRefAllele());
                goodCount++;
            }
        }
        Assert.assertEquals(count, expectedLociCount);
        Assert.assertEquals(goodCount, expectedNumPass);
    }
}
