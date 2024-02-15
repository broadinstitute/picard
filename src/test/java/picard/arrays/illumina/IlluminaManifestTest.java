package picard.arrays.illumina;

import htsjdk.tribble.annotation.Strand;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

public class IlluminaManifestTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_ILLUMINA_MANIFEST_FILE = new File(TEST_DATA_DIR, "HumanExome-12v1-1_A.csv");
    private static final File TEST_BAD_ILLUMINA_MANIFEST_FILE = new File(TEST_DATA_DIR, "HumanExome-12v1-1_A.1.3.extended.csv");

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
            if (count == 1) {
                Assert.assertEquals(record.getIlmnId(), "exm-rs1000026-131_T_F_1990486741");
                Assert.assertEquals(record.getName(), "exm-rs1000026");
                Assert.assertEquals(record.getIlmnStrand(), IlluminaManifestRecord.IlluminaStrand.TOP);
                Assert.assertEquals(record.getSnp(), "[A/G]");
                Assert.assertTrue(record.isSnp());
                Assert.assertFalse(record.isIndel());
                Assert.assertFalse(record.isAmbiguous());
                Assert.assertEquals(record.getAddressAId(), "81790917");
                Assert.assertEquals(record.getAlleleAProbeSeq(), "TAAGTCAAAGGAAAACAAGTCAATAAATCCACTATCTATGGCTCCAAGGA");
                Assert.assertNull(record.getAddressBId());
                Assert.assertNull(record.getAlleleBProbeSeq());
                Assert.assertEquals(record.getGenomeBuild(), "37.1");
                Assert.assertEquals(record.getMajorGenomeBuild(), "37");
                Assert.assertEquals(record.getChr(), "21");
                Assert.assertEquals(record.getPosition(), 38934599);
                Assert.assertEquals(record.getPloidy(), "diploid");
                Assert.assertEquals(record.getSpecies(), "Homo sapiens");
                Assert.assertEquals(record.getSource(), "dbSNP");
                Assert.assertEquals(record.getSourceVersion(), "131");
                Assert.assertEquals(record.getSourceStrand(), IlluminaManifestRecord.IlluminaStrand.TOP);
                Assert.assertEquals(record.getSourceSeq(), "GAAAGAGCCATAAGTCAAAGGAAAACAAGTCAATAAATCCACTATCTATGGCTCCAAGGA[A/G]TAGAGGAAGCACCCAAAGTGATATTATTGTGAAACATTATTATTAATATGGGAAAGCCGC");
                Assert.assertEquals(record.getTopGenomicSeq(), "GAAAGAGCCATAAGTCAAAGGAAAACAAGTCAATAAATCCACTATCTATGGCTCCAAGGA[A/G]TAGAGGAAGCACCCAAAGTGATATTATTGTGAAACATTATTATTAATATGGGAAAGCCGC");
                Assert.assertEquals(record.getBeadSetId(), 662);
                Assert.assertEquals(record.getExpClusters(), "3");
                Assert.assertEquals(record.getRefStrand(), Strand.NEGATIVE);
                Assert.assertFalse(record.getIntensityOnly());
            }
        }
        Assert.assertEquals(count, 4);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testInvalidIlluminaManifest() throws IOException {
        new IlluminaManifest(TEST_BAD_ILLUMINA_MANIFEST_FILE);
    }
}
