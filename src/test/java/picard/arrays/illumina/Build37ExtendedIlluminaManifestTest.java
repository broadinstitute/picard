package picard.arrays.illumina;

import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

public class Build37ExtendedIlluminaManifestTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final String TEST_EXOME_EXT_MANIFEST_FILENAME = "HumanExome-12v1-1_A.2.0.extended.csv";

    @DataProvider(name = "testBuild37ExtendedIlluminaManifestDataProvider")
    public Object[][] testBuild37ExtendedIlluminaManifestDataProvider() {
        return new Object[][]{
                { TEST_EXOME_EXT_MANIFEST_FILENAME, "HumanExome-12v1-1_A.bpm", "4/23/2012", "Infinium HD Ultra", 8, "2.0", 7 },
                { "HumanExome-12v1-1_A.1.3.extended.csv", "HumanExome-12v1-1_A.bpm", "4/23/2012", "Infinium HD Ultra", 4, "1.3", 4 },
                { "MEG_AllofUs_20002558X351448_A2.1.4.extended.csv", "MEG_AllofUs_20002558X351448_A2.bpm", "2/13/2019", "Infinium LCG", 5, "1.4", 1 },
                { "GDA-8v1-0_A5.2.0.extended.csv", "GDA-8v1-0_A5.bpm", "8/28/2020", "Infinium LCG", 6, "2.0", 6 }
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
            if (!record.isFail()) {
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

    @Test
    public void testBuild37ExtendedIlluminaManifestContent() throws IOException {
        final Build37ExtendedIlluminaManifest manifest = new Build37ExtendedIlluminaManifest(new File(TEST_DATA_DIR, TEST_EXOME_EXT_MANIFEST_FILENAME));

        final Iterator<Build37ExtendedIlluminaManifestRecord> iterator = manifest.extendedIterator();
        int count = 0;
        int goodCount = 0;
        while (iterator.hasNext()) {
            count++;
            final Build37ExtendedIlluminaManifestRecord record = iterator.next();
            Assert.assertNotNull(record.getName());
            Assert.assertNotNull(record.getIlmnId());
            Assert.assertNotNull(record.getSnp());
            if (!record.isFail()) {
                Assert.assertNotNull(record.getB37Chr());
                Assert.assertNotNull(record.getB37Pos());
                Assert.assertNotNull(record.getAlleleA());
                Assert.assertNotNull(record.getAlleleB());
                Assert.assertNotNull(record.getRefAllele());
                goodCount++;
            }
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
                Assert.assertEquals(record.getB37Chr(), "21");
                Assert.assertEquals(record.getB37Pos().intValue(), 38934599);
                Assert.assertEquals(record.getSnpRefAllele(), "C");
                Assert.assertEquals(record.getRefAllele(), Allele.REF_C);
                Assert.assertEquals(record.getSnpAlleleA(), "T");
                Assert.assertEquals(record.getAlleleA(), Allele.ALT_T);
                Assert.assertEquals(record.getSnpAlleleB(), "C");
                Assert.assertEquals(record.getAlleleB(), Allele.REF_C);
                Assert.assertEquals(record.getRsId(), "rs1000026");
                Assert.assertFalse(record.isFail());
                Assert.assertFalse(record.isDupe());
            }
        }
        Assert.assertEquals(count, 8);
        Assert.assertEquals(goodCount, 7);
    }
}
