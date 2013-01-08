package org.broad.tribble.dbsnp;

import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.TestUtils;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;


/**
 * @author aaron
 *         <p/>
 *         Class OldDbSNPCodecTest
 *         <p/>
 *         tests out the basics of a dbsnp feature
 */
public class OldDbSNPCodecTest {

    public static final File testFile = new File(TestUtils.DATA_DIR + "basicDbSNP.dbsnp");
    public static Index index;
    private FeatureReader<OldDbSNPFeature> reader;

    // setup a new source before each class
    @BeforeTest
    public void beforeTest() {
        index = IndexFactory.createLinearIndex(testFile, new OldDbSNPCodec());
        reader = AbstractFeatureReader.getFeatureReader(testFile.getAbsolutePath(), new OldDbSNPCodec(), index);

    }

    @Test
    public void testReadAllLines() {
        // Query
        try {
            Iterator<OldDbSNPFeature> iter = reader.query("1", 0, 35000000);
            int count = 0;
            while (iter.hasNext()) {
                OldDbSNPFeature feat = iter.next();
                count++;
            }
            Assert.assertEquals(count, 42);
        } catch (IOException e) {
            Assert.fail("failed to generate iterator from feature source");
        }
    }

    @Test
    public void testReadAllLinesAlternateQueryString1() {
        // Query
        try {
            Iterator<OldDbSNPFeature> iter = reader.query("1", 433, 1937);
            int count = 0;
            while (iter.hasNext()) {
                OldDbSNPFeature feat = iter.next();
                count++;
            }
            Assert.assertEquals(count, 41);
        } catch (IOException e) {
            Assert.fail("failed to generate iterator from feature source");
        }
    }

    @Test
    public void testReadAllLinesAlternateQueryString2() {
        // Query
        try {
            Iterator<OldDbSNPFeature> iter = reader.query("1", 1936, 1937);
            int count = 0;
            while (iter.hasNext()) {
                OldDbSNPFeature feat = iter.next();
                count++;
            }
            Assert.assertEquals(count, 1);
        } catch (IOException e) {
            Assert.fail("failed to generate iterator from feature source");
        }
    }

    @Test
    public void testReturnedDBSNPEntry() {
        // Query
        try {
            Iterator<OldDbSNPFeature> iter = reader.query("1", 492, 492);
            if (!iter.hasNext()) Assert.fail("expected at least one entry");

            // check all the fields from the file
            OldDbSNPFeature feat = iter.next();
            Assert.assertTrue(feat.getChr().equals("1"));
            Assert.assertEquals(feat.getStart(), 492);
            Assert.assertEquals(feat.getEnd(), 492);
            Assert.assertTrue(feat.getRsID().equals("rs55998931"));
            Assert.assertEquals(feat.getScore(), 0);
            Assert.assertEquals(feat.getStrand(), Strand.POSITIVE);
            Assert.assertEquals(feat.getNCBIRefBase(), "C");
            Assert.assertEquals(feat.getUCSCRefBase(), "C");
            Assert.assertEquals(2, feat.getObserved().length);
            Assert.assertTrue("C".equals(feat.getObserved()[0]));
            Assert.assertTrue("T".equals(feat.getObserved()[1]));
            Assert.assertTrue("genomic".equals(feat.getMolType()));
            Assert.assertTrue("single".equals(feat.getVariantType()));
            Assert.assertTrue("unknown".equals(feat.getValidationStatus()));
            Assert.assertEquals(feat.getAvHet(), 0.0, 0.0001);
            Assert.assertEquals(feat.getAvHetSE(), 0.0, 0.0001);
            Assert.assertTrue("unknown".equals(feat.getFunction()));
            Assert.assertTrue("exact".equals(feat.getLocationType()));
            Assert.assertEquals(feat.getWeight(), 1);


        } catch (IOException e) {
            Assert.fail("failed to generate iterator from feature source");
        }
    }
}
