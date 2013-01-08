package org.broad.tribble.gelitext;

import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.TestUtils;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;


/**
 * @author aaron
 *         <p/>
 *         Class GeliTextTest
 *         <p/>
 *         test out the geli text source codec and feature
 */
public class GeliTextTest {
    public static final File testFile = new File(TestUtils.DATA_DIR + "testGeliText.txt");
    public static Index index;
    private FeatureReader<GeliTextFeature> source;

    // setup a new source before each class

    @BeforeSuite
    public void beforeTest() {
        index = IndexFactory.createLinearIndex(testFile, new GeliTextCodec());
        source = AbstractFeatureReader.getFeatureReader(testFile.getAbsolutePath(), new GeliTextCodec(), index);
    }

    @Test
    public void testReadAllLines() {
        // Query
        try {
            Iterator<GeliTextFeature> iter = source.query("22", 14438070, 14592250);
            int count = 0;
            while (iter.hasNext()) {
                GeliTextFeature feat = iter.next();
                count++;
            }
            Assert.assertEquals(count, 50);
        } catch (IOException e) {
            Assert.fail("failed to generate iterator from feature source");
        }
    }

    @Test
    public void testGetSubRegion() {
        // Query
        try {
            Iterator<GeliTextFeature> iter = source.query("22", 14438070, 14539060); // should be the first 41 records
            int count = 0;
            while (iter.hasNext()) {
                GeliTextFeature feat = iter.next();
                count++;
            }
            Assert.assertEquals(count, 41);
        } catch (IOException e) {
            Assert.fail("failed to generate iterator from feature source");
        }
    }

    @Test
    public void testFirstRecord() {
        // Query
        try {
            Iterator<GeliTextFeature> iter = source.query("22", 14438070, 14592250);
            int count = 0;

            GeliTextFeature feat = iter.next();
            // check the first records contents
            // 22 14438070 A   0 0     GG      33.2618 33.2618 0       0       0       0     0 0       0       33.2618 0       0
            Assert.assertTrue("22".equals(feat.getChr()));
            Assert.assertEquals(feat.getStart(), 14438070);
            Assert.assertEquals('A', feat.getRefBase());
            Assert.assertEquals(feat.getDepthOfCoverage(), 0.0, 0.0001);
            Assert.assertEquals(feat.getMaximumMappingQual(), 0.0, 0.0001);
            Assert.assertTrue(DiploidGenotype.GG.equals(feat.getGenotype()));
            Assert.assertEquals(feat.getDepthOfCoverage(), 0.0, 0.0001);
            Assert.assertEquals(feat.getLODBestToReference(), 33.2618, 0.0001);
            Assert.assertEquals(feat.getLODBestToNext(), 33.2618, 0.0001);
            for (int x = 0; x < feat.getLikelihoods().length; x++) {
                if (x == DiploidGenotype.GG.ordinal())
                    Assert.assertEquals(feat.getLikelihoods()[x], 33.2618, 0.0001);
                else
                    Assert.assertEquals(feat.getLikelihoods()[x], 0, 0.0001);
            }

        } catch (IOException e) {
            Assert.fail("failed to generate iterator from feature source");
        }
    }
}
