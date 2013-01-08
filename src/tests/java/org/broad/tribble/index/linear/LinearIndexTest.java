/*
 * Copyright (c) 2009-2010 by The Broad Institute, Inc.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

package org.broad.tribble.index.linear;

import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.TestUtils;
import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.bed.BEDFeature;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class LinearIndexTest {
    private static final File RANDOM_FILE = new File("notMeaningful");

    private final static Block CHR1_B1 = new Block(1, 10);
    private final static Block CHR1_B2 = new Block(10, 20);
    private final static Block CHR1_B3 = new Block(20, 30);
    private final static Block CHR2_B1 = new Block(1, 100);
    private final static Block CHR2_B2 = new Block(100, 200);

    private LinearIndex idx;

    @BeforeTest
    public void setup() {
        idx = createTestIndex();
    }

    // chr1 (0, 10]
    // chr1 (10, 20]
    // chr1 (20, 30]
    // chr2 (0, 100]
    // chr2 (100, 200]
    private static LinearIndex createTestIndex() {
        LinearIndex.ChrIndex chr1 = new LinearIndex.ChrIndex("chr1", 10);
        chr1.addBlock(CHR1_B1);
        chr1.addBlock(CHR1_B2);
        chr1.addBlock(CHR1_B3);
        chr1.updateLongestFeature(1);

        LinearIndex.ChrIndex chr2 = new LinearIndex.ChrIndex("chr2", 100);
        chr2.addBlock(CHR2_B1);
        chr2.addBlock(CHR2_B2);
        chr2.updateLongestFeature(50);

        List<LinearIndex.ChrIndex> indices = Arrays.asList(chr1, chr2);
        return new LinearIndex(indices, RANDOM_FILE);
    }

    @Test()
    public void testBasicFeatures() {
        Assert.assertEquals(idx.getChrIndexClass(), LinearIndex.ChrIndex.class);
        Assert.assertEquals(idx.getType(), IndexFactory.IndexType.LINEAR.getHeaderValue());
        Assert.assertFalse(idx.hasFileSize());
        Assert.assertFalse(idx.hasTimestamp());
        Assert.assertFalse(idx.hasMD5());
        Assert.assertTrue(idx.isCurrentVersion());

        Assert.assertNotNull(idx.getSequenceNames());
        Assert.assertEquals(idx.getSequenceNames().size(), 2);
        Assert.assertTrue(idx.getSequenceNames().contains("chr1"));
        Assert.assertTrue(idx.getSequenceNames().contains("chr2"));
        Assert.assertTrue(idx.containsChromosome("chr1"));
        Assert.assertTrue(idx.containsChromosome("chr2"));
        Assert.assertFalse(idx.containsChromosome("chr3"));

        Assert.assertEquals(idx.getIndexedFile(), new File(RANDOM_FILE.getAbsolutePath()));

        Assert.assertNotNull(idx.getBlocks("chr1"));
        Assert.assertEquals(idx.getBlocks("chr1").size(), 3);

        Assert.assertNotNull(idx.getBlocks("chr2"));
        Assert.assertEquals(idx.getBlocks("chr2").size(), 2);
    }

    @Test()
    public void testEquals() {
        LinearIndex idx2 = createTestIndex();

        Assert.assertEquals(idx, idx, "Identical indices are equal");
        Assert.assertTrue(idx.equalsIgnoreProperties(idx), "Identical indices are equalIgnoreTimeStamp");
        Assert.assertTrue(idx.equalsIgnoreProperties(idx2), "Indices constructed the same are equalIgnoreTimeStamp");

        idx2.setTS(123456789);
        Assert.assertNotSame(idx, idx2, "Indices with different timestamps are not the same");
        Assert.assertTrue(idx.equalsIgnoreProperties(idx2), "Indices with different timestamps are equalIgnoreTimeStamp");
    }


    // chr1 (0, 10]
    // chr1 (10, 20]
    // chr1 (20, 30]
    // chr2 (0, 100]
    // chr2 (100, 200]
    //@Test()
    // TODO -- this is not a useful test as written -- the linear index always returns a single block since by
    // TODO -- definition they are contiguous and can be collapsed to a single block.
    public void testBasicQuery() {
        testQuery("chr1", 1, 1, CHR1_B1);
        testQuery("chr1", 1, 2, CHR1_B1);
        testQuery("chr1", 1, 9, CHR1_B1);
        testQuery("chr1", 10, 10, CHR1_B1);

        testQuery("chr1", 10, 11, CHR1_B1, CHR1_B2);
        testQuery("chr1", 11, 11, CHR1_B2);
        testQuery("chr1", 11, 12, CHR1_B2);
        testQuery("chr1", 11, 19, CHR1_B2);

        testQuery("chr1", 10, 19, CHR1_B1, CHR1_B2);
        testQuery("chr1", 10, 21, CHR1_B1, CHR1_B2, CHR1_B3);
        testQuery("chr1", 25, 30, CHR1_B3);
        testQuery("chr1", 35, 40);

        testQuery("chr2", 1, 1, CHR2_B1);
        testQuery("chr2", 100, 100, CHR2_B1);
        testQuery("chr2", 125, 125, CHR2_B1, CHR2_B2); // because of the 50 bp events
        testQuery("chr2", 151, 151, CHR2_B2); // because of the 50 bp events
        testQuery("chr2", 249, 249, CHR2_B2); // because of the 50 bp events
        testQuery("chr2", 251, 251); // just escaping the 50 bp longest event
    }

    private final void testQuery(String chr, int start, int stop, Block... expectedBlocksArray) {
        List<Block> qBlocks = idx.getBlocks(chr, start, stop);
        List<Block> eBlocks = Arrays.asList(expectedBlocksArray);

        Assert.assertEquals(qBlocks.size(), eBlocks.size(),
                String.format("Query %s:%d-%d returned %d blocks but we only expected %d.", chr, start, stop, qBlocks.size(), eBlocks.size()));
        for (int i = 0; i < qBlocks.size(); i++)
            Assert.assertEquals(qBlocks.get(i), eBlocks.get(i));
    }

    File fakeBed = new File(TestUtils.DATA_DIR + "fakeBed.bed");

    @Test
    public void oneEntryFirstChr() {
        BEDCodec code = new BEDCodec();
        Index index = IndexFactory.createLinearIndex(fakeBed, code);
        AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(fakeBed.getAbsolutePath(), code, index);

        try {
            CloseableTribbleIterator it = reader.iterator();
            int count = 0;
            while (it.hasNext()) {
                it.next();
                count++;
            }
            Assert.assertEquals(51, count);
        } catch (IOException e) {
            Assert.fail("Unable to get iterator due to " + e.getMessage());
        }
    }


    @Test
    /**
     *
     * chr2	1	200000000	LONG_FEATURE
     * ...
     * chr2	179098961	179380395	Hs.134602
     * chr2	179209546	179287210	Hs.620337
     * chr2	179266309	179266748	Hs.609465
     * chr2	179296428	179300012	Hs.623987
     * chr2	179302952	179303488	Hs.594545

     */
    public void testOverlappingFeatures() throws Exception {
        //chr2:179,222,066-179,262,059<- CONTAINS TTN

        Set<String> names = new HashSet<String>(Arrays.asList("Hs.134602", "Hs.620337", "Hs.609465", "Hs.623987",
                "Hs.594545", "LONG_FEATURE"));

        String bedFile = TestUtils.DATA_DIR + "bed/Unigene.sample.bed";
        String chr = "chr2";
        int start = 179266309;
        int end = 179303488;
        int expectedCount = 6;


        // Linear binned index
        LinearIndex.enableAdaptiveIndexing = false;
        int binSize = 1000;
        Index idx = IndexFactory.createLinearIndex(new File(bedFile), new BEDCodec(), binSize);

        FeatureReader<BEDFeature> bfr = AbstractFeatureReader.getFeatureReader(bedFile, new BEDCodec(), idx);
        CloseableTribbleIterator<BEDFeature> iter = bfr.query(chr, start, end);
        int countInterval = 0;
        while (iter.hasNext()) {
            BEDFeature feature = iter.next();
            Assert.assertTrue(feature.getEnd() >= start && feature.getStart() <= end);
            Assert.assertTrue(names.contains(feature.getName()));
            countInterval++;
        }

        Assert.assertEquals(countInterval, expectedCount);

        //Repeat with adaptive indexing
        LinearIndex.enableAdaptiveIndexing = true;
        idx = IndexFactory.createLinearIndex(new File(bedFile), new BEDCodec(), binSize);

        bfr = AbstractFeatureReader.getFeatureReader(bedFile, new BEDCodec(), idx);
        iter = bfr.query(chr, start, end);
        countInterval = 0;
        while (iter.hasNext()) {
            BEDFeature feature = iter.next();
            Assert.assertTrue(feature.getEnd() >= start && feature.getStart() <= end);
            Assert.assertTrue(names.contains(feature.getName()));
            countInterval++;
        }

        Assert.assertEquals(countInterval, expectedCount);


    }

}
