/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
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

package org.broad.tribble.index.interval;

import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.TestUtils;
import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.bed.BEDFeature;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * User: jrobinso
 * Date: Mar 24, 2010
 */
public class IntervalTreeTest {

    static IntervalTree tree;

    @BeforeClass
    public static void setupTree() {
        tree = new IntervalTree();
        tree.insert(new Interval(0, 3, null));
        tree.insert(new Interval(5, 8, null));
        tree.insert(new Interval(6, 10, null));
        tree.insert(new Interval(8, 9, null));
        tree.insert(new Interval(15, 23, null));
        tree.insert(new Interval(16, 21, null));
        tree.insert(new Interval(17, 19, null));
        tree.insert(new Interval(19, 20, null));
        tree.insert(new Interval(25, 30, null));
        tree.insert(new Interval(26, 27, null));
    }

    @Test
    public void testSearch() {

        final Interval queryInterval = new Interval(1, 2);
        List<Interval> intervals = tree.findOverlapping(queryInterval);
        Assert.assertNotNull(intervals);

        for (Interval iv : intervals) {
            Assert.assertTrue(queryInterval.overlaps(iv));
        }
    }

    @Test
    public void testBed() throws Exception {
        String bedFile = TestUtils.DATA_DIR + "/index/chrY_Y4_small.bed";
        tree = new IntervalTree();
        Assert.assertTrue(tree.isValid());

        BufferedReader br = new BufferedReader(new FileReader(bedFile));
        String nextLine = "";
        while ((nextLine = br.readLine()) != null) {
            if (!(nextLine.startsWith("#") || nextLine.startsWith("track"))) {
                String[] tokens = nextLine.split("\t");
                if (tokens.length > 2) {
                    int start = Integer.parseInt(tokens[1]);
                    int end = Integer.parseInt(tokens[2]);
                    tree.insert(new Interval(start, end));
                }

            }
        }

//        List iv = (List) tree.findOverlapping(new Interval(2770226, 2770300));
        Interval searchInterval = new Interval(2782632, 2782732);
        List<Interval> iv = tree.findOverlapping(searchInterval);
        for (Interval i : iv) {
            Assert.assertTrue(i.overlaps(searchInterval));
        }

        br.close();

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

        String bedFile = TestUtils.DATA_DIR + "/bed/Unigene.sample.bed";
        String chr = "chr2";
        int start = 179266309;
        int end = 179303488 ;
        int expectedCount = 6;


        // Interval tree index
        int batchSize = 1;
        Index idx = IndexFactory.createIntervalIndex(new File(bedFile), new BEDCodec(), batchSize);

        FeatureReader<BEDFeature> bfr = AbstractFeatureReader.getFeatureReader(bedFile, new BEDCodec(), idx);
        CloseableTribbleIterator<BEDFeature>iter = bfr.query(chr, start, end);
        int countInterval = 0;
        while (iter.hasNext()) {
            BEDFeature feature = iter.next();
            Assert.assertTrue(feature.getEnd() >= start && feature.getStart() <= end);
            Assert.assertTrue(names.contains(feature.getName()));
            countInterval++;
        }

        Assert.assertEquals(countInterval, expectedCount);


    }


}

