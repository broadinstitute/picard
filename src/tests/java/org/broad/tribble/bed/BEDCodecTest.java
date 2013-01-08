/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broad.tribble.bed;

import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.TestUtils;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.bed.FullBEDFeature.Exon;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.index.linear.LinearIndex;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.awt.*;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

public class BEDCodecTest {

    @Test
    public void testSimpleDecode() {
        BEDCodec codec = new BEDCodec();

        BEDFeature feature;

        feature = codec.decode("chr1 1");
        Assert.assertEquals(feature.getChr(), "chr1");
        Assert.assertEquals(feature.getStart(), 2);
        Assert.assertEquals(feature.getEnd(), 2);

        feature = codec.decode("chr1 1 2");
        Assert.assertEquals(feature.getChr(), "chr1");
        Assert.assertEquals(feature.getStart(), 2);
        Assert.assertEquals(feature.getEnd(), 2);

        feature = codec.decode("chr1 1 3");
        Assert.assertEquals(feature.getChr(), "chr1");
        Assert.assertEquals(feature.getStart(), 2);
        Assert.assertEquals(feature.getEnd(), 3);
    }

    @Test
    public void testFullDecode() {
        BEDCodec codec = new BEDCodec();

        FullBEDFeature feature;
        List<Exon> exons;

        // Borrowed samples from Example: on http://genome.ucsc.edu/FAQ/FAQformat#format1

        feature = (FullBEDFeature) codec.decode("chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512");
        Assert.assertEquals(feature.getChr(), "chr22");
        Assert.assertEquals(feature.getStart(), 1001);
        Assert.assertEquals(feature.getEnd(), 5000);
        Assert.assertEquals(feature.getName(), "cloneA");
        Assert.assertEquals(feature.getScore(), 960f);
        Assert.assertEquals(feature.getStrand(), Strand.POSITIVE);
        Assert.assertEquals(feature.getColor(), new Color(0));

        exons = feature.getExons();
        Assert.assertEquals(exons.size(), 2);

        Assert.assertEquals(exons.get(0).getNumber(), 1);
        Assert.assertEquals(exons.get(0).start, 1001);
        Assert.assertEquals(exons.get(0).end, 1567);
        Assert.assertEquals(exons.get(0).getCdStart(), 1001);
        Assert.assertEquals(exons.get(0).getCdEnd(), 1567);
        Assert.assertEquals(exons.get(0).getCodingLength(), 567);

        Assert.assertEquals(exons.get(1).getNumber(), 2);
        Assert.assertEquals(exons.get(1).start, 4513);
        Assert.assertEquals(exons.get(1).end, 5000);
        Assert.assertEquals(exons.get(1).getCdStart(), 4513);
        Assert.assertEquals(exons.get(1).getCdEnd(), 5000);
        Assert.assertEquals(exons.get(1).getCodingLength(), 488);

        feature = (FullBEDFeature) codec.decode("chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601");
        Assert.assertEquals(feature.getChr(), "chr22");
        Assert.assertEquals(feature.getStart(), 2001);
        Assert.assertEquals(feature.getEnd(), 6000);
        Assert.assertEquals(feature.getName(), "cloneB");
        Assert.assertEquals(feature.getScore(), 900f);
        Assert.assertEquals(feature.getStrand(), Strand.NEGATIVE);
        Assert.assertEquals(feature.getColor(), new Color(0));

        exons = feature.getExons();
        Assert.assertEquals(exons.size(), 2);

        Assert.assertEquals(exons.get(0).getNumber(), 2);
        Assert.assertEquals(exons.get(0).start, 2001);
        Assert.assertEquals(exons.get(0).end, 2433);
        Assert.assertEquals(exons.get(0).getCdStart(), 2001);
        Assert.assertEquals(exons.get(0).getCdEnd(), 2433);
        Assert.assertEquals(exons.get(0).getCodingLength(), 433);

        Assert.assertEquals(exons.get(1).getNumber(), 1);
        Assert.assertEquals(exons.get(1).start, 5602);
        Assert.assertEquals(exons.get(1).end, 6000);
        Assert.assertEquals(exons.get(1).getCdStart(), 5602);
        Assert.assertEquals(exons.get(1).getCdEnd(), 6000);
        Assert.assertEquals(exons.get(1).getCodingLength(), 399);
    }

    @Test
    public void testDecodeBEDFile_good() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/NA12878.deletions.10kbp.het.gq99.hand_curated.hg19_fixed.bed";
        int expected_lines = 34;
        /*
        Line 0:
        1	25592413	25657872
        Line 3:
        1	152555536	152587611
        Line 28:
        14	73996607	74025282
        Remember tribble increments numbers by 1
         */

        BEDCodec codec = new BEDCodec();

        AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(filepath, codec, false);

        Iterable<Feature> iter = reader.iterator();
        int count = 0;
        for (Feature feat : iter) {
            Assert.assertTrue(feat.getChr().length() > 0);
            Assert.assertTrue(feat.getEnd() >= feat.getStart());

            if (count == 0) {
                Assert.assertEquals("1", feat.getChr());
                Assert.assertEquals(25592413 + 1, feat.getStart());
                Assert.assertEquals(25657872, feat.getEnd());
            }

            if (count == 3) {
                Assert.assertEquals("1", feat.getChr());
                Assert.assertEquals(152555536 + 1, feat.getStart());
                Assert.assertEquals(152587611, feat.getEnd());
            }

            if (count == 28) {
                Assert.assertEquals("14", feat.getChr());
                Assert.assertEquals(73996607 + 1, feat.getStart());
                Assert.assertEquals(74025282, feat.getEnd());
            }

            count += 1;
        }

        Assert.assertEquals(expected_lines, count);

        reader.close();

    }

    /**
     * Test reading a BED file which is malformed.
     *
     * @throws Exception
     */
    @Test(expectedExceptions = RuntimeException.class)
    public void testDecodeBEDFile_bad() throws Exception {
        //This file has an extra tab in the second to last line
        String filepath = TestUtils.DATA_DIR + "bed/NA12878.deletions.10kbp.het.gq99.hand_curated.hg19.bed";
        //The iterator implementation next() actually performs a get / read_next. The bad line is number 32,
        //so we actually will only get 31 lines before reading that line.
        int expected_count = 31;
        BEDCodec codec = new BEDCodec();

        AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(filepath, codec, false);

        Iterable<Feature> iter = reader.iterator();
        int count = 0;
        for (Feature feat : iter) {
            count += 1;
        }
        reader.close();
    }

    private void createIndex(File testFile, File idxFile) throws IOException {
        // Create an index if missing
        if (idxFile.exists()) {
            idxFile.delete();
        }
        LinearIndex idx = (LinearIndex) IndexFactory.createLinearIndex(testFile, new BEDCodec());

        LittleEndianOutputStream stream = null;
        try {
            stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(idxFile)));
            idx.write(stream);
        } finally {
            if (stream != null) {
                stream.close();
            }
        }

    }
}
