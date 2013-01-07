/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broad.tribble.index;

import org.broad.tribble.FeatureCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.bed.BEDCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Aug-23
 */
public class IndexFactoryTest {

    final File sortedBedFile = new File("test/data/bed/Unigene.sample.bed");
    final File unsortedBedFile = new File("test/data/bed/unsorted.bed");
    final File discontinuousFile = new File("test/data/bed/disconcontigs.bed");
    final FeatureCodec bedCodec = new BEDCodec();

    @Test
    public void testCreateLinearIndex() throws Exception {
        Index index = IndexFactory.createLinearIndex(sortedBedFile, bedCodec);
        String chr = "chr2";

        Assert.assertTrue(index.getSequenceNames().contains(chr));
        Assert.assertTrue(index.containsChromosome(chr));
        Assert.assertEquals(1, index.getSequenceNames().size());
        List<Block> blocks = index.getBlocks(chr, 1, 50);
        Assert.assertEquals(1, blocks.size());

        Block block = blocks.get(0);
        Assert.assertEquals(78, block.getSize());
    }

    @Test(expectedExceptions = TribbleException.MalformedFeatureFile.class, dataProvider = "indexFactoryProvider")
    public void testCreateIndexUnsorted(IndexFactory.IndexType type) throws Exception{
        Index index = IndexFactory.createIndex(unsortedBedFile, bedCodec, type);
    }

    @Test(expectedExceptions = TribbleException.MalformedFeatureFile.class, dataProvider = "indexFactoryProvider")
    public void testCreateIndexDiscontinuousContigs(IndexFactory.IndexType type) throws Exception{
        Index index = IndexFactory.createIndex(discontinuousFile, bedCodec, type);
    }

    @DataProvider(name = "indexFactoryProvider")
    public Object[][] getIndexFactoryTypes(){
        return new Object[][] {
                new Object[] { IndexFactory.IndexType.LINEAR },
                new Object[] { IndexFactory.IndexType.INTERVAL_TREE }
        };
    }

}
