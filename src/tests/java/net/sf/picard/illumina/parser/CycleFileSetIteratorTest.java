/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.illumina.parser;

import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.util.Iterator;

/**
 * @author alecw@broadinstitute.org
 */
public class CycleFileSetIteratorTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/rta");

    @Test
    public void testBasic() {
        final int lane = 1;
        final File laneDir = new File(TEST_DATA_DIR, "L00" + lane);
        final CycleFileSetIterator iterator = new CycleFileSetIterator(laneDir, lane, ClusterIntensityFileReader.FileType.cif, 4, null);
        final int expectedCycles = 4;
        final int expectedTiles = 2;
        final int firstTileNum = 1;
        checkExpectedTilesAndCycles(iterator, expectedCycles, expectedTiles, firstTileNum);
    }

    @Test
    public void testIgnoringCycles() {
        final int lane = 1;
        final File laneDir = new File(TEST_DATA_DIR, "L00" + lane);
        final CycleFileSetIterator iterator = new CycleFileSetIterator(laneDir, lane, ClusterIntensityFileReader.FileType.cif, 3, null);
        final int expectedCycles = 3;
        final int expectedTiles = 2;
        final int firstTileNum = 1;
        checkExpectedTilesAndCycles(iterator, expectedCycles, expectedTiles, firstTileNum);
    }

    @Test
    public void testSeekToTile() {
        final int lane = 1;
        final File laneDir = new File(TEST_DATA_DIR, "L00" + lane);
        final CycleFileSetIterator iterator = new CycleFileSetIterator(laneDir, lane, ClusterIntensityFileReader.FileType.cif, 4, null);
        final int expectedCycles = 4;
        final int expectedTiles = 1;
        iterator.seekToTile(2);
        final int firstTileNum = 2;
        checkExpectedTilesAndCycles(iterator, expectedCycles, expectedTiles, firstTileNum);
    }

    private void checkExpectedTilesAndCycles(final CycleFileSetIterator iterator, final int expectedCycles,
                                             final int expectedTiles, final int firstTileNum) {
        int numTiles;
        for (numTiles = 0; iterator.hasNext(); ++numTiles) {
            int numCycles;
            final Iterator<TiledIlluminaFile> innerIterator = iterator.next();
            for (numCycles = 0; innerIterator.hasNext(); ++numCycles) {
                final TiledIlluminaFile file = innerIterator.next();
                Assert.assertEquals(file.tile, numTiles + firstTileNum);
            }
            Assert.assertEquals(numCycles, expectedCycles);
        }
        Assert.assertEquals(numTiles, expectedTiles);
    }
}
