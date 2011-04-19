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

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Tests CnfParser and CifParser implementations of IlluminaCycleFileSetParserTest.
 *
 * @author alecw@broadinstitute.org
 */
public class IlluminaCycleFileSetParserTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/rta");
    private static final int LANE = 1;
    private static final int TILE_1_READS = 153010;
    private static final int TILE_2_READS = 149900;
    private static final int READS = TILE_1_READS + TILE_2_READS;

    /**
     * Read noise files, configured as single-end, non-barcoded run with length=4.
     */
    @Test
    public void testSingleEnd() {
        final ReadConfiguration readConfiguration = new ReadConfiguration();
        readConfiguration.setFirstStart(1);
        readConfiguration.setFirstEnd(4);
        readConfiguration.assertValid();

        final CnfParser cnfParser = new CnfParser(readConfiguration, TEST_DATA_DIR, LANE, null);

        int numReads;
        for (numReads = 0; cnfParser.hasNext(); ++numReads) {
            final IlluminaReadData read = new IlluminaReadData();
            read.setFirstEnd(new IlluminaEndData());
            cnfParser.next(read);
            Assert.assertEquals(read.getFirstEnd().getNoise().getA().length, readConfiguration.getFirstLength());
            Assert.assertEquals(read.getFirstEnd().getNoise().getC().length, readConfiguration.getFirstLength());
            Assert.assertEquals(read.getFirstEnd().getNoise().getG().length, readConfiguration.getFirstLength());
            Assert.assertEquals(read.getFirstEnd().getNoise().getT().length, readConfiguration.getFirstLength());
        }
        Assert.assertEquals(numReads, READS);
    }

    /**
     * Read raw intensity files, configured as paired-end, barcoded run.
     */
    @Test
    public void testPairedEndWithBarcode() {
        final ReadConfiguration readConfiguration = new ReadConfiguration();
        readConfiguration.setFirstStart(1);
        readConfiguration.setFirstEnd(1);
        readConfiguration.setBarcoded(true);
        readConfiguration.setBarcodeStart(2);
        readConfiguration.setBarcodeEnd(3);
        readConfiguration.setPairedEnd(true);
        readConfiguration.setSecondStart(4);
        readConfiguration.setSecondEnd(4);
        readConfiguration.assertValid();

        final CifParser cifParser = new CifParser(readConfiguration, TEST_DATA_DIR, LANE, null);

        int numReads;
        for (numReads = 0; cifParser.hasNext(); ++numReads) {
            final IlluminaReadData read = new IlluminaReadData();
            read.setFirstEnd(new IlluminaEndData());
            read.setSecondEnd(new IlluminaEndData());
            read.setBarcodeRead(new IlluminaEndData());
            cifParser.next(read);
            Assert.assertEquals(read.getFirstEnd().getRawIntensities().getA().length, readConfiguration.getFirstLength());
            Assert.assertEquals(read.getSecondEnd().getRawIntensities().getC().length, readConfiguration.getSecondLength());
            Assert.assertEquals(read.getBarcodeRead().getRawIntensities().getG().length, readConfiguration.getBarcodeLength());
        }
        Assert.assertEquals(numReads, READS);
    }

    @Test
    public void testSeekToTileAndCompareWithNse() {
        final ReadConfiguration readConfiguration = new ReadConfiguration();
        readConfiguration.setFirstStart(1);
        readConfiguration.setFirstEnd(4);
        readConfiguration.assertValid();

        final CnfParser cnfParser = new CnfParser(readConfiguration, TEST_DATA_DIR, LANE, null);
        cnfParser.seekToTile(2);
        final NseParser nseParser = new NseParser(readConfiguration, TEST_DATA_DIR, 2, null);
        nseParser.seekToTile(2);
        int numReads;
        for (numReads = 0; cnfParser.hasNext(); ++numReads) {
            final IlluminaReadData read = new IlluminaReadData();
            read.setFirstEnd(new IlluminaEndData());
            cnfParser.next(read);
            final IlluminaReadData nseParserRead = new IlluminaReadData();
            nseParserRead.setFirstEnd(new IlluminaEndData());
            if (nseParser.hasNext()) {
                // the nse file only has 100 reads in it.
                nseParser.next(nseParserRead);
                Assert.assertEquals(read.getFirstEnd().getNoise(), nseParserRead.getFirstEnd().getNoise());
            }
        }
        Assert.assertEquals(numReads, TILE_2_READS);
    }
}
