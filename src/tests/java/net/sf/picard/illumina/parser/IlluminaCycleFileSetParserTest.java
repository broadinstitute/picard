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
    public void testSeekCnf() {
        final ReadConfiguration readConfiguration = new ReadConfiguration();
        readConfiguration.setFirstStart(1);
        readConfiguration.setFirstEnd(4);
        readConfiguration.assertValid();

        final CnfParser cnfParser = new CnfParser(readConfiguration, TEST_DATA_DIR, LANE, null);
        cnfParser.seekToTile(2);
        int numReads;
        int noiseIndex = 0;
        for (numReads = 0; cnfParser.hasNext(); ++numReads) {
            final IlluminaReadData read = new IlluminaReadData();
            read.setFirstEnd(new IlluminaEndData());
            cnfParser.next(read);

            if(noiseIndex != noiseData.length && noiseData[noiseIndex].index == numReads) {
                final CnfTestReads testReads = noiseData[noiseIndex];
                //The cycles are arranged in channels but the comparisons below are done on a per cycle basis
                final FourChannelIntensityData fcid = read.getFirstEnd().getNoise();
                short [] a = fcid.getA();
                short [] c = fcid.getC();
                short [] g = fcid.getG();
                short [] t = fcid.getT();

                for(int i = 0; i < a.length; i++) {
                    Assert.assertEquals(testReads.acgtNoiseReads[i][0], a[i]);
                    Assert.assertEquals(testReads.acgtNoiseReads[i][1], c[i]);
                    Assert.assertEquals(testReads.acgtNoiseReads[i][2], g[i]);
                    Assert.assertEquals(testReads.acgtNoiseReads[i][3], t[i]);
                }
            }
        }
        Assert.assertEquals(numReads, TILE_2_READS);
    }

    private static class CnfTestReads {
        public CnfTestReads(int index, short [][] acgtNoiseReads) {
            this.index = index;
            this.acgtNoiseReads = acgtNoiseReads;
        }
        public int index;
        public short[][] acgtNoiseReads;
    }

    public static CnfTestReads [] noiseData = new CnfTestReads[]{
            new CnfTestReads(0,   new short[][]{ new short[]{0,  9, 0, 0},
                                                 new short[]{0, 10, 0, 0},
                                                 new short[]{0,  0, 0, 0},
                                                 new short[]{0,  0, 0, 0}}),

            new CnfTestReads(19,  new short[][]{ new short[]{0,  21, 0, 11},
                                                 new short[]{0,  13, 0, 0},
                                                 new short[]{0,  0,  0, 0},
                                                 new short[]{0,  0,  0, 0}}),

            new CnfTestReads(39,  new short[][]{ new short[]{10,  26, 0, 11},
                                                 new short[]{0,   20, 0, 0},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   12, 0, 0}}),

            new CnfTestReads(74,  new short[][]{ new short[]{9,    9, 0, 11},
                                                 new short[]{0,   10, 0, 10},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   13, 0, 0}}),

            new CnfTestReads(90,  new short[][]{ new short[]{8,   11, 0, 11},
                                                 new short[]{0,   15, 0, 16},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   14, 0, 0}}),

            new CnfTestReads(98,  new short[][]{ new short[]{6,   10, 0, 11},
                                                 new short[]{0,   10, 0, 9},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   8, 0, 0}}),

            new CnfTestReads(99,  new short[][]{ new short[]{9,   9, 0, 10},
                                                 new short[]{0,   10, 0, 9},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   8, 0, 0}})
    };
}
