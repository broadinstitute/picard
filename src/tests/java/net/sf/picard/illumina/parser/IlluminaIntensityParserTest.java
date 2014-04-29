package net.sf.picard.illumina.parser;

import java.io.File;

import org.testng.Assert;
import org.testng.annotations.Test;

public class IlluminaIntensityParserTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/L001");
    private static final int LANE = 1;
    private static final int TILE_1_READS = 153010;
    private static final int TILE_2_READS = 149900;
    private static final int READS = TILE_1_READS + TILE_2_READS;
    private static final int [] CYCLES = new int[]{1,2,3,4};
    private static final int [] OUTPUT_LENGTHS = new int[]{4};

    private static CycleIlluminaFileMap makeFileMap(final String extension) {
        final CycleIlluminaFileMap fileMap = new CycleIlluminaFileMap();
        fileMap.put(1, new CycleFilesIterator(TEST_DATA_DIR, 1, 1, CYCLES, extension));
        fileMap.put(2, new CycleFilesIterator(TEST_DATA_DIR, 1, 2, CYCLES, extension));
        return fileMap;
    }

    private static CifParser makeCifParser(final String readStructure) {
        final OutputMapping outMap = new OutputMapping(new ReadStructure(readStructure));
        return new CifParser(TEST_DATA_DIR, LANE, makeFileMap(".cif"), outMap);
    }

    /**
     * Read noise files, configured as single-end, non-barcoded run with length=4.
     */
    @Test
    public void testSingleEnd() {
        final CifParser cifParser = makeCifParser("4T");

        int numReads;
        for (numReads = 0; cifParser.hasNext(); ++numReads) {
            final RawIntensityData rda = cifParser.next();
            for(final IntensityChannel iChan : IntensityChannel.values()) {
                Assert.assertEquals(rda.getRawIntensities()[0].getChannel(iChan).length, 4);
            }
        }

        Assert.assertEquals(numReads, READS);
    }

    /**
     * Read raw intensity files, configured as paired-end, barcoded run.
     */
    @Test
    public void testPairedEndWithBarcode() {

        final int [] outputLengths = new int[]{1, 2, 1};
        final CifParser cifParser = makeCifParser("1T2B1T");

        int numReads;
        for (numReads = 0; cifParser.hasNext(); ++numReads) {
            final RawIntensityData rda = cifParser.next();
            final FourChannelIntensityData [] fcids = rda.getRawIntensities();
            for(int i = 0; i <outputLengths.length; i++) {
                for(final IntensityChannel iChan : IntensityChannel.values()) {
                    Assert.assertEquals(fcids[i].getChannel(iChan).length, outputLengths[i]);
                }
            }
        }
        Assert.assertEquals(numReads, READS);
    }

    @Test
    public void testSeekCnf() {

        final CycleIlluminaFileMap fileMap = makeFileMap(".cnf");
        final OutputMapping outMap = new OutputMapping(new ReadStructure("4T"));
        final CnfParser cnfParser = new CnfParser(TEST_DATA_DIR, LANE, fileMap, outMap);
        cnfParser.seekToTile(2);
        int numReads;
        int noiseIndex = 0;
        for (numReads = 0; cnfParser.hasNext(); ++numReads) {
            final FourChannelIntensityData [] fcids = cnfParser.next().getNoise();

            if(noiseIndex != noiseData.length && noiseData[noiseIndex].index == numReads) {
                final CnfTestReads testReads = noiseData[noiseIndex];
                //The cycles are arranged in channels but the comparisons below are done on a per cycle basis
                final FourChannelIntensityData fcid = fcids[0];
                final short [] a = fcid.getA();
                final short [] c = fcid.getC();
                final short [] g = fcid.getG();
                final short [] t = fcid.getT();

                for(int i = 0; i < a.length; i++) {
                    Assert.assertEquals(a[i], testReads.acgtNoiseReads[i][0], "A values for read " + numReads + " on cycle " + (i+1) + " are not the same!");
                    Assert.assertEquals(c[i], testReads.acgtNoiseReads[i][1], "C values for read " + numReads + " on cycle " + (i+1) + " are not the same!");
                    Assert.assertEquals(g[i], testReads.acgtNoiseReads[i][2], "G values for read " + numReads + " on cycle " + (i+1) + " are not the same!");
                    Assert.assertEquals(t[i], testReads.acgtNoiseReads[i][3], "T values for read " + numReads + " on cycle " + (i+1) + " are not the same!");
                }
                ++noiseIndex;
            }
        }
        Assert.assertEquals(numReads, TILE_2_READS);
    }

    private static class CnfTestReads {
        public CnfTestReads(final int index, final short [][] acgtNoiseReads) {
            this.index = index;
            this.acgtNoiseReads = acgtNoiseReads;
        }
        public int index;
        public short[][] acgtNoiseReads;
    }

    public static CnfTestReads [] noiseData = new CnfTestReads[]{
            new CnfTestReads(0,   new short[][]{ new short[]{0,  9, 0, 0}, //cycle 1 (a,c,g,t)
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

            new CnfTestReads(74,  new short[][]{ new short[]{9,   25, 0, 10},
                                                 new short[]{0,   9,  0, 10},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   12, 0, 0}}),

            new CnfTestReads(90,  new short[][]{ new short[]{8,   11, 0, 11},
                                                 new short[]{0,   15, 0, 16},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   14, 0, 0}}),

            new CnfTestReads(98,  new short[][]{ new short[]{6,   10, 0, 11},
                                                 new short[]{0,   10, 0, 9},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   8, 0, 0}}),

            new CnfTestReads(99,  new short[][]{ new short[]{9,   9,  0, 10},
                                                 new short[]{10,  16, 0, 10},
                                                 new short[]{0,   0,  0, 0},
                                                 new short[]{0,   12, 0, 0}})
    };
}
