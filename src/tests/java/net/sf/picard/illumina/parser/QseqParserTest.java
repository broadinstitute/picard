package net.sf.picard.illumina.parser;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;

import static net.sf.picard.util.CollectionUtil.makeList;
import static net.sf.picard.illumina.parser.TestDataUtil.*;

public class QseqParserTest {

    //The actual testdata copied into source is in QseqTdUtil which provides test data (Td) for this an other tests

    @DataProvider(name = "qSeqSingleTile")
    public Object [][] qSeqOneEnds() {
        return new Object[][] {
            //using file s_1_1_0001.txt
            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)), 0, new int[]{76},getReadData(s_1_1_0001)},
                
            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)), 0, new int[]{38,38},
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{38,38}, 0) //get read data for lane 1 end 1 tile 1, expect it to be 76 bases but split QseqReadData.getBases() into two arrays of size 38
            },

            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)), 0, new int[]{0,38,38},
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{0,38,38}, 0)
            },

            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)), 26, new int[]{32,38,32},
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{32,38,32}, 26)
            },

            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)), 0, new int[]{32,38,32},
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{32,38,32}, 0)
            },

            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)), 0, new int[]{76,8,76},
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{76,8,76}, 0)
            },

            //using file s_1_1_0002_qseq.txt
            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0002)), 84, new int[]{76,8,76},
                getSplitOffsetReadData(s_1_1_0002, 76, new int[]{76,8,76}, 84)
            },

            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0002)), 0, new int[]{5,32,39},
                getSplitOffsetReadData(s_1_1_0002, 76,  new int[]{5,32,39}, 0)
            }
        };
    }

    @DataProvider(name = "qSeqMultiTile")
    public Object [][] qSeqMultiTile() {
        return new Object[][] {
            new Object[]{1, new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001, s_1_1_0002)), 0, new int[]{76},
                getTiledReadData(s_1_1, makeList(1,2))
            },
            //using file s_1_1_0001.txt
            new Object[]{1, new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001, s_1_1_0002)), 0, new int[]{38, 38},
                getSplitOffsetReadData(s_1_1, makeList(1,2), 76, new int[]{38,38}, 0)
            },

            //Using s_1_2_0001_qseq.txt as if it were s_1_1_0003_qseq.txt so we don't have to define more test data as static constants
            new Object[]{1, new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_1_0001, s_1_1_0002, s_1_2_0003)), 0, new int[]{38, 38},
                getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_2_0003), 76, new int[]{38,38}, 0)
            },

            new Object[]{1, new IlluminaFileMap(Arrays.asList(1,2,3), getQseqs(s_1_1_0001, s_1_1_0002, s_1_2_0003)), 0, new int[]{0,21,39,16},
                getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_2_0003), 76, new int[]{0,21,39, 16}, 0)
            }
        };
    }

    //SIMULATE THE FOLLOWING
    //read 1st end, middle end, last end, total 
    //1-4 Output arrays, different positions
    @Test(dataProvider="qSeqSingleTile")
    public void qSeqReadParserTest(final int lane, final IlluminaFileMap tilesToReadFiles, final int writeOffset, final int[] outputLengths, final Map<Integer, QseqReadData> testAgainst) {
        final QseqReadParser filesParser = new QseqReadParser(lane, tilesToReadFiles, writeOffset, outputLengths);

        int currentRead = 0;
        while(filesParser.hasNext()) {
            final QseqReadData qseqRead = new QseqReadData(outputLengths);
            filesParser.next(qseqRead);
            final QseqReadData testRead = testAgainst.get(currentRead);
            if(testRead != null) {
               testQSeqs(qseqRead, testRead);
            }
            currentRead++;
        }
    }

    @Test(dataProvider="qSeqMultiTile")
    public void qSeqReadParserTest_multiTile(final int lane, final IlluminaFileMap tilesToReadFiles, final int writeOffset, final int[] allWriteLengths, final Map<Integer, QseqReadData> testAgainst) {
        qSeqReadParserTest(lane, tilesToReadFiles, writeOffset, allWriteLengths, testAgainst);
    }

    @DataProvider(name="qSeqIntegrations")
    public Object[][] qSeqIntegrations() {
        return new Object[][] {
            //1 end 1 tile
            new Object[] {1, new int[]{76}, makeList(new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001))), getReadData(s_1_1_0001)},

            //1 end 2 tiles
            new Object[] {1, new int[]{38,38}, makeList(new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001, s_1_1_0002))),
                 getSplitOffsetReadData(s_1_1, makeList(1,2), 76, new int[]{38,38}, 0)},

            //2 ends 1 tile
            new Object[] {1, new int[]{152}, makeList(new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)), new IlluminaFileMap(makeList(1), getQseqs(s_1_2_0001))),
                combineReads(getReadData(s_1_1_0001), getReadData(s_1_2_0001))},

            //In this example tile 3 is tile 1 reused and therefore the resultant QSeqReadData will not read tile 3 but tile 1 TODO: Perhaps actually check the lane/tile number in the QSeqFileParser which would not allow us to do this
            //2 ends 3 tiles
            new Object[] {1, new int[]{76,76},
                makeList(new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_1_0001,s_1_1_0002, s_1_1_0003)), new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_2_0001, s_1_2_0002, s_1_2_0003))),
                combineReads(152, new int[]{76,76}, 0, getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_1_0003)), getTiledReadData(makeList(s_1_2_0001, s_1_2_0002, s_1_2_0003)))},

            //3 ends 2 tiles 3rd end is the same as the first
            new Object[] {1, new int[]{0, 38, 38, 76, 76},
                makeList(new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001,s_1_1_0002)), new IlluminaFileMap(makeList(1,2), getQseqs(s_1_2_0001, s_1_2_0002)), new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001,s_1_1_0002))),
               combineReads(228, new int[]{0, 38, 38, 76, 76}, 0, getTiledReadData(makeList(s_1_1_0001, s_1_1_0002)), getTiledReadData(makeList(s_1_2_0001, s_1_2_0002)), getTiledReadData(makeList(s_1_1_0001, s_1_1_0002)))},

            //4 ends 3 tiles
            new Object[] {1, new int[]{0, 38, 38, 76, 76, 0, 76},
                makeList(new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_1_0001,s_1_1_0002, s_1_1_0003)), new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_2_0001, s_1_2_0002, s_1_2_0003)), new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_1_0001,s_1_1_0002, s_1_1_0003)),  new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_2_0001, s_1_2_0002, s_1_2_0003))),
               combineReads(304, new int[]{0, 38, 38, 76, 76, 0, 76}, 0, getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_1_0003)), getTiledReadData(makeList(s_1_2_0001, s_1_2_0002, s_1_2_0003)), getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_1_0003)), getTiledReadData(makeList(s_1_2_0001, s_1_2_0002, s_1_2_0003)))},
        };
    }

    @Test(dataProvider="qSeqIntegrations")
    public void qSeqParserTest(final int lane, final int [] writeLengths, final List<IlluminaFileMap> readNumberToTile, final Map<Integer, QseqReadData> testAgainst) {
        final QseqParser qseqParser = new QseqParser(lane, writeLengths, readNumberToTile);

        int currentRead = 0;
        while(qseqParser.hasNext()) {
            final QseqReadData qseqRead = qseqParser.next();
            final QseqReadData testRead = testAgainst.get(currentRead);
            if(testRead != null) {
               testQSeqs(qseqRead, testRead);
            }
            currentRead++;
        }
    }

    public void testQSeqs(final QseqReadData actual, final QseqReadData expected){
        if(actual == expected) return;
        final byte [][] actualBases = actual.getBases();
        final byte [][] expectedBases = expected.getBases();
        final byte [][] actualQualities = actual.getQualities();
        final byte [][] expectedQualities = expected.getQualities();

        Assert.assertEquals(actualBases.length, expectedBases.length);
        Assert.assertEquals(actualQualities.length, expectedQualities.length);

        for(int i = 0; i < actualBases.length; i++) {
            for(int j = 0; j < actualBases[i].length; j++) {
                Assert.assertEquals(actualBases[i][j], expectedBases[i][j]);
                Assert.assertEquals(actualQualities[i][j], expectedQualities[i][j]);
            }
        }

        Assert.assertEquals(actual.isPf(), expected.isPf());
        Assert.assertEquals(actual.getXCoordinate(), expected.getXCoordinate());
        Assert.assertEquals(actual.getYCoordinate(), expected.getYCoordinate());
        Assert.assertEquals(actual.getLane(), expected.getLane());
        Assert.assertEquals(actual.getTile(), expected.getTile());
    }
}

