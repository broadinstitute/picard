package net.sf.picard.illumina.parser;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;

import static net.sf.picard.util.CollectionUtil.makeList;
import static net.sf.picard.illumina.parser.QSeqTdUtil.*;
import static net.sf.picard.illumina.parser.OutputMapping.TwoDIndex;

public class QseqParserTest {

    //The actual testdata copied into source is in QseqTdUtil which provides test data (Td) for this an other tests

    @DataProvider(name = "qSeqSingleTile")
    public Object [][] qSeqOneEnds() {
        return new Object[][] {
            //using file s_1_1_0001.txt
            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)),
                new Range[]{new Range(0,75)},
                new TwoDIndex[]{new TwoDIndex(0,0)},
                new OutputMapping(new ReadStructure("76T")),
                getReadData(s_1_1_0001)},
                
            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)),
                new Range[]{new Range(0,37), new Range(38,75)},
                new TwoDIndex[]{new TwoDIndex(0,0), new TwoDIndex(1,0)},
                new OutputMapping(new ReadStructure("38T38T")),
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{38,38}, 0) //get read data for lane 1 end 1 tile 1, expect it to be 76 bases but split QseqReadData.getBases() into two arrays of size 38
            },

            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)),
                new Range[]{ new Range(0, 5),   new Range(6,  43),  new Range(44, 75)},
                new TwoDIndex[]{ new TwoDIndex(0, 26), new TwoDIndex(1, 0),  new TwoDIndex(2, 0)},
                new OutputMapping(new ReadStructure("32T38B32T")),
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{32,38,32}, 26)
            },


            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)),
                new Range[]{ new Range(0,  31), new Range(32, 69), new Range(70, 75)},
                new TwoDIndex[]{ new TwoDIndex(0,  0), new TwoDIndex(1, 0), new TwoDIndex(2, 0)},
                new OutputMapping(new ReadStructure("32T38B32T")),
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{32,38,32}, 0)
            },

            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)),
                new Range[]{ new Range(0,  75)},
                new TwoDIndex[]{ new TwoDIndex(0,  0)},
                new OutputMapping(new ReadStructure("76B8T76B")),
                getSplitOffsetReadData(s_1_1_0001, 76, new int[]{76,8,76}, 0)
            },

            //using file s_1_1_0002_qseq.txt
            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0002)),
                new Range[]{ new Range(0,  75)},
                new TwoDIndex[]{ new TwoDIndex(2, 0)},
                new OutputMapping(new ReadStructure("76B8T76B")),
                getSplitOffsetReadData(s_1_1_0002, 76, new int[]{76,8,76}, 84)
            },

            new Object[]{1, new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0002)),
                new Range[]{ new Range(0,  4), new Range(5, 36), new Range(37, 75)},
                new TwoDIndex[]{ new TwoDIndex(0,  0), new TwoDIndex(1, 0), new TwoDIndex(2, 0)},
                new OutputMapping(new ReadStructure("5T32B39T")),
                getSplitOffsetReadData(s_1_1_0002, 76,  new int[]{5,32,39}, 0)
            }
        };
    }

    @DataProvider(name = "qSeqMultiTile")
    public Object [][] qSeqMultiTile() {
        return new Object[][] {
            new Object[]{1, new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001, s_1_1_0002)),
                new Range[]{ new Range(0,75)},
                new TwoDIndex[]{ new TwoDIndex(0,0)},
                new OutputMapping(new ReadStructure("76T")),
                getTiledReadData(s_1_1, makeList(1,2))
            },
            //using file s_1_1_0001.txt
            new Object[]{1, new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001, s_1_1_0002)),
                new Range[]{ new Range(0,37), new Range(38,75)},
                new TwoDIndex[]{ new TwoDIndex(0,0), new TwoDIndex(1, 0)},
                new OutputMapping(new ReadStructure("38T38T")),
                getSplitOffsetReadData(s_1_1, makeList(1,2), 76, new int[]{38,38}, 0)
            },

            //Using s_1_2_0001_qseq.txt as if it were s_1_1_0003_qseq.txt so we don't have to define more test data as static constants
            new Object[]{1, new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_1_0001, s_1_1_0002, s_1_2_0003)),
                new Range[]{ new Range(0,37), new Range(38,75)},
                new TwoDIndex[]{ new TwoDIndex(0,0), new TwoDIndex(1, 0)},
                new OutputMapping(new ReadStructure("38T38T")),
                getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_2_0003), 76, new int[]{38,38}, 0)
            },

            new Object[]{1, new IlluminaFileMap(Arrays.asList(1,2,3), getQseqs(s_1_1_0001, s_1_1_0002, s_1_2_0003)),
                new Range[]{ new Range(0,20), new Range(21,59), new Range(60, 75)},
                new TwoDIndex[]{ new TwoDIndex(0,0), new TwoDIndex(1,0), new TwoDIndex(2,0)},
                new OutputMapping(new ReadStructure("21B39T16T")),
                getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_2_0003), 76, new int[]{21,39, 16}, 0)
            }
        };
    }

    //SIMULATE THE FOLLOWING
    //read 1st end, middle end, last end, total 
    //1-4 Output arrays, different positions
    @Test(dataProvider="qSeqSingleTile")
    public void qSeqReadParserTest(final int lane, final IlluminaFileMap tilesToReadFiles, final Range [] srcRanges, final TwoDIndex[] dstRanges, final OutputMapping om, final Map<Integer, QseqReadData> testAgainst) {
        final QseqReadParser filesParser = new QseqReadParser(lane, tilesToReadFiles, srcRanges, dstRanges, om);
        final Map<Integer, QseqReadData> filteredTestAgainst = filterAllSkips(testAgainst, om);

        final int [] outputLengths = om.getOutputReadLengths();
        int currentRead = 0;
        while(filesParser.hasNext()) {
            final QseqReadData qseqRead = new QseqReadData(outputLengths);
            filesParser.next(qseqRead);
            final QseqReadData testRead = filteredTestAgainst.get(currentRead);
            if(testRead != null) {
               testQSeqs(qseqRead, testRead);
            }
            currentRead++;
        }
    }

    @Test(dataProvider="qSeqMultiTile")
    public void qSeqReadParserTest_multiTile(final int lane, final IlluminaFileMap tilesToReadFiles, final Range [] srcRanges, final TwoDIndex[] dstRanges, final OutputMapping om, final Map<Integer, QseqReadData> testAgainst) {
        qSeqReadParserTest(lane, tilesToReadFiles, srcRanges, dstRanges, om, testAgainst);
    }

    @DataProvider(name="qSeqIntegrations")
    public Object[][] qSeqIntegrations() {
        return new Object[][] {
            //1 end 1 tile
            new Object[] {1, "76T", makeList(new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001))), getReadData(s_1_1_0001)},

            //1 end 2 tiles
            new Object[] {1, "38B38B", makeList(new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001, s_1_1_0002))),
                 getSplitOffsetReadData(s_1_1, makeList(1,2), 76, new int[]{38,38}, 0)},

            //2 ends 1 tile
            new Object[] {1, "152T", makeList(new IlluminaFileMap(makeList(1), getQseqs(s_1_1_0001)), new IlluminaFileMap(makeList(1), getQseqs(s_1_2_0001))),
                combineReads(getReadData(s_1_1_0001), getReadData(s_1_2_0001))},

            //In this example tile 3 is tile 1 reused and therefore the resultant QSeqReadData will not read tile 3 but tile 1
            //2 ends 3 tiles
            new Object[] {1, "76T76T",
                makeList(new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_1_0001,s_1_1_0002, s_1_1_0003)), new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_2_0001, s_1_2_0002, s_1_2_0003))),
                combineReads(152, new int[]{76,76}, 0, getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_1_0003)), getTiledReadData(makeList(s_1_2_0001, s_1_2_0002, s_1_2_0003)))},

            //3 ends 2 tiles 3rd end is the same as the first
            new Object[] {1, "38T38B76T76B",
                makeList(new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001,s_1_1_0002)), new IlluminaFileMap(makeList(1,2), getQseqs(s_1_2_0001, s_1_2_0002)), new IlluminaFileMap(makeList(1,2), getQseqs(s_1_1_0001,s_1_1_0002))),
               combineReads(228, new int[]{38, 38, 76, 76}, 0, getTiledReadData(makeList(s_1_1_0001, s_1_1_0002)), getTiledReadData(makeList(s_1_2_0001, s_1_2_0002)), getTiledReadData(makeList(s_1_1_0001, s_1_1_0002)))},

            //4 ends 3 tiles
            new Object[] {1, "38T38T76B76B76B",
                makeList(new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_1_0001,s_1_1_0002, s_1_1_0003)), new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_2_0001, s_1_2_0002, s_1_2_0003)), new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_1_0001,s_1_1_0002, s_1_1_0003)),  new IlluminaFileMap(makeList(1,2,3), getQseqs(s_1_2_0001, s_1_2_0002, s_1_2_0003))),
               combineReads(304, new int[]{38, 38, 76, 76, 76}, 0, getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_1_0003)), getTiledReadData(makeList(s_1_2_0001, s_1_2_0002, s_1_2_0003)), getTiledReadData(makeList(s_1_1_0001, s_1_1_0002, s_1_1_0003)), getTiledReadData(makeList(s_1_2_0001, s_1_2_0002, s_1_2_0003)))},
        };
    }

    @Test(dataProvider="qSeqIntegrations")
    public void qSeqParserTest(final int lane, final String readStructure, final List<IlluminaFileMap> readNumberToTile, final Map<Integer, QseqReadData> testAgainst) {
        final OutputMapping om      = new OutputMapping(new ReadStructure(readStructure));
        final QseqParser qseqParser = new QseqParser(lane, readNumberToTile, new OutputMapping(new ReadStructure(readStructure)));
        final Map<Integer, QseqReadData> filteredTestAgainst = filterAllSkips(testAgainst, om);

        int currentRead = 0;
        while(qseqParser.hasNext()) {
            final QseqReadData qseqRead = qseqParser.next();
            final QseqReadData testRead = filteredTestAgainst.get(currentRead);
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

    @DataProvider(name = "validSplitRangeData")
    public static Object [][] validSplitRangeData() {
        return new Object [][] {
          //how output should be structure, lengths of qseqs, output ranges as a result
          {"76T", new int[]{76},     makeList(new Range(0, 75))},
          {"76T", new int[]{38, 38}, makeList(new Range(0, 37), new Range(38, 75))},
          {"76T76S", new int[]{76,76},              makeList(new Range(0, 75))},
          {"76T76S", new int[]{8, 30, 38, 38, 38},  makeList(new Range(0, 7), new Range(8, 37), new Range(38, 75))},
          {"8S76T8S", new int[]{8,76,8},           makeList(new Range(8,83))},
          {"8S76T8S", new int[]{76,8,8},           makeList(new Range(8,75), new Range(76, 83))},
          {"8S76T8S", new int[]{30, 8, 8, 30, 16}, makeList(new Range(8,29), new Range(30,37), new Range(38, 45), new Range(46, 75), new Range(76, 83))},

          {"101T101T",    new int[]{11, 90, 101},     makeList(new Range(0,10),  new Range(11, 100), new Range(101, 201))},
          {"101T16S101T", new int[]{218},             makeList(new Range(0,100), new Range(117, 217))},
          {"101T16S101S", new int[]{51, 25, 25, 117}, makeList(new Range(0,50),  new Range(51, 75), new Range(76, 100))},
          {"101T16S101T", new int[]{51, 25, 25, 117}, makeList(new Range(0,50),  new Range(51, 75), new Range(76, 100), new Range(117, 217))},
          {"101S16T",     new int[]{50, 51, 16},      makeList(new Range(101, 116))},

          {"25T25T25T",         new int[]{60, 15}, makeList(new Range(0,24), new Range(25,49), new Range(50, 59), new Range(60, 74))},
          {"25T30S25T25T",      new int[]{30,75},  makeList(new Range(0,24), new Range(55,79), new Range(80, 104))},
          {"25T30S25T10S25T",   new int[]{125},    makeList(new Range(0,24), new Range(55,79), new Range(90, 114))}
        };
    }

    @Test(dataProvider = "validSplitRangeData")
    public void testSplitOutputRangesOnQseqBoundaries(final String readStructure, final int [] readLengths, final List<Range> expectedRanges) {
        final List<Range> outputRanges = QseqParser.splitOutputRangesOnQseqBoundaries(readLengths, new OutputMapping(new ReadStructure(readStructure)));
        Assert.assertEquals(outputRanges, expectedRanges);
    }

    //TODO: NEED TO DO SOME INVALID ONES

    @DataProvider(name = "outputRangesToInputRangesData")
    public static Object [][] outputRangesToInputRangesData() {
        return new Object [][] {
            //QSeqRead Lengths, OutputRanges
            //  InputRanges per Read
            {
                new int[]{75, 8, 75},      makeList(new Range(0, 74), new Range(75, 82), new Range(83, 157)),
                makeList(
                    makeList(new Range(0, 74)),
                    makeList(new Range(0, 7)),
                    makeList(new Range(0, 74))
                )
            },
            {
                new int[]{75, 8, 75},      makeList(new Range(0, 74), new Range(75, 82), new Range(91, 157)),
                makeList(
                    makeList(new Range(0,  74)),
                    makeList(new Range(0,  7)),
                    makeList(new Range(8,  74))
                )
            },

            {
                new int[]{101, 8, 8, 101}, makeList(new Range(15, 20),  new Range(90, 100), new Range(109,113), new Range(117, 130), new Range(150, 180)),
                makeList(
                    makeList(new Range(15, 20), new Range(90, 100)),
                    new ArrayList<Range>(),
                    makeList(new Range(0, 4)),
                    makeList(new Range(0, 13), new Range(33,63))
                )
            },
            {
                new int[]{250, 250},       makeList(new Range(25, 249), new Range(250, 276), new Range(480, 499)),
                makeList(
                    makeList(new Range(25, 249)),
                    makeList(new Range(0, 26), new Range(230, 249))
                )
            },
            {
                new int[]{310},            makeList(new Range(0, 99), new Range(300, 309)),
                makeList(
                    makeList(new Range(0,99), new Range(300,309))
                )
            }
        };
    }

    @Test(dataProvider = "outputRangesToInputRangesData")
    public void testOutputRangesToInputRanges(final int [] readLengths, final List<Range> outputRanges, final List<List<Range>> expectedResults) {
        final List<List<Range>> actualResults = QseqParser.outputRangesToInputRanges(readLengths, outputRanges);
        Assert.assertEquals(actualResults, expectedResults);
    }

    @DataProvider(name = "outputRangesTo2DTargetsPerReadData")
    public static Object [][] outputRangesTo2DTargetsPerReadData() {
        return new Object [][] {
            //QSeqRead Lengths, OutputRanges
            //  InputRanges per Read
            {
                "76T8B76T",
                makeList(1,1,1),
                makeList(new Range(0, 75), new Range(76, 83), new Range(84, 157)),
                makeList(
                    makeList(new TwoDIndex(0, 0)),
                    makeList(new TwoDIndex(1, 0)),
                    makeList(new TwoDIndex(2, 0))
                )
            },
            {
                "76T8B8S68T",
                makeList(1,1,1),
                makeList(new Range(0, 75), new Range(76, 83), new Range(92, 158)),
                makeList(
                    makeList(new TwoDIndex(0,  0)),
                    makeList(new TwoDIndex(1,  0)),
                    makeList(new TwoDIndex(2,  0))
                )
            },

            { //less realistic scenario, but useful for testing
                "101T8B8B101T",
                makeList(2,0,1, 2),
                makeList(new Range(15, 20),  new Range(90, 100), new Range(109,113), new Range(117, 130), new Range(150, 180)),
                makeList(
                    makeList(new TwoDIndex(0, 15), new TwoDIndex(0, 90)),
                    new ArrayList<TwoDIndex>(),
                    makeList(new TwoDIndex(2, 0)),
                    makeList(new TwoDIndex(3, 0), new TwoDIndex(3,33))
                )
            },
            {  //Another less realistic scenario
                "250T27T190S43T",
                makeList(1,1,1),
                makeList(new Range(25, 249), new Range(250, 276), new Range(480, 499)),
                makeList(
                    makeList(new TwoDIndex(0, 25)),
                    makeList(new TwoDIndex(1, 0)),
                    makeList(new TwoDIndex(2, 13))
                )
            },
            {
                "100T200S10T8S",
                makeList(1,1),
                makeList(new Range(0, 99), new Range(300, 309)),
                makeList(
                    makeList(new TwoDIndex(0,0)),
                    makeList(new TwoDIndex(1,0))
                )
            }
        };
    }

    @Test(dataProvider = "outputRangesTo2DTargetsPerReadData")
    public void testOutputRangesTo2DTargetsPerRead(String readStructure, final List<Integer> rangesPerRead, final List<Range> splitOutputRanges, final List<List<TwoDIndex>> expectedResults) {
        final List<List<TwoDIndex>> actualResults = QseqParser.outputRangesTo2DTargetsPerRead(rangesPerRead, splitOutputRanges, new OutputMapping(new ReadStructure(readStructure)));
        Assert.assertEquals(actualResults, expectedResults);
    }
}

