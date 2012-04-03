package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import static net.sf.picard.util.CollectionUtil.makeList;
import org.testng.Assert;
import org.testng.annotations.*;

import net.sf.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class IlluminaFileUtilTest {
    private static final int DEFAULT_LANE            = 7;
    private static final List<Integer> DEFAULT_TILES            = makeList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
    private static final List<Integer> DEFAULT_TILE_TEST_SUBSET = makeList(1, 4, 5, 6, 9, 10);
    private static final int DEFAULT_NUM_ENDS    = 4;
    private static final int [] DEFAULT_ENDS     = {1,2,3,4};
    private static final int [] DEFAULT_CYCLES   = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    private static final int DEFAULT_LAST_CYCLE  = 20;

    private File intensityDir;
    private File basecallDir;

    @BeforeMethod
    private void setUp() throws Exception {
        intensityDir = IoUtil.createTempDir("ift_test", "Intensities");
        basecallDir = new File(intensityDir, "BaseCalls");
        if(!basecallDir.mkdir()) {
            throw new RuntimeException("Couldn't make basecalls dir " + basecallDir.getAbsolutePath());
        }
    }

    @AfterMethod
    private void tearDown() {
        IoUtil.deleteDirectoryTree(intensityDir);
    }

    @DataProvider(name="validLanes")
    public Object[][] validLanes() {
        return new Object[][] {
            {0,  "s_0_1111", "s_0_1_1111_qseq.txt"},
            {1,  "s_1_23",   "s_1_9_0023_qseq.txt.gz"},
            {10, "s_10_1",   "s_10_7_0001_qseq.txt.bz2"}
        };
    }
    public void regexMatches(final String regex, final String toMatch) {
        regexMatches(regex, toMatch, true);
    }

    public void regexMatches(final String regex, final String toMatch, boolean expectedResult) {
        final Pattern pt = Pattern.compile(regex);
        final Matcher ma = pt.matcher(toMatch);
        Assert.assertEquals(ma.matches(), expectedResult);
    }

    @Test(dataProvider="validLanes")
    public void regexTests(final int lane, final String ltExample, final String qseqExample) {
        regexMatches(IlluminaFileUtil.makeParameterizedLaneAndTileRegex(lane), ltExample);
        regexMatches(IlluminaFileUtil.makeParameterizedQseqRegex(lane),        qseqExample);
    }

    @DataProvider(name="validLanesInvalidRegexes")
    public Object[][] validLanesInvalidRegexes() {
        return new Object[][] {
                {0,  "s_-0_111", "  s_-1_1_1111_qseq.txt"},
                {1,  "s_1_A3",     "s_1_9_0C23_qseq.txt.gz"},
                {10, "s_-100_1",   "s_10_20_0001_qseq.txt.bz2"},
                {20, "s_21_1",     "s_20_7_0001_qseq.txt.bz3"}
        };
    }
    @Test(dataProvider="validLanesInvalidRegexes")
    public void notMatchingRegexTest(final int lane, final String ltExample, final String qseqExample) {
        regexMatches(IlluminaFileUtil.makeParameterizedLaneAndTileRegex(lane), ltExample,   false);
        regexMatches(IlluminaFileUtil.makeParameterizedQseqRegex(lane),        qseqExample, false);
    }

    @DataProvider(name="invalidLanes")
    public Object[][] invalidLanes() {
        return new Object[][] {
                {-1000},
                {-10},
                {-1}
        };
    }

    @Test(dataProvider="invalidLanes", expectedExceptions = PicardException.class)
    public void invalidLaneForLTRegex(final int lane) {
        IlluminaFileUtil.makeParameterizedLaneAndTileRegex(lane);
    }

    @Test(dataProvider="invalidLanes", expectedExceptions = PicardException.class)
    public void invalidLaneForQseqRegex(final int lane) {
        IlluminaFileUtil.makeParameterizedQseqRegex(lane);
    }

    public void assertDefaults(final IlluminaFileUtil fileUtil, final Integer lane) {
        if(lane == null) {
            Assert.assertEquals(fileUtil.getLane(), DEFAULT_LANE);
        } else {
            Assert.assertEquals(new Integer(fileUtil.getLane()), lane);
        }

        Assert.assertEquals(fileUtil.barcode(), fileUtil.getUtil(SupportedIlluminaFormat.Barcode));
        Assert.assertEquals(fileUtil.barcode().getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.bcl(),     fileUtil.getUtil(SupportedIlluminaFormat.Bcl));
        Assert.assertEquals(fileUtil.bcl().getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.cif(),     fileUtil.getUtil(SupportedIlluminaFormat.Cif));
        Assert.assertEquals(fileUtil.cif().getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.cnf(),     fileUtil.getUtil(SupportedIlluminaFormat.Cnf));
        Assert.assertEquals(fileUtil.cnf().getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.pos(),     fileUtil.getUtil(SupportedIlluminaFormat.Pos));
        Assert.assertEquals(fileUtil.pos().getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.locs(),    fileUtil.getUtil(SupportedIlluminaFormat.Locs));
        Assert.assertEquals(fileUtil.locs().getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.clocs(),   fileUtil.getUtil(SupportedIlluminaFormat.Clocs));
        Assert.assertEquals(fileUtil.clocs().getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.filter(),  fileUtil.getUtil(SupportedIlluminaFormat.Filter));
        Assert.assertEquals(fileUtil.filter().getTiles(), DEFAULT_TILES);

        Assert.assertEquals(fileUtil.qseq(),    fileUtil.getUtil(SupportedIlluminaFormat.Qseq));
        Assert.assertEquals(fileUtil.qseq().getTiles(),   DEFAULT_TILES);

        Assert.assertEquals(fileUtil.qseq().numberOfEnds(),      DEFAULT_NUM_ENDS);

        final int [] detectedCycles = fileUtil.bcl().getDetectedCycles();
        Assert.assertEquals(detectedCycles.length, DEFAULT_CYCLES.length);
        for(int i = 0; i < DEFAULT_CYCLES.length; i++) {
            Assert.assertEquals(detectedCycles[i],  DEFAULT_CYCLES[i], "Elements differ at index " + i);
        }

        Assert.assertEquals(fileUtil.getActualTiles(Arrays.asList(SupportedIlluminaFormat.values())), DEFAULT_TILES);
    }

    @Test
    public void passNewUtilTest() {
        for(final SupportedIlluminaFormat format : SupportedIlluminaFormat.values()) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES,   DEFAULT_NUM_ENDS);
            makeFiles(format, intensityDir, DEFAULT_LANE+1, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_NUM_ENDS, ".gz");
            makeFiles(format, intensityDir, DEFAULT_LANE+2, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_NUM_ENDS, ".bz2");
        }

        for(int i = 0; i < 3; i++) {
            final IlluminaFileUtil fileUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), DEFAULT_LANE + i);
            Assert.assertEquals(fileUtil.getActualTiles(Arrays.asList(SupportedIlluminaFormat.values())), DEFAULT_TILES);
            assertDefaults(fileUtil, DEFAULT_LANE+i);
        }
    }

    @Test
    public void passingVerifyTest_noQseq() {
        for(final SupportedIlluminaFormat format : SupportedIlluminaFormat.values()) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES,   DEFAULT_NUM_ENDS);
            makeFiles(format, intensityDir, DEFAULT_LANE+1, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_NUM_ENDS, ".gz");
            makeFiles(format, intensityDir, DEFAULT_LANE+2, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_NUM_ENDS, ".bz2");
        }

        for(int i = 0; i < 3; i++) {
            final IlluminaFileUtil fileUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), DEFAULT_LANE + i);


            for(final SupportedIlluminaFormat format : SupportedIlluminaFormat.values()) {
                if(format != SupportedIlluminaFormat.Qseq) {
                    Assert.assertEquals(new ArrayList<String>(), fileUtil.getUtil(format).verify(DEFAULT_TILES, DEFAULT_CYCLES));
                }
            }
        }
    }

    @DataProvider(name="passingQseqVerifyTestData")
    public Object[][] passingQseqVerifyTestData() {
        return new Object[][] {
            {4, makeList(1,2,3), cycleRange(1,   152)},
            {6, makeList(1),     cycleRange(1,   144)},
            {6, makeList(1),     cycleRange(100, 110)},
        };
    }

    @Test(dataProvider="passingQseqVerifyTestData")
    public void passingVerifyTest_qseq(final int lane, List<Integer> tiles, final int [] cycles) {
        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(new File("testdata/net/sf/picard/illumina/IlluminaTests/", "BasecallsDir"), lane);
        Assert.assertEquals(new ArrayList<String>(), fileUtil.qseq().verify(tiles, cycles));
    }

    //Need one where the cycle range asked for is greater than that available
    //Need one where the various reads are missing

    @DataProvider(name="failingQseqVerifyTestData")
    public Object[][] failingQseqVerifyTestData() {
        return new Object[][] {
             {4, makeList(1,2,3), cycleRange(1,   153), new ArrayList<String>(),                   1},
             {4, makeList(1,2,3), cycleRange(1,   152), makeList("BaseCalls/s_4_2_0002_qseq.txt"), 1},
             {4, makeList(1,2,3), cycleRange(1,   10),  makeList("BaseCalls/s_4_1_0002_qseq.txt"), 1},
             {4, makeList(1,2,3), cycleRange(1,   10),  makeList("BaseCalls/s_4_2_0002_qseq.txt"), 1},
             {6, makeList(1),     cycleRange(1,   146), makeList("BaseCalls/s_6_2_0001_qseq.txt"), 78},
             {6, makeList(1),     cycleRange(100, 110), makeList("BaseCalls/s_6_3_0001_qseq.txt"), 11},
        };
    }
    @Test(dataProvider="failingQseqVerifyTestData")
    public void failingVerifyTest_qseq(final int lane, List<Integer> tiles, final int [] cycles, final List<String> toDelete, int numErrors) {
        IoUtil.copyDirectoryTree(new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir"), basecallDir);

        deleteRelativeFiles(toDelete);

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, lane);
        Assert.assertEquals(numErrors, fileUtil.qseq().verify(tiles, cycles).size());
    }

    @DataProvider(name="failingVerifyNoQSeq")
    public Object[][] failingVerify_noQseq() {
        return new Object[][] {
            { makeList("BaseCalls/L007/C4.1/s_7_9.bcl"),                                 1},
            { makeList("BaseCalls/L007/C20.1/s_7_11.bcl", "BaseCalls/L007/C20.1/s_7_12.bcl"), 2},
            { makeList("BaseCalls/L007/C10.1/"), 1},
            { makeList("L007/C1.1/s_7_1.cif", "L007/C20.1/s_7_12.cnf"), 2},

            { makeList("L007/s_7_2.locs", "L007/s_7_3.locs", "L007/s_7_4.locs"),                                 3},
            { makeList("BaseCalls/L007/C20.1/s_7_11.bcl", "BaseCalls/L007/C20.1/s_7_12.bcl", "L007/s_7_4.locs"), 3},
            { makeList("BaseCalls/L007/s_7_0005.filter", "BaseCalls/L007/s_7_0011.filter", "BaseCalls/s_7_0002_barcode.txt"), 3},
        };
    }

    @Test(dataProvider="failingVerifyNoQSeq")
    public void failingVerifyTest_noQseq(final List<String> filesToDelete, final int expectedNumDifferences) {
        for(final SupportedIlluminaFormat format : SupportedIlluminaFormat.values()) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES,   DEFAULT_NUM_ENDS);
        }

        deleteRelativeFiles(filesToDelete);

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), DEFAULT_LANE);

        int totalDifferences = 0;
        List<String> differences = new ArrayList<String>();
        for(final SupportedIlluminaFormat format : SupportedIlluminaFormat.values()) {
            if(format != SupportedIlluminaFormat.Qseq) {
                final List<String> curDiffs = fileUtil.getUtil(format).verify(DEFAULT_TILES, DEFAULT_CYCLES);
                differences.addAll(curDiffs);

                totalDifferences += curDiffs.size();
            }
        }

        Assert.assertEquals(expectedNumDifferences, totalDifferences);

    }

    @DataProvider(name="missingTileFormats")
    public Object[][] missingTileFormats() {
        return new Object[][]{
            {
                1,
                makeList(SupportedIlluminaFormat.Bcl, SupportedIlluminaFormat.Barcode),
                makeList(SupportedIlluminaFormat.Bcl, SupportedIlluminaFormat.Barcode),
                makeList("BaseCalls/s_1_0007_barcode.txt.gz"),
                ".gz"
            },

            {
                2,
                Arrays.asList(SupportedIlluminaFormat.values()),
                Arrays.asList(SupportedIlluminaFormat.values()),
                makeCycleFileList(new File("BaseCalls"), ".bcl", 2, DEFAULT_CYCLES, 2),
                ".gz"
            },
            {
                3,
                Arrays.asList(SupportedIlluminaFormat.values()),
                Arrays.asList(SupportedIlluminaFormat.values()),
                makeList("BaseCalls/L003/C1.1/s_3_2.bcl"),
                ".bz2"
            },
            {
                4,
                Arrays.asList(SupportedIlluminaFormat.values()),
                Arrays.asList(SupportedIlluminaFormat.Pos, SupportedIlluminaFormat.Locs),
                makeList("s_4_10_pos.txt", "L004/s_4_2.locs"),
                null
            },
            {
                5,
                makeList(SupportedIlluminaFormat.Cif, SupportedIlluminaFormat.Cnf, SupportedIlluminaFormat.Qseq),
                makeList(SupportedIlluminaFormat.Cif, SupportedIlluminaFormat.Qseq),
                makeList("BaseCalls/s_5_1_0003_qseq.txt.gz", "BaseCalls/s_5_2_0003_qseq.txt.gz",
                         "BaseCalls/s_5_3_0003_qseq.txt.gz", "BaseCalls/s_5_4_0003_qseq.txt.gz",
                         "BaseCalls/s_5_1_0004_qseq.txt.gz", "BaseCalls/s_5_2_0004_qseq.txt.gz",
                         "BaseCalls/s_5_3_0004_qseq.txt.gz", "BaseCalls/s_5_4_0004_qseq.txt.gz"),
                ".gz"
            }
        };
    }

    public static final void emptyRelativeFiles(final File baseFile, final List<String> relativeFilesToDelete) {
        for(final String relativeFile : relativeFilesToDelete) {
            final File actualFile = new File(baseFile, relativeFile);


            if(!actualFile.exists()) {
                throw new RuntimeException("Trying to empty a non-existent file" + actualFile.getAbsolutePath());
            }

            if(actualFile.isDirectory()) {
                throw new RuntimeException("Trying to empty a directory(" + actualFile.getAbsolutePath() +")");
            } else {
                if(!actualFile.delete()) {
                    throw new RuntimeException("Couldn't remove previous file when emptying(" + actualFile.getAbsolutePath() + ")");
                } else {
                    try {
                        if(!actualFile.createNewFile()) {
                            throw new RuntimeException("Couldn't create empty file: " + actualFile.getAbsolutePath() + ")");
                        }
                    } catch(IOException ioe) {
                        throw new RuntimeException(ioe);
                    }
                }
            }
            if(!actualFile.exists()) {
                throw new PicardException("File should exist: " + actualFile);
            }
        }
    }

    public static final void deleteRelativeFiles(final File baseFile, final List<String> relativeFilesToDelete) {
        for(final String relativeFile : relativeFilesToDelete) {
            final File actualFile = new File(baseFile, relativeFile);


            if(!actualFile.exists()) {
                throw new RuntimeException("Trying to delete a non-existent file" + actualFile.getAbsolutePath());
            }

            if(actualFile.isDirectory()) {
                IoUtil.deleteDirectoryTree(actualFile);
            } else {
                actualFile.delete();
            }
            if(actualFile.exists()) {
                throw new RuntimeException("File still exists after calling delete: " + actualFile);
            }
        }
    }

    public final void deleteRelativeFiles(final List<String> relativeFilesToDelete) {
        deleteRelativeFiles(intensityDir, relativeFilesToDelete);
    }

    @Test(dataProvider="missingTileFormats")
    public void missingTileTest(final int lane,
                                final List<SupportedIlluminaFormat> formats,
                                final List<SupportedIlluminaFormat> formatsToGetTiles,
                                final List<String> relativeFilesToDelete,
                                final String compression) {
        for(final SupportedIlluminaFormat format : formats) {
            makeFiles(format, intensityDir, lane, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_NUM_ENDS, compression);
        }

        deleteRelativeFiles(relativeFilesToDelete);

        PicardException pExc = null;
        try{
            final IlluminaFileUtil fUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), lane);
            fUtil.getActualTiles(formatsToGetTiles);
        } catch(final PicardException exception) {
            pExc = exception;
        }

        Assert.assertNotNull(pExc, "Didn't raise a Picard Exception for missing tile!");
        Assert.assertTrue(pExc.getMessage().contains("Formats do not have the same number of tiles! "), "Wrong exception thrown for missing tile!");
    }

    @DataProvider(name="perTileFileFormats")
    public Object[][] perTileFileUtils() {
        return new Object[][]{
            {SupportedIlluminaFormat.Locs,    null,   false, laneDir(DEFAULT_LANE)},
            {SupportedIlluminaFormat.Clocs,   null,   false, laneDir(DEFAULT_LANE)},
            {SupportedIlluminaFormat.Pos,     ".gz",  false, null},
            {SupportedIlluminaFormat.Pos,     null,   false, null},
            {SupportedIlluminaFormat.Filter,  null,   true,  "BaseCalls/" + laneDir(DEFAULT_LANE)},
            {SupportedIlluminaFormat.Barcode, ".bz2", true,  "BaseCalls"}
        };
    }

    public File makePerTileFile(final File parentDir, final int lane, final int tile, final String extension, final String compression, final boolean longFormat) {
        return new File(parentDir, "s_" + lane + "_" + longTile(tile, longFormat) + extension + (compression != null ? compression : ""));
    }

    public void testDefaultPerTileUtil(final IlluminaFileUtil.PerTileFileUtil ptfu, final String compression, final boolean longFormat, final File parentDir) {
        final IlluminaFileMap fm = ptfu.getFiles();
        final IlluminaFileMap fmWTiles = ptfu.getFiles(DEFAULT_TILES);

        Assert.assertEquals(fm.size(), DEFAULT_TILES.size());

        for(final Integer tile : DEFAULT_TILES) {
            final File tFile  = fm.get(tile);
            final File tFile2 = fmWTiles.get(tile);
            Assert.assertEquals(tFile.getAbsolutePath(), tFile2.getAbsolutePath());
            Assert.assertEquals(tFile, makePerTileFile(parentDir, DEFAULT_LANE, tile, ptfu.extension, compression, longFormat));
            Assert.assertTrue(tFile.exists());
            Assert.assertTrue(tFile.length() > 0);
        }

        final List<Integer> tiles = new ArrayList<Integer>(DEFAULT_TILE_TEST_SUBSET);
        final IlluminaFileMap subsetMap = ptfu.getFiles(DEFAULT_TILE_TEST_SUBSET);
        for(final Integer tile : subsetMap.keySet()) {
            tiles.remove(tile);
            Assert.assertTrue(DEFAULT_TILE_TEST_SUBSET.contains(tile));
            final File tFile  = subsetMap.get(tile);
            Assert.assertEquals(tFile, makePerTileFile(parentDir, DEFAULT_LANE, tile, ptfu.extension, compression, longFormat));
            Assert.assertTrue(tFile.exists());
            Assert.assertTrue(tFile.length() > 0);
        }

        Assert.assertTrue(tiles.isEmpty());
    }

    @Test(dataProvider="perTileFileFormats")
    public void perTileFileUtilsTest(final SupportedIlluminaFormat format, final String compression, final boolean longFormat, final String parentDir) {
        makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_NUM_ENDS, compression);

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, DEFAULT_LANE);
        final IlluminaFileUtil.PerTileFileUtil ptfu = (IlluminaFileUtil.PerTileFileUtil)fileUtil.getUtil(format);

        Assert.assertTrue(ptfu.filesAvailable());
        testDefaultPerTileUtil(ptfu, compression, longFormat, (parentDir == null) ? intensityDir : new File(intensityDir, parentDir));

        final IlluminaFileUtil noFilesFu = new IlluminaFileUtil(basecallDir, DEFAULT_LANE+20);
        final IlluminaFileUtil.PerTileFileUtil noFilesPtfu = (IlluminaFileUtil.PerTileFileUtil)noFilesFu.getUtil(format);
        Assert.assertFalse(noFilesPtfu.filesAvailable());
        Assert.assertTrue(noFilesPtfu.getFiles().isEmpty());
        Assert.assertTrue(noFilesPtfu.getFiles(DEFAULT_TILES).isEmpty());
    }

    public File makePerTilePerCycleFilePath(final File parentDir, final int lane, final int tile, final int cycle, final String extension) {
        return new File(parentDir, "C" + cycle + ".1/s_" + lane + "_" + tile + extension);
    }

    public void testDefaultPerTilePerCycleUtil(final IlluminaFileUtil.PerTilePerCycleFileUtil pcfu, final File parentDir, final int [] cycles) {
        final CycleIlluminaFileMap cfm         = pcfu.getFiles(cycles);
        final CycleIlluminaFileMap cfmWTiles   = pcfu.getFiles(DEFAULT_TILES, cycles);
        final CycleIlluminaFileMap cfmNoCycles;
        if(Arrays.equals(cycles, DEFAULT_CYCLES)) {
            cfmNoCycles = pcfu.getFiles();
        } else {
            cfmNoCycles = null;
        }

        Assert.assertEquals(cfm.size(), DEFAULT_TILES.size());

        for(final Integer tile : DEFAULT_TILES) {
            final CycleFilesIterator tFileIter  = cfm.get(tile);
            final CycleFilesIterator tFileIter2 = cfmWTiles.get(tile);
            final CycleFilesIterator tFileIter3;
            if(cfmNoCycles != null) {
                tFileIter3 = cfmNoCycles.get(tile);
            } else {
                tFileIter3 = null;
            }

            for(final int cycle : cycles) {
                final File tcFile = tFileIter.next();
                final File tcFile2 = tFileIter2.next();

                Assert.assertEquals(tcFile.getAbsolutePath(), tcFile2.getAbsolutePath());
                if(tFileIter3 != null) {
                    final File tfFile3 = tFileIter3.next();
                    Assert.assertEquals(tcFile.getAbsolutePath(), tfFile3.getAbsolutePath());
                }

                Assert.assertEquals(tcFile, makePerTilePerCycleFilePath(parentDir, DEFAULT_LANE, tile, cycle, pcfu.extension));
                Assert.assertTrue(tcFile.exists());
                Assert.assertTrue(tcFile.length() > 0);
            }
        }
    }


    public void testSubsetDefaultPerTilePerCycleUtil(final IlluminaFileUtil.PerTilePerCycleFileUtil pcfu, final File parentDir, int [] cycles) {
        final List<Integer> tiles = new ArrayList<Integer>(DEFAULT_TILE_TEST_SUBSET);
        final CycleIlluminaFileMap subsetMap = pcfu.getFiles(DEFAULT_TILE_TEST_SUBSET, cycles);
        final CycleIlluminaFileMap cfmNoCycles;
        if(Arrays.equals(cycles, DEFAULT_CYCLES)) {
            cfmNoCycles = pcfu.getFiles(DEFAULT_TILE_TEST_SUBSET);
        } else {
            cfmNoCycles = null;
        }

        for(final Integer tile : subsetMap.keySet()) {
            tiles.remove(tile);
            Assert.assertTrue(DEFAULT_TILE_TEST_SUBSET.contains(tile));
            final CycleFilesIterator tFileIter  = subsetMap.get(tile);
            final CycleFilesIterator tFileIter2;
            if(cfmNoCycles != null) {
                tFileIter2 = cfmNoCycles.get(tile);
            } else {
                tFileIter2 = null;
            }


            for(final int cycle : cycles) {
                final File tcFile = tFileIter.next();
                if(tFileIter2 != null) {
                    Assert.assertEquals(tcFile, tFileIter2.next());
                }
                Assert.assertEquals(tcFile, makePerTilePerCycleFilePath(parentDir, DEFAULT_LANE, tile, cycle, pcfu.extension));
                Assert.assertTrue(tcFile.exists());
                Assert.assertTrue(tcFile.length() > 0);
            }

            Assert.assertFalse(tFileIter.hasNext());
        }

        Assert.assertTrue(tiles.isEmpty());
    }

    public static int [] cycleRange(final Range range) {
        return cycleRange(range.start, range.end);
    }

    public static int [] cycleRange(final int start, final int end) {
        final int [] cycles = new int[end - start + 1];
        for(int i = 0; i < cycles.length; i++) {
            cycles[i] = start + i;
        }

        return cycles;
    }

    public static int [] cycleRange(final int end) {
        return cycleRange(1, end);
    }

    @DataProvider(name="perTilePerCycleFileFormats")
    public Object[][] perTilePerCycleFileFormats() {
        return new Object[][]{
            {SupportedIlluminaFormat.Bcl,    "BaseCalls/" + laneDir(DEFAULT_LANE),  DEFAULT_CYCLES,     false, false },
            {SupportedIlluminaFormat.Bcl,    "BaseCalls/" + laneDir(DEFAULT_LANE),  cycleRange(4),      true,  true  },
            {SupportedIlluminaFormat.Cif,    laneDir(DEFAULT_LANE)               ,  DEFAULT_CYCLES,     false, false },
            {SupportedIlluminaFormat.Cif,    laneDir(DEFAULT_LANE)               ,  cycleRange(8,12),   true,  false },
            {SupportedIlluminaFormat.Cif,    laneDir(DEFAULT_LANE)               ,  cycleRange(8,12),   false, false },
            {SupportedIlluminaFormat.Cnf,    laneDir(DEFAULT_LANE)               ,  DEFAULT_CYCLES,     false, false },
            {SupportedIlluminaFormat.Cnf,    laneDir(DEFAULT_LANE)               ,  cycleRange(15,DEFAULT_LAST_CYCLE), false, false }
        };
    }

    @Test(dataProvider="perTilePerCycleFileFormats")
    public void perTilePerCycleFileUtilsTest(final SupportedIlluminaFormat format, final String parentDir, final int [] cycles, boolean createEarlySkippedCycles, boolean createLateSkippedCycles) {
        if(createEarlySkippedCycles) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, cycleRange(1, cycles[0]), DEFAULT_NUM_ENDS, null);
        }

        makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, cycles, DEFAULT_NUM_ENDS, null);

        if(createLateSkippedCycles) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, cycleRange(cycles[cycles.length-1] + 1, DEFAULT_LAST_CYCLE), DEFAULT_NUM_ENDS, null);
        }

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, DEFAULT_LANE);
        final IlluminaFileUtil.PerTilePerCycleFileUtil pcfu = (IlluminaFileUtil.PerTilePerCycleFileUtil)fileUtil.getUtil(format);

        Assert.assertTrue(pcfu.filesAvailable());
        testDefaultPerTilePerCycleUtil(pcfu, (parentDir == null) ? intensityDir : new File(intensityDir, parentDir), cycles);
        testSubsetDefaultPerTilePerCycleUtil(pcfu, (parentDir == null) ? intensityDir : new File(intensityDir, parentDir), cycles);

        final IlluminaFileUtil noFilesFu = new IlluminaFileUtil(basecallDir, DEFAULT_LANE+20);
        final IlluminaFileUtil.PerTilePerCycleFileUtil noFilesPcfu = (IlluminaFileUtil.PerTilePerCycleFileUtil)noFilesFu.getUtil(format);

        Assert.assertFalse(noFilesPcfu.filesAvailable());
        Assert.assertTrue(noFilesPcfu.getFiles().isEmpty());
        Assert.assertTrue(noFilesPcfu.getFiles(DEFAULT_TILES).isEmpty());
    }

    @DataProvider(name="missingCycleDataRanges")
    public Object[][] missingCycleDataRanges() {
        return new Object[][] {
            {makeList(new Range(10,15))},
            {makeList(new Range(9,12), new Range(14,15))}
        };
    }

    @Test(expectedExceptions = PicardException.class, dataProvider="missingCycleDataRanges")
    public void perTilePerCycleFileUtilsMissingCycleTest(final List<Range> cycleRangesToMake) {
        final SupportedIlluminaFormat format = SupportedIlluminaFormat.Bcl;
        final String parentDir = "BaseCalls/" + laneDir(DEFAULT_LANE);

        for(final Range range : cycleRangesToMake) {
            makeFiles(format, intensityDir, DEFAULT_LANE, DEFAULT_TILES, cycleRange(range), DEFAULT_NUM_ENDS, null);
        }

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, DEFAULT_LANE);
        final IlluminaFileUtil.PerTilePerCycleFileUtil pcfu = (IlluminaFileUtil.PerTilePerCycleFileUtil)fileUtil.getUtil(format);

        Assert.assertTrue(pcfu.filesAvailable());
        int [] cycles = cycleRange(9,16);
        CycleIlluminaFileMap cfm = pcfu.getFiles(cycles);
        cfm.assertValid(DEFAULT_TILES, cycles);
    }

    @DataProvider(name="qseqTestData")
    public Object[][] qseqTestData() {
        return new Object[][]{
            //lane, compression
            {1, ".gz" },
            {8, ".bz2"},
            {10, null }
        };
    }

    @Test(dataProvider="qseqTestData")
    public void qseqFileUtilTest(final int lane, final String compression) {
        makeFiles(SupportedIlluminaFormat.Qseq, intensityDir, lane, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_NUM_ENDS, compression);

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, lane);
        final IlluminaFileUtil.QSeqIlluminaFileUtil qseq = fileUtil.qseq();

        Assert.assertTrue(qseq.filesAvailable());
        Assert.assertEquals(qseq.numberOfEnds(), DEFAULT_NUM_ENDS);
        testQseqUtil(qseq, lane, compression);

        final IlluminaFileUtil noFilesFu = new IlluminaFileUtil(basecallDir, DEFAULT_LANE+20);
        Assert.assertFalse(noFilesFu.qseq().filesAvailable());
        Assert.assertEquals(qseq.numberOfEnds(), DEFAULT_NUM_ENDS);
        Assert.assertTrue(noFilesFu.qseq().getFiles().isEmpty());
        Assert.assertTrue(noFilesFu.qseq().getFiles(DEFAULT_TILES).isEmpty());
    }

    public void testQseqUtil(final IlluminaFileUtil.QSeqIlluminaFileUtil qseq, final int lane, final String compression) {
        final List<IlluminaFileMap> listOfFm       = qseq.getFiles();
        final List<IlluminaFileMap> listOfFmWTiles = qseq.getFiles(DEFAULT_TILES);

        Assert.assertEquals(listOfFm.size(), DEFAULT_NUM_ENDS);

        for(int i = 0; i < DEFAULT_NUM_ENDS; i++) {
            final int currentEnd = i+1;
            final IlluminaFileMap fm = listOfFm.get(i);
            final IlluminaFileMap fmWTiles = listOfFmWTiles.get(i);

            for(final Integer tile : DEFAULT_TILES) {
                final File tFile  = fm.get(tile);
                final File tFile2 = fmWTiles.get(tile);
                Assert.assertEquals(tFile.getAbsolutePath(), tFile2.getAbsolutePath());
                Assert.assertEquals(tFile, makeQSeqFilePath(basecallDir, lane, tile, currentEnd, compression));
                Assert.assertTrue(tFile.exists());
                Assert.assertTrue(tFile.length() > 0);
            }
        }

        final List<IlluminaFileMap> listOfsubsetMap = qseq.getFiles(DEFAULT_TILE_TEST_SUBSET);

        for(int i = 0; i < DEFAULT_NUM_ENDS; i++) {
            final int currentEnd = i+1;
            final IlluminaFileMap subsetMap = listOfsubsetMap.get(i);

            final List<Integer> tiles = new ArrayList<Integer>(DEFAULT_TILE_TEST_SUBSET);

            for(final Integer tile : subsetMap.keySet()) {
                tiles.remove(tile);
                final File tFile  = subsetMap.get(tile);
                Assert.assertEquals(tFile, makeQSeqFilePath(basecallDir, lane, tile, currentEnd, compression));
                Assert.assertTrue(tFile.exists());
                Assert.assertTrue(tFile.length() > 0);
            }
            Assert.assertTrue(tiles.isEmpty());
        }
    }

    private static File makeQSeqFilePath(final File basecallDir, final int lane, final int tile, final int end, final String compression) {
        return new File(basecallDir, "s_" + lane +"_" + end + "_" + longTile(tile,true) + "_qseq.txt" + (compression != null ? compression : ""));
    }

    public static void makeFiles(final SupportedIlluminaFormat format, final File intensityDir, int lane, List<Integer> tiles, final int [] cycles, final int ends) {
        makeFiles(format, intensityDir, lane, tiles, cycles, ends, null);
    }
    public static void makeFiles(final SupportedIlluminaFormat format, final File intensityDir, int lane, List<Integer> tiles, final int [] cycles, final int ends, final String compression) {
        String laneDir = String.valueOf(lane);
        while(laneDir.length() < 3) {
            laneDir = "0" + laneDir;
        }
        laneDir = "L" + laneDir;


        final File basecallDir      = new File(intensityDir, "BaseCalls");
        final File basecallLaneDir  = new File(basecallDir, laneDir);
        final File intensityLaneDir = new File(intensityDir, laneDir);

        switch(format) {
            //per tile formats
            case Barcode:
                makePerTileFiles(basecallDir,      lane, tiles, maybeAddExt("_barcode.txt", compression), true);
                break;

            case Pos:
                makePerTileFiles(intensityDir,     lane, tiles, maybeAddExt("_pos.txt",     compression), false);
                break;

            case Locs:
                makePerTileFiles(intensityLaneDir, lane, tiles, maybeAddExt(".locs", null),   false);
                break;

            case Clocs:
                makePerTileFiles(intensityLaneDir, lane, tiles, maybeAddExt(".clocs", null),  false);
                break;

            case Filter:
                makePerTileFiles(basecallLaneDir,  lane, tiles, maybeAddExt(".filter", null), true);
                break;

            //per tile per cycle formats
            case Bcl:
                makePerTilePerCycleFiles(basecallLaneDir,  lane, tiles, cycles, ".bcl");
                break;

            case Cif:
                makePerTilePerCycleFiles(intensityLaneDir, lane, tiles, cycles,  ".cif");
                break;

            case Cnf:
                makePerTilePerCycleFiles(intensityLaneDir, lane, tiles, cycles,  ".cnf");
                break;

            //Qseq is based on lane, tile, and end
            case Qseq:
                makeQSeqFiles(basecallDir, lane, tiles, ends, compression);
                break;
        }
    }

    private static void makeQSeqFiles(final File basecallDir, final int lane, final List<Integer> tiles, final int ends, final String compressionExt) {
        if(!basecallDir.exists()) {
            if(!basecallDir.mkdir()) {
                throw new RuntimeException("Couldn't create directory " + basecallDir.getAbsolutePath());
            }
        }

        for(final Integer tile : tiles) {
            for(int end = 1; end <= ends; end++) {
                writeNonEmptyFile(new File(basecallDir, "s_" + lane + "_" + end + "_" + longTile(tile, true) + maybeAddExt("_qseq.txt", compressionExt)));
            }
        }
    }

    private static void makePerTileFiles(final File parentDir, final int lane, final List<Integer> tiles, final String ext, final boolean longName) {
        if(!parentDir.exists()) {
            if(!parentDir.mkdir()) {
                throw new RuntimeException("Couldn't create directory " + parentDir.getAbsolutePath());
            }
        }

        for(final Integer tile : tiles) {
                writeNonEmptyFile(new File(parentDir, "s_" + lane + "_" + longTile(tile, longName) + ext));
        }
    }

    private static void makePerTilePerCycleFiles(final File parentDir, final int lane, final List<Integer> tiles, final int [] cycles, final String ext) {
        if(!parentDir.exists()) {
            if(!parentDir.mkdir()) {
                throw new RuntimeException("Couldn't create directory " + parentDir.getAbsolutePath());
            }
        }

        for(final int cycle : cycles) {
            final File cycleDir = new File(parentDir, "C" + cycle + ".1");
            if(!cycleDir.exists()) {
                if(!cycleDir.mkdir()) {
                    throw new RuntimeException("Couldn't create directory " + cycleDir.getAbsolutePath());
                }
            }

            for(final Integer tile : tiles) {
                writeNonEmptyFile(new File(cycleDir, "s_" + lane + "_" + tile + ext));
            }
        }
    }

    private static List<String> makeCycleFileList(final File dir, final String ext, final int lane, final int [] cycles, final int ... tiles) {
        return makeCycleFileList(dir, ext, lane, cycles, false, tiles);
    }

    private static List<String> makeCycleFileList(final File dir, final String ext, final int lane, final int [] cycles, final boolean longFmt, final int ... tiles) {
        final List<String> files = new ArrayList<String>();
        final File laneDir = new File(dir, laneDir(lane));

        for(final int cycle : cycles) {
            final File cycleDir = new File(laneDir, "C" + cycle + ".1");
            for(final Integer tile : tiles) {
                files.add(cycleDir + "/s_" + lane + "_" + longTile(tile, longFmt) + ext);
            }
        }

        return files;
    }

    private static void writeNonEmptyFile(final File file) {
        try {
            final FileWriter fw = new FileWriter(file);
            fw.write("This is just so that the file doesn't appear to be of 0 length");
            fw.flush();
            fw.close();
        } catch (IOException e) {
            throw new RuntimeException("Exception trying to create non-empty file!", e);
        }
    }

    private static String laneDir(final int lane) {
        String ldir = String.valueOf(lane);
        while(ldir.length() < 3) {
            ldir = "0" + ldir;
        }
        return "L" + ldir;
    }

    private static String longTile(final int tile, final boolean makeLong) {
        if(makeLong) {
            String lt = String.valueOf(tile);
            while(lt.length() < 4) {
                lt = "0" + lt;
            }
            return lt;
        } else {
            return String.valueOf(tile);
        }
    }
    
    private static String maybeAddExt(final String fileExt, final String compressionExt) {
        if(compressionExt != null) {
            return fileExt + compressionExt;
        } else {
            return fileExt;
        }
    }
}
