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
    private static final int DEFAULT_ENDS            = 4;
    private static final int DEFAULT_CYCLES          = 20;

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

        Assert.assertEquals(fileUtil.qseq().numberOfEnds(), DEFAULT_ENDS);
        Assert.assertEquals(fileUtil.bcl().getNumCycles(),  DEFAULT_CYCLES);

        Assert.assertEquals(fileUtil.getTiles(Arrays.asList(SupportedIlluminaFormat.values())), DEFAULT_TILES);
    }

    @Test
    public void passNewUtilTest() {
        for(final SupportedIlluminaFormat format : SupportedIlluminaFormat.values()) {
            makeFiles(format, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_ENDS);
            makeFiles(format, DEFAULT_LANE+1, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_ENDS, ".gz");
            makeFiles(format, DEFAULT_LANE+2, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_ENDS, ".bz2");
        }

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), DEFAULT_LANE);
        assertDefaults(fileUtil, null);

        final IlluminaFileUtil fileUtil2 = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), DEFAULT_LANE+1);
        Assert.assertEquals(fileUtil.getTiles(Arrays.asList(SupportedIlluminaFormat.values())), DEFAULT_TILES);
        assertDefaults(fileUtil2, DEFAULT_LANE+1);

        final IlluminaFileUtil fileUtil3 = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), DEFAULT_LANE+2);
        Assert.assertEquals(fileUtil.getTiles(Arrays.asList(SupportedIlluminaFormat.values())), DEFAULT_TILES);
        assertDefaults(fileUtil3, DEFAULT_LANE+2);
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

    @Test(dataProvider="missingTileFormats")
    public void missingTileTest(final int lane,
                                final List<SupportedIlluminaFormat> formats,
                                final List<SupportedIlluminaFormat> formatsToGetTiles,
                                final List<String> relativeFilesToDelete,
                                final String compression) {
        for(final SupportedIlluminaFormat format : formats) {
            makeFiles(format, lane, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_ENDS, compression);
        }

        for(final String relativeFile : relativeFilesToDelete) {
            final File actualFile = new File(intensityDir, relativeFile);
            if(!actualFile.exists()) {
                throw new PicardException("Trying to delete a non-existent file" + actualFile.getAbsolutePath());
            }
            actualFile.delete();
            if(actualFile.exists()) {
                throw new PicardException("File still exists after calling delete: " + actualFile);
            }
        }

        PicardException pExc = null;
        try{
            final IlluminaFileUtil fUtil = new IlluminaFileUtil(new File(intensityDir, "BaseCalls"), lane);
            fUtil.getTiles(formatsToGetTiles);
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
        makeFiles(format, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_ENDS, compression);

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

    public void testDefaultPerTilePerCycleUtil(final IlluminaFileUtil.PerTilePerCycleFileUtil pcfu, final File parentDir) {
        final CycleIlluminaFileMap cfm       = pcfu.getFiles();
        final CycleIlluminaFileMap cfmWTiles = pcfu.getFiles(DEFAULT_TILES);

        Assert.assertEquals(cfm.size(), DEFAULT_TILES.size());

        for(final Integer tile : DEFAULT_TILES) {
            final CycleFilesIterator tFileIter  = cfm.get(tile);
            final CycleFilesIterator tFileIter2 = cfmWTiles.get(tile);

            for(int cycle = 1; cycle <= DEFAULT_CYCLES; cycle++) {
                final File tcFile = tFileIter.next();
                final File tcFile2 = tFileIter2.next();

                Assert.assertEquals(tcFile.getAbsolutePath(), tcFile2.getAbsolutePath());
                Assert.assertEquals(tcFile, makePerTilePerCycleFilePath(parentDir, DEFAULT_LANE, tile, cycle, pcfu.extension));
                Assert.assertTrue(tcFile.exists());
                Assert.assertTrue(tcFile.length() > 0);
            }
        }

        final List<Integer> tiles = new ArrayList<Integer>(DEFAULT_TILE_TEST_SUBSET);
        final CycleIlluminaFileMap subsetMap = pcfu.getFiles(DEFAULT_TILE_TEST_SUBSET);
        for(final Integer tile : subsetMap.keySet()) {
            tiles.remove(tile);
            Assert.assertTrue(DEFAULT_TILE_TEST_SUBSET.contains(tile));
            final CycleFilesIterator tFileIter  = subsetMap.get(tile);

            for(int cycle = 1; cycle <= DEFAULT_CYCLES; cycle++) {
                final File tcFile = tFileIter.next();
                Assert.assertEquals(tcFile, makePerTilePerCycleFilePath(parentDir, DEFAULT_LANE, tile, cycle, pcfu.extension));
                Assert.assertTrue(tcFile.exists());
                Assert.assertTrue(tcFile.length() > 0);
            }
        }

        Assert.assertTrue(tiles.isEmpty());
    }

    @DataProvider(name="perTilePerCycleFileFormats")
    public Object[][] perTilePerCycleFileFormats() {
        return new Object[][]{
            {SupportedIlluminaFormat.Bcl,    "BaseCalls/" + laneDir(DEFAULT_LANE) },
            {SupportedIlluminaFormat.Cif,    laneDir(DEFAULT_LANE)                },
            {SupportedIlluminaFormat.Cnf,    laneDir(DEFAULT_LANE)                }
        };
    }

    @Test(dataProvider="perTilePerCycleFileFormats")
    public void perTilePerCycleFileUtilsTest(final SupportedIlluminaFormat format, final String parentDir) {
        makeFiles(format, DEFAULT_LANE, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_ENDS, null);

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, DEFAULT_LANE);
        final IlluminaFileUtil.PerTilePerCycleFileUtil pcfu = (IlluminaFileUtil.PerTilePerCycleFileUtil)fileUtil.getUtil(format);

        Assert.assertTrue(pcfu.filesAvailable());
        testDefaultPerTilePerCycleUtil(pcfu, (parentDir == null) ? intensityDir : new File(intensityDir, parentDir));

        final IlluminaFileUtil noFilesFu = new IlluminaFileUtil(basecallDir, DEFAULT_LANE+20);
        final IlluminaFileUtil.PerTilePerCycleFileUtil noFilesPcfu = (IlluminaFileUtil.PerTilePerCycleFileUtil)noFilesFu.getUtil(format);

        Assert.assertFalse(noFilesPcfu.filesAvailable());
        Assert.assertTrue(noFilesPcfu.getFiles().isEmpty());
        Assert.assertTrue(noFilesPcfu.getFiles(DEFAULT_TILES).isEmpty());
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
        makeFiles(SupportedIlluminaFormat.Qseq, lane, DEFAULT_TILES, DEFAULT_CYCLES, DEFAULT_ENDS, compression);

        final IlluminaFileUtil fileUtil = new IlluminaFileUtil(basecallDir, lane);
        final IlluminaFileUtil.QSeqIlluminaFileUtil qseq = fileUtil.qseq();

        Assert.assertTrue(qseq.filesAvailable());
        Assert.assertEquals(qseq.numberOfEnds(), DEFAULT_ENDS);
        testQseqUtil(qseq, lane, compression);

        final IlluminaFileUtil noFilesFu = new IlluminaFileUtil(basecallDir, DEFAULT_LANE+20);
        Assert.assertFalse(noFilesFu.qseq().filesAvailable());
        Assert.assertEquals(qseq.numberOfEnds(), DEFAULT_ENDS);
        Assert.assertTrue(noFilesFu.qseq().getFiles().isEmpty());
        Assert.assertTrue(noFilesFu.qseq().getFiles(DEFAULT_TILES).isEmpty());
    }

    public void testQseqUtil(final IlluminaFileUtil.QSeqIlluminaFileUtil qseq, final int lane, final String compression) {
        final List<IlluminaFileMap> listOfFm       = qseq.getFiles();
        final List<IlluminaFileMap> listOfFmWTiles = qseq.getFiles(DEFAULT_TILES);

        Assert.assertEquals(listOfFm.size(), DEFAULT_ENDS);

        for(int i = 0; i < DEFAULT_ENDS; i++) {
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

        for(int i = 0; i < DEFAULT_ENDS; i++) {
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

    private void makeFiles(final SupportedIlluminaFormat format, int lane, List<Integer> tiles, final int cycles, final int ends) {
        makeFiles(format, lane, tiles, cycles, ends, null);
    }
    private void makeFiles(final SupportedIlluminaFormat format, int lane, List<Integer> tiles, final int cycles, final int ends, final String compression) {
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
                makePerTilePerCycleFiles(basecallLaneDir,      lane, tiles, cycles, ".bcl");
                break;

            case Cif:
                makePerTilePerCycleFiles(intensityLaneDir, lane, tiles, cycles, ".cif");
                break;

            case Cnf:
                makePerTilePerCycleFiles(intensityLaneDir, lane, tiles, cycles, ".cnf");
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

    private static void makePerTilePerCycleFiles(final File parentDir, final int lane, final List<Integer> tiles, final int numCycles, final String ext) {
        if(!parentDir.exists()) {
            if(!parentDir.mkdir()) {
                throw new RuntimeException("Couldn't create directory " + parentDir.getAbsolutePath());
            }
        }

        for(int cycle = 1; cycle <= numCycles; cycle++) {
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

    private static List<String> makeCycleFileList(final File dir, final String ext, final int lane, final int numCycles, final int ... tiles) {
        return makeCycleFileList(dir, ext, lane, numCycles, false, tiles);
    }

    private static List<String> makeCycleFileList(final File dir, final String ext, final int lane, final int numCycles, final boolean longFmt, final int ... tiles) {
        final List<String> files = new ArrayList<String>();
        final File laneDir = new File(dir, laneDir(lane));

        for(int cycle = 1; cycle <= numCycles; cycle++) {
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
