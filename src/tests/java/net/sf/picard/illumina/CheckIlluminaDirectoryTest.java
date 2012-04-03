package net.sf.picard.illumina;

import java.io.*;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.*;

import static net.sf.picard.util.CollectionUtil.makeList;
import static net.sf.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat;
import static net.sf.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat.*;
import static net.sf.picard.illumina.parser.IlluminaDataType.*;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.illumina.parser.IlluminaDataType;
import net.sf.picard.illumina.parser.IlluminaFileUtil;
import net.sf.picard.illumina.parser.IlluminaFileUtilTest;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.CloserUtil;
import org.testng.Assert;
import org.testng.annotations.*;


public class CheckIlluminaDirectoryTest {

    private File illuminaDir;
    private File dataDir;
    private File interopDir;
    private File intensityDir;
    private File basecallDir;

    @BeforeMethod
    private void setUp() throws Exception {
        illuminaDir = IoUtil.createTempDir("ift_test", "IlluminaDir");

        interopDir   = new File(illuminaDir, "InterOp");
        if(!interopDir.mkdir()) {
            throw new RuntimeException("Couldn't make interop dir " + interopDir.getAbsolutePath());
        }

        dataDir = new File(illuminaDir, "Data");
        if(!dataDir.mkdir()) {
            throw new RuntimeException("Couldn't make data dir " + dataDir.getAbsolutePath());
        }

        intensityDir = new File(dataDir, "Intensities");
        if(!intensityDir.mkdir()) {
            throw new RuntimeException("Couldn't make intensity dir " + intensityDir.getAbsolutePath());
        }

        basecallDir  = new File(intensityDir, "BaseCalls");
        if(!basecallDir.mkdir()) {
            throw new RuntimeException("Couldn't make basecalls dir " + basecallDir.getAbsolutePath());
        }
    }

    @AfterMethod
    private void tearDown() {
        IoUtil.deleteDirectoryTree(intensityDir);
    }

    public void makeFiles(SupportedIlluminaFormat [] formats, int lane, List<Integer> tiles, int [] cycles) {
        for(final IlluminaFileUtil.SupportedIlluminaFormat format : formats) {
            IlluminaFileUtilTest.makeFiles(format, intensityDir, lane, tiles, cycles, 0);
        }
    }

    public String [] makeCheckerArgs(final File basecallDir, final int lane, final String readStructure, final IlluminaDataType [] dataTypes) {
        final String [] dataTypeArgs = new String[dataTypes.length + 3];

        dataTypeArgs[0] = "B=" + basecallDir;
        dataTypeArgs[1] = StandardOptionDefinitions.LANE_SHORT_NAME + "=" + lane;
        dataTypeArgs[2] = "RS=" + readStructure;

        for(int i = 0; i < dataTypes.length; i++) {
            dataTypeArgs[i+3] = "DT=" + dataTypes[i];
        }

        return dataTypeArgs;
    }

    public File writeTileMetricsOutFile(Map<Integer, List<Integer>> lanesToTiles) {
        return writeTileMetricsOutFile(interopDir, (byte)2, (byte)10, lanesToTiles);
    }

    public File writeTileMetricsOutFile(final File interopDir, final byte versionNumber, final byte recordSize, Map<Integer, List<Integer>> lanesToTiles) {
        final File tileMetricsOut = new File(interopDir, "TileMetricsOut.bin");
        if(!tileMetricsOut.exists()) {
            try {
                if(!tileMetricsOut.createNewFile()) {
                    throw new PicardException("Could not create tileMetricsOut file(" + tileMetricsOut.getAbsolutePath() + ")");
                }
            } catch (IOException e) {
                throw new PicardException("IOException creating tileMetricsOut file (" + tileMetricsOut + ") for writing!", e);
            }
        }

        int totalEntries = 0;
        for(final Map.Entry<Integer, List<Integer>> l2t : lanesToTiles.entrySet()) {
            totalEntries += l2t.getValue().size();
        }

        final MappedByteBuffer buf;
        try {
            final RandomAccessFile raf = new RandomAccessFile(tileMetricsOut, "rw");
            final FileChannel channel = raf.getChannel();
            buf = channel.map(FileChannel.MapMode.READ_WRITE, 0, 2 + 10 * totalEntries);
            buf.order(ByteOrder.LITTLE_ENDIAN);

            buf.put(versionNumber);
            buf.put(recordSize);

            for(final int lane : lanesToTiles.keySet()) {
                for(final int tile : lanesToTiles.get(lane)) {
                    buf.putShort((short)lane);
                    buf.putShort((short)tile);
                    buf.putShort((short)0);
                    buf.putFloat(0F);
                }
            }

            buf.force();
            CloserUtil.close(channel);
            CloserUtil.close(raf);
        } catch (IOException e) {
            throw new PicardException("IOException writing tileMetricsOut file (" + tileMetricsOut + ")", e);
        }

        return tileMetricsOut;
    }

    public static Map<Integer, List<Integer>> makeMap(final List<Integer> lanes, final List<List<Integer>> tiles) {
        final Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();

        if(lanes.size() != tiles.size()) {
            throw new IllegalArgumentException("Number of lanes (" + lanes + ") does not equal number of tiles!");
        }

        for(int i = 0; i < lanes.size(); i++) {
            map.put(lanes.get(i), tiles.get(i));
        }

        return map;
    }

    @DataProvider(name="positiveTestData")
    public Object [][] positiveTestData() {
        return new Object[][] {
            {
                new SupportedIlluminaFormat[]{Bcl, Locs, Pos, Filter, Qseq}, new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Position, IlluminaDataType.PF},
                3, makeList(1101,1201,1301, 2101,2201,2301), IlluminaFileUtilTest.cycleRange(1,50), "25T25T"
            },
            {
                new SupportedIlluminaFormat[]{Bcl, Locs, Filter}, new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Position, IlluminaDataType.PF},
                2, makeList(1101,1201,1301, 2101,2201,2301), IlluminaFileUtilTest.cycleRange(1,50), "8S15T8S"
            },
            {
                new SupportedIlluminaFormat[]{Bcl, Filter}, new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF},
                2, makeList(1101,1201,1301, 2101,2201,2301), IlluminaFileUtilTest.cycleRange(1,152), "68T8B68T"
            },
            {
                new SupportedIlluminaFormat[]{Bcl, Pos, Filter}, new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Position, IlluminaDataType.PF},
                5, makeList(1101,1201,1301, 2101,2201,2301), IlluminaFileUtilTest.cycleRange(1,50), "25T25T"
            }
        };
    }

    //Note: The positiveTest and negativeTests don't actually test Qseqs (the Qseq in the first test case above is there to make sure
    //BCLs are preferred over Qseqs)

    @Test(dataProvider="positiveTestData")
    public void positiveTests(final IlluminaFileUtil.SupportedIlluminaFormat[] formats,
                              final IlluminaDataType[] dataTypes,
                              final int lane,
                              final List<Integer> tiles,
                              final int[] cycles,
                              final String readStructure) {
        makeFiles(formats, lane, tiles, cycles);
        writeTileMetricsOutFile(makeMap(makeList(lane-1, lane + 1, lane),
                                        makeList(makeList(1,2,3), tiles, tiles)));

        String [] args = makeCheckerArgs(basecallDir, lane, readStructure, dataTypes);
        int result = new CheckIlluminaDirectory().instanceMain(args);
        Assert.assertEquals(result, 0);
    }

    @DataProvider(name="negativeTestData")
    public Object[][] negativeTestData() {
        return new Object[][] {
            { //Completely missing data types
                new SupportedIlluminaFormat[]{Bcl, Filter}, new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF, IlluminaDataType.Position, IlluminaDataType.RawIntensities},
                new ArrayList<String>(),
                new ArrayList<String>(),
                2, makeList(1101,1201,1301, 2101,2201,2301), IlluminaFileUtilTest.cycleRange(1,152), "68T8B68T",
                2
            },
            {
                new SupportedIlluminaFormat[]{Bcl, Filter}, new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF},
                makeList("BaseCalls/L002/C13.1/s_2_1201.bcl", "BaseCalls/L002/C13.1/s_2_2101.bcl"),
                makeList("BaseCalls/L002/s_2_2101.filter"),
                2, makeList(1101,1201,1301, 2101,2201,2301), IlluminaFileUtilTest.cycleRange(1,20), "13T",
                3
            },
            {
                new SupportedIlluminaFormat[]{Bcl, Filter}, new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF},
                new ArrayList<String>(),
                new ArrayList<String>(),
                5, makeList(1101,1201,1301, 2101,2201,2301), IlluminaFileUtilTest.cycleRange(1,152), "250T",
                98
            },
        };
    }

    @Test(dataProvider="negativeTestData")
    public void negativeTests(final IlluminaFileUtil.SupportedIlluminaFormat[] formats,
                              final IlluminaDataType[] dataTypes,
                              final List<String> filesToDelete,
                              final List<String> filesToEmpty,
                              final int lane,
                              final List<Integer> tiles,
                              final int[] cycles,
                              final String readStructure,
                              final int expectedNumErrors) {
        makeFiles(formats, lane, tiles, cycles);
        IlluminaFileUtilTest.deleteRelativeFiles(intensityDir, filesToDelete);
        IlluminaFileUtilTest.emptyRelativeFiles(intensityDir,  filesToEmpty);
        writeTileMetricsOutFile(makeMap(makeList(lane - 1, lane + 1, lane), makeList(makeList(1, 2, 3), tiles, tiles)));

        String [] args = makeCheckerArgs(basecallDir, lane, readStructure, dataTypes);
        int result = new CheckIlluminaDirectory().instanceMain(args);
        Assert.assertEquals(expectedNumErrors, result);
    }

    public void writeFileOfSize(final File file, final int size) {
        try {
            final BufferedWriter writer = new BufferedWriter(new FileWriter(file));
            for(int i = 0; i < size; i++) {
                int toWrite = Math.min(1000, size);
                char [] writeBuffer = new char[toWrite];
                for(int j = 0; j < writeBuffer.length; j++) {
                    writeBuffer[j] = (char)(Math.random() * 150);
                }

                writer.write(writeBuffer);
            }
            writer.flush();
            writer.close();
        } catch(Exception exc) {
            throw new RuntimeException(exc);
        }
    }

    @Test
    public void differentSizedBclTest() {
        final SupportedIlluminaFormat [] formats = new SupportedIlluminaFormat[]{Bcl, Filter};
        final int lane = 5;
        final List<Integer> tiles = makeList(1,2,3,4);
        final int [] cycles = IlluminaFileUtilTest.cycleRange(1, 50);
        final IlluminaDataType [] dataTypes = new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores};

        makeFiles(new SupportedIlluminaFormat[]{Bcl, Filter}, lane, tiles, cycles);
        writeTileMetricsOutFile(makeMap(makeList(lane-1, lane + 1, lane),
                                        makeList(makeList(1,2,3), tiles, tiles)));

        final File cycleDir = new File(basecallDir, "L005/C9.1");
        writeFileOfSize(new File(cycleDir, "s_5_3.bcl"), 222);

        String [] args = makeCheckerArgs(basecallDir, lane, "50T", dataTypes);
        int result = new CheckIlluminaDirectory().instanceMain(args);
        Assert.assertEquals(1, result);
    }

    @Test(expectedExceptions = PicardException.class)
    public void basedirDoesntExistTest() {
        String [] args = makeCheckerArgs(new File("a_made_up_file/in_some_weird_location"), 1, "76T76T", new IlluminaDataType[]{IlluminaDataType.Position});
        
        final int result = new CheckIlluminaDirectory().instanceMain(args);
    }

    @Test
    public void qseqTest() {
        final IlluminaDataType [] dataTypes = new IlluminaDataType[]{BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF, IlluminaDataType.Position};
        final int lane = 4;
        final List<Integer> tiles = makeList(1,2,3);

        IoUtil.copyDirectoryTree(new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir"), basecallDir);
        writeTileMetricsOutFile(makeMap(makeList(lane), makeList(tiles)));

        String [] args = makeCheckerArgs(basecallDir, lane, "76T76T", dataTypes);
        int result = new CheckIlluminaDirectory().instanceMain(args);
        Assert.assertEquals(result, 0);

        args = makeCheckerArgs(basecallDir, lane, "76T77T", dataTypes);
        result = new CheckIlluminaDirectory().instanceMain(args);
        Assert.assertEquals(result, 1);

        IlluminaFileUtilTest.deleteRelativeFiles(basecallDir, makeList("s_4_1_0002_qseq.txt"));

        args = makeCheckerArgs(basecallDir, lane, "76T76T", dataTypes);
        result = new CheckIlluminaDirectory().instanceMain(args);
        Assert.assertEquals(result, 1);
    }
}
