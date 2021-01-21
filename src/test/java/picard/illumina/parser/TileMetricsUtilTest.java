package picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.IntStream;

public class TileMetricsUtilTest {

    private Path testDir;

    @BeforeClass
    public void setupTestDir() throws IOException {
        this.testDir = Files.createTempDirectory("TileMetricsUtilTests");
    }

    @AfterClass
    public void cleanupTestDir() {
        IOUtil.deleteDirectoryTree(testDir.toFile());
    }

    @DataProvider
    public Object[][] testCases() {
        List<Integer> oneCycleFile = Collections.singletonList(1);
        List<Integer> allCycleFiles = Arrays.asList(1, 2, 3, 4, 5);
        List<Integer> someCycleFiles = Arrays.asList(2, 4);

        return new Object[][]{
            //specificTestDir, numCycleMetricsFiles, populatedCycleDirs, baseMetricsFile
            {"notNovaSeq", 0, Collections.emptyList(), true},
            {"novaSeqBaseFileOnly", 0, Collections.emptyList(), true},
            {"novaSeqSingleCycleFile", 1, oneCycleFile, false},
            {"novaSeqBaseAndSingleCycleFile", 1, oneCycleFile, true},
            {"novaSeqAllCycleFiles", 5, allCycleFiles, false},
            {"novaSeqBaseAndAllCycleFiles", 5, allCycleFiles, true},
            {"novaSeqSomeCycleFiles", 5, someCycleFiles, false},
            {"novaSeqBaseAndSomeCycleFiles", 5, someCycleFiles, true}
        };
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMissingTileMetrics() throws IOException {
        List<Integer> populatedCycleDirs = Collections.emptyList();
        Path baseDir = createTestDirs("missingTileMetrics", 0, populatedCycleDirs, false);
        TileMetricsUtil.findTileMetricsFiles(
                baseDir.toFile(), 0);
    }

    @Test(dataProvider = "testCases")
    public void testFilePresencePermutation(
            String specificTestDir,
            int numCycles,
            List<Integer> populatedCycleDirs,
            boolean baseMetricsFile
    ) throws IOException {

        Path baseDir = createTestDirs(specificTestDir, numCycles, populatedCycleDirs, baseMetricsFile);
        List<File> tileMetricsFiles = TileMetricsUtil.findTileMetricsFiles(
            baseDir.toFile(), numCycles);

        List<File> expected = generateFindTileMetricsExpected(baseDir, populatedCycleDirs, baseMetricsFile);
        Assert.assertEquals(tileMetricsFiles, expected);
    }

    private Path getBaseOpDir(Path baseDir) {
        return baseDir.resolve(TileMetricsUtil.INTEROP_SUBDIRECTORY_NAME);
    }

    private Path createTestDirs(
            String specificTestDir,
            int numCycleMetricsFiles,
            List<Integer> populatedCycleDirs,
            boolean baseMetricsFile
    ) throws IOException {
        Path baseDir = Files.createDirectory(testDir.resolve(specificTestDir));
        // create numCycleMetricsFiles cycle dirs put tile metrics in those indicated by populatedCycleDirs
        Path baseOpDir = Files.createDirectory(getBaseOpDir(baseDir));
        IntStream.range(1, numCycleMetricsFiles + 1).forEach(i -> {
            try {
                Path cycleDir = Files.createDirectory(baseOpDir.resolve("C" + i + ".1"));
                if (populatedCycleDirs.contains(i)) {
                    Files.createFile(cycleDir.resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME));
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
        if (baseMetricsFile) {
            Files.createFile(baseOpDir.resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME));
        }
        return baseDir;
    }

    private List<File> generateFindTileMetricsExpected(Path baseDir, List<Integer> populatedCycleDirs, boolean baseMetricsFile) {
        List<File> expected = new ArrayList<>();
        if (baseMetricsFile) {
            expected.add(
                baseDir
                    .resolve(TileMetricsUtil.INTEROP_SUBDIRECTORY_NAME)
                    .resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME).toFile()
            );
        }
        if (populatedCycleDirs.size() > 0) {
            Path baseOpDir = getBaseOpDir(baseDir);
            populatedCycleDirs.stream().sorted(Comparator.reverseOrder()).forEach(i ->
                expected.add(
                    baseOpDir
                        .resolve("C" + i + ".1")
                        .resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME)
                        .toFile()
                )
            );
        }
        return expected;
    }
}
