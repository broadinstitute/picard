package picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

public class TileMetricsUtilTest {
    @Test
    public void testFindTileMetrics() throws IOException {
        List<Integer> populatedCycleDirs = new ArrayList<>();
        Path baseDir = createTestDirs(0, populatedCycleDirs, true);
        List<File> tileMetricsFiles = TileMetricsUtil.renderTileMetricsFilesFromBasecallingDirectory(
                baseDir.toFile(), 0, false);

        List<File> expected = generateFindTileMetricsExpected(baseDir, populatedCycleDirs, true);
        Assert.assertEquals(tileMetricsFiles, expected);
        IOUtil.deleteDirectoryTree(baseDir.toFile());
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMissingTileMetrics() throws IOException {
        List<Integer> populatedCycleDirs = new ArrayList<>();
        Path baseDir = createTestDirs(0, populatedCycleDirs, false);
        TileMetricsUtil.renderTileMetricsFilesFromBasecallingDirectory(
                baseDir.toFile(), 0, false);

        IOUtil.deleteDirectoryTree(baseDir.toFile());
    }

    // put tile metrics file in several permutations of random cycle dirs
    // and make sure we get all metricsFiles each time
    @Test
    public void testFindTileMetricsNovaSeq() throws IOException {
        int numCycles = 15;
        int numRandomDirectoryPermutations = 5;
        for (int i=0; i<numRandomDirectoryPermutations; i++) {
            List<Integer> populatedCycleDirs = randomPopulatedCycleDirs(numCycles);
            testPermutation(populatedCycleDirs, numCycles, true);
            testPermutation(populatedCycleDirs, numCycles, false);
        }
    }

    private void testPermutation(
            List<Integer> populatedCycleDirs,
            int numCycles,
            boolean baseMetricsFile
    ) throws IOException {

        Path baseDir = createTestDirs(numCycles, populatedCycleDirs, baseMetricsFile);
        List<File> tileMetricsFiles = TileMetricsUtil.renderTileMetricsFilesFromBasecallingDirectory(
            baseDir.toFile(), numCycles, true);

        List<File> expected = generateFindTileMetricsExpected(baseDir, populatedCycleDirs, baseMetricsFile);
        Assert.assertEquals(tileMetricsFiles, expected);
        IOUtil.deleteDirectoryTree(baseDir.toFile());
    }

    private List<Integer> randomPopulatedCycleDirs(int maxDirNum) {
        List<Integer> populatedCycleDirs = new ArrayList<>();
        Random rng = new Random();
        IntStream.range(1, maxDirNum).forEach(i -> {
            if (rng.nextBoolean()) {
                populatedCycleDirs.add(i);
            }
        });
        return populatedCycleDirs;
    }

    private Path getBaseOpDir(Path baseDir) {
        return baseDir.resolve(TileMetricsUtil.INTEROP_SUBDIRECTORY_NAME);
    }

    private Path createTestDirs(
            int numCycleMetricsFiles,
            List<Integer> populatedCycleDirs,
            boolean baseMetricsFile
    ) throws IOException {
        Path baseDir = Files.createTempDirectory("TileMetricsUtilTest");
        // create numCylcleMetricsFiles cycle dirs put tile metrics in those indicated by populatedCycleDirs
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
            expected.add(new File(
                baseDir
                    .resolve(TileMetricsUtil.INTEROP_SUBDIRECTORY_NAME)
                    .resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME).toUri())
            );
        }
        if (populatedCycleDirs.size() > 0) {
            Path baseOpDir = getBaseOpDir(baseDir);
            populatedCycleDirs.stream().sorted(Comparator.reverseOrder()).forEach(i ->
                expected.add(new File(
                    baseOpDir
                        .resolve("C" + i + ".1")
                        .resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME)
                        .toUri()
                ))
            );
        }
        return expected;
    }
}
