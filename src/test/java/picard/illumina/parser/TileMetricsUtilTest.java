package picard.illumina.parser;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.IntStream;

public class TileMetricsUtilTest {
    @Test
    public void testFindTileMetrics() throws IOException {
        Path baseDir = createTestDirs(0, true);
        File tileMetricsFile = TileMetricsUtil.renderTileMetricsFileFromBasecallingDirectory(
                baseDir.toFile(), 0, false);

        Assert.assertEquals(tileMetricsFile.getAbsolutePath(),
                baseDir.resolve(TileMetricsUtil.INTEROP_SUBDIRECTORY_NAME)
                        .resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME).toString());

        IOUtil.deleteDirectoryTree(baseDir.toFile());
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMissingTileMetrics() throws IOException {
        Path baseDir = createTestDirs(0, false);
        TileMetricsUtil.renderTileMetricsFileFromBasecallingDirectory(
                baseDir.toFile(), 0, false);

        IOUtil.deleteDirectoryTree(baseDir.toFile());
    }

    // put tile metrics file in 1 random cycle dir and interop based ensure we get the one in base
    @Test
    public void testFindTileMetricsNovaSeqInBaseInterOp() throws IOException {
        Path baseDir = createTestDirs(1, true);
        File tileMetricsFile = TileMetricsUtil.renderTileMetricsFileFromBasecallingDirectory(
                baseDir.toFile(), 1, true);

        Assert.assertEquals(tileMetricsFile.getAbsolutePath(),
                baseDir.resolve(TileMetricsUtil.INTEROP_SUBDIRECTORY_NAME)
                        .resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME).toString());

        IOUtil.deleteDirectoryTree(baseDir.toFile());
    }

    // put a file named tile metrics in 2 cycle directories and ensure that we get the one in the highest cycle
    @Test
    public void testFindTileMetricsNovaSeqInCycleDirs() throws IOException {
        Path baseDir = createTestDirs(2, false);
        File tileMetricsFile = TileMetricsUtil.renderTileMetricsFileFromBasecallingDirectory(
                baseDir.toFile(), 2, true);

        Assert.assertEquals(tileMetricsFile.getAbsolutePath(),
                baseDir.resolve(TileMetricsUtil.INTEROP_SUBDIRECTORY_NAME).resolve("C2.1")
                        .resolve(TileMetricsUtil.TILE_METRICS_OUT_FILE_NAME).toString());

        IOUtil.deleteDirectoryTree(baseDir.toFile());
    }

    private Path createTestDirs(int numCycleMetricsFiles, boolean baseMetricsFile) throws IOException {
        Path baseDir = Files.createTempDirectory("TileMetricsUtilTest");
        // create 10 cycle dirs put tile metrics in the last one.
        Path baseOpDir = Files.createDirectory(baseDir.resolve(TileMetricsUtil.INTEROP_SUBDIRECTORY_NAME));
        IntStream.range(1, numCycleMetricsFiles + 1).forEach(i -> {
            try {
                Path cycleDir = Files.createDirectory(baseOpDir.resolve("C" + i + ".1"));
                if (i == numCycleMetricsFiles) {
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
}
