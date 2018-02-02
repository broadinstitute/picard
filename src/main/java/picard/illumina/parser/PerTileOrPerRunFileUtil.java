package picard.illumina.parser;

import picard.illumina.parser.fakers.FileFaker;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Pattern;

public class PerTileOrPerRunFileUtil extends PerTileFileUtil {
    private final File runFile;

    public PerTileOrPerRunFileUtil(String extension, File base, FileFaker faker, int lane) {
        super(extension, base, faker, lane);
        Pattern runFileMatchPattern = Pattern.compile("^s" + extension + "$");
        runFile = getRunFile(base.getParentFile(), runFileMatchPattern);
    }

    @Override
    public boolean filesAvailable() {
        return super.filesAvailable() || runFile != null;
    }

    @Override
    public void setTilesForPerRunFile(List<Integer> tiles) {
        if (runFile != null) {
            tiles.forEach(i -> fileMap.put(i, runFile));
            this.tiles = tiles;
        }
    }

    @Override
    public boolean checkTileCount() {
        return runFile == null;
    }

    @Override
    public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
        // if we have a runFile the base should be the parent directory of that file and not
        // the lane base directory of a per tile file.
        final File baseDir = runFile != null ? runFile.getParentFile() : this.base;
        return verifyPerTile(baseDir, expectedTiles);
    }
}
