package picard.illumina.parser;

import picard.illumina.parser.fakers.FileFaker;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class PerTileFileUtil extends ParameterizedFileUtil {
    private final IlluminaFileMap fileMap;

    public PerTileFileUtil(final String extension, final File base,
                           final FileFaker faker, final int lane) {
        this(extension, base, faker, lane, DefaultSkipEmptyFiles);
    }

    public PerTileFileUtil(final String extension, final File base,
        final FileFaker faker, final int lane, final boolean skipEmptyFiles) {
        super(true, extension, base, faker, lane, skipEmptyFiles);
        this.fileMap = getTiledFiles(base, matchPattern);
        if (!fileMap.isEmpty()) {
            this.tiles = new ArrayList<>(this.fileMap.keySet());
        } else {
            this.tiles = new ArrayList<>();
        }
    }

    @Override
    public boolean filesAvailable() {
        return !fileMap.isEmpty();
    }

    public IlluminaFileMap getFiles() {
        return fileMap;
    }

    public IlluminaFileMap getFiles(final List<Integer> tiles) {
        return fileMap.keep(tiles);
    }

    @Override
    public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
        final List<String> failures = new LinkedList<String>();

        if (!base.exists()) {
            failures.add("Base directory(" + base.getAbsolutePath() + ") does not exist!");
        } else {
                if (!tiles.containsAll(expectedTiles)) {
                    final List<Integer> missing = new ArrayList<Integer>(expectedTiles);
                    missing.removeAll(tiles);
                    failures.add("Missing tile " + missing + " for file type " + extension + ".");
                }
        }
        return failures;
    }

    @Override
    public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] cycles,
                                  final IlluminaFileUtil.SupportedIlluminaFormat format) {
        final List<String> failures = new LinkedList<String>();
        if (!base.exists()) {
            failures.add("Base directory(" + base.getAbsolutePath() + ") does not exist!");
        } else {
            for (final Integer tile : expectedTiles) {
                if (!tiles.contains(tile) || fileMap.get(tile).length() == 0) {
                    //create a new file of this type
                    try {
                        faker.fakeFile(base, tile, lane, extension);
                    } catch (final IOException e) {
                        failures.add(String.format("Could not create fake file %s: %s", fileMap.get(tile),
                                e.getMessage()));
                    }

                }
            }
        }
        return failures;
    }

}
