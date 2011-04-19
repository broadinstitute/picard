/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.picard.illumina.parser;

import net.sf.picard.io.IoUtil;
import net.sf.picard.PicardException;
import net.sf.picard.util.Log;

import java.io.File;
import java.util.Iterator;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * An iterator of iterators for cluster intensity format files, used to store raw intensities or noise.
 * The outer iteration is over all the tiles in the directory, in order.  The inner iteration is over
 * the cycle files for the given tile and file type.
 *
 * @author alecw@broadinstitute.org
 */
class CycleFileSetIterator implements Iterator<Iterator<TiledIlluminaFile>> {
    private static Log log = Log.getInstance(CycleFileSetIterator.class);
    // Regex for cycle subdirectories
    private static final Pattern CYCLE_SUBDIRECTORY_PATTERN = Pattern.compile("C(\\d+)\\..+");

    private final File directory;
    private final int lane;
    private final ClusterIntensityFileReader.FileType fileType;
    // cycle directories, in cycle order.
    private final File[] cycleDirs;
    // List of files for the first cycle, in order in which they should be iterated over.
    private final TiledIlluminaFile[] firstCycleFiles;
    // Index of next tile to iterate over.
    private int tileIndex = 0;


    /**
     * Prepare to iterate over all tiles for this lane.
     * @param directory Typically a subdirectory of the form Data/Intensities/L00<lane> .
     * @param lane
     * @param fileType Either cif or cnf.
     * @param numCycles The expected number of cycles.  If the run died prematurely, there may be more cycle files
     * than this, and the remainder are ignored.
     */
    CycleFileSetIterator(final File directory, final int lane, final ClusterIntensityFileReader.FileType fileType,
                         final int numCycles) {
        this(directory, lane, fileType, numCycles, null);
    }

    /**
     * Prepare to iterate
     * @param directory Typically a subdirectory of the form Data/Intensities/L00<lane>
     * @param lane
     * @param fileType Either cif or cnf.
     * @param numCycles The expected number of cycles.  If the run died prematurely, there may be more cycle files
     * than this, and the remainder are ignored.
     * @param tiles Requested tiles.  Iteration over tiles is in this order.  If null, all tiles are iterated,
     * in numeric order.  @throws PicardException If one of the requested tiles does not exist.
     */
    CycleFileSetIterator(final File directory, final int lane, final ClusterIntensityFileReader.FileType fileType,
                         final int numCycles, final List<Integer> tiles) {
        this.directory = directory;
        this.fileType = fileType;
        this.lane = lane;
        final File[] tempCycleDirs = IoUtil.getFilesMatchingRegexp(this.directory, CYCLE_SUBDIRECTORY_PATTERN);
        if (tempCycleDirs == null || tempCycleDirs.length == 0) {
            throw new IlluminaFileNotFoundException(directory, directory + " has no cycle subdirectories");
        }
        final CycleDirectoryComparator comparator = new CycleDirectoryComparator();
        Arrays.sort(tempCycleDirs, comparator);
        for (int i = 0; i < tempCycleDirs.length - 1; ++i) {
            if (comparator.compare(tempCycleDirs[i], tempCycleDirs[i+1]) == 0) {
                throw new PicardException("More than one directory with same cycle number: " + tempCycleDirs[i] +
                "; " + tempCycleDirs[i+1]);
            }
        }
        if (tempCycleDirs.length < numCycles) {
            throw new PicardException("Requested more cycle directories (" + numCycles + ") than available (" +
            tempCycleDirs.length + ") for " + directory);
        } else if (tempCycleDirs.length > numCycles) {
            log.warn(directory + " contains " + tempCycleDirs.length + " cycle directories, but only " +
                    numCycles + " cycles were requested.");
            cycleDirs = Arrays.copyOf(tempCycleDirs, numCycles);
        } else {
            cycleDirs = tempCycleDirs;
        }
        if (tiles == null) {
            final String regexp = "s_" + this.lane + "_(\\d+)." + this.fileType;
            firstCycleFiles = IlluminaFileUtil.getTiledIlluminaBasecallFiles(cycleDirs[0], regexp);
        } else {
            firstCycleFiles = new TiledIlluminaFile[tiles.size()];
            for (int i = 0; i < tiles.size(); ++i) {
                final File f = new File(cycleDirs[0], "s_" + this.lane + "_" + tiles.get(i) + "." + this.fileType);
                if (!f.exists()) {
                    throw new IlluminaFileNotFoundException(f, "Requested tile not found");
                }
                firstCycleFiles[i] = new TiledIlluminaFile(f, tiles.get(i));
            }
        }
    }

    /**
     * Move iteration to the given tile number.
     * @param tileNumber
     */
    void seekToTile(final int tileNumber) {
        for (tileIndex = 0; tileIndex < firstCycleFiles.length; ++tileIndex) {
            if (firstCycleFiles[tileIndex].tile == tileNumber) {
                break;
            }
        }
    }

    int getNumberOfCycleFiles() {
        return cycleDirs.length;
    }

    @Override
    public boolean hasNext() {
        return tileIndex < firstCycleFiles.length;
    }

    @Override
    public Iterator<TiledIlluminaFile> next() {
        return new TileCycleFileIterator(tileIndex++);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Iterator over all the cycle files for a tile, in cycle order.
     */
    class TileCycleFileIterator implements Iterator<TiledIlluminaFile> {
        private final int tileIndex;
        private int cycleIndex = 0;

        TileCycleFileIterator(final int tileIndex) {
            this.tileIndex = tileIndex;
        }

        @Override
        public boolean hasNext() {
            return cycleIndex < cycleDirs.length;
        }

        @Override
        public TiledIlluminaFile next() {
            final File f = new File(cycleDirs[cycleIndex++], firstCycleFiles[tileIndex].file.getName());
            return new TiledIlluminaFile(f, firstCycleFiles[tileIndex].tile);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    private static class CycleDirectoryComparator implements Comparator<File> {

        @Override
        public int compare(final File file1, final File file2) {
            Matcher matcher = CYCLE_SUBDIRECTORY_PATTERN.matcher(file1.getName());
            if (!matcher.matches()) {
                throw new PicardException("unpossible");
            }
            final int cycle1 = Integer.parseInt(matcher.group(1));
            matcher = CYCLE_SUBDIRECTORY_PATTERN.matcher(file2.getName());
            if (!matcher.matches()) {
                throw new PicardException("unpossible");
            }
            final int cycle2 = Integer.parseInt(matcher.group(1));
            return cycle1 - cycle2;
        }
    }

}
