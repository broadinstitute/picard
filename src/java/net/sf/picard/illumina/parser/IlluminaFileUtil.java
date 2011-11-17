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

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Misc utilities for parsing Illumina output files.
 *
 * @author jburke@broadinstitute.org
 */
public class IlluminaFileUtil {
    private static final Pattern CYCLE_SUBDIRECTORY_PATTERN = Pattern.compile("C(\\d+)\\..+");

    /**
     * For creating tile number strings that are 4 digits wide with zero padding.
     */
    private static final ThreadLocal<NumberFormat> tileNumberFormatter = new ThreadLocal<NumberFormat>() {
        @Override
        protected NumberFormat initialValue() {
            final NumberFormat formatter = NumberFormat.getNumberInstance();
            formatter.setMinimumIntegerDigits(4);
            formatter.setGroupingUsed(false);
            return formatter;
        }
    };

    /**
     * See if files exist for a particular file type, lane and end number.
     * @param directory Where to look for the files.
     * @param fileType (e.g. "qseq")
     * @param lane
     * @param end Should be 1, 2 or 3.
     * @return true if any files matching the specifications exist.
     */
    public static boolean endedIlluminaBasecallFilesExist(final File directory, final String fileType, final int lane, final int end) {
        return getEndedIlluminaBasecallFiles(directory, fileType, lane, end).size() > 0;
    }

     /**
     * Find Illumina files that are end-specific.  Finds both .txt and .txt.gz files.
     *
     * @param directory Where to look for the files.
     * @param fileType (e.g. "qseq")
     * @param lane
     * @return searches for a list of files of the form s_<lane>_<end>_<tile>_<fileType>.txt(.gz)? and finds the highest
     * <end> (where end is a positive integer).
     */
    public static int getNumberOfIlluminaEnds(final File directory, final String fileType, final int lane) {
        final String regexp = "s_" + lane + "_(\\d+)_\\d{4}_" + fileType + ".txt(.gz)?";
        final Pattern pattern = Pattern.compile(regexp);

        final File[] files = getNonEmptyFilesMatchingRegexp(directory, pattern);

        int highestMatch = 0;
        for (int i = 0; i < files.length; ++i) {
            final Matcher m = pattern.matcher(files[i].getName());
            m.matches();
            final String tileString = m.group(1);
            final int end = Integer.parseInt(tileString);
            if(end > highestMatch)
                highestMatch = end;
        }

        return highestMatch;
    }

    /**
     * Find Illumina files that are end-specific.  Finds both .txt and .txt.gz files.
     *
     * @param directory Where to look for the files.
     * @param fileType (e.g. "qseq")
     * @param lane
     * @param end Should be 1, 2 or 3.
     * @return a list of files of the form s_<lane>_<end>_<tile>_<fileType>.txt(.gz)?, sorted by tile number, along with the
     * tile number for each file.
     */
    public static IlluminaFileMap getEndedIlluminaBasecallFiles(final File directory, final String fileType, final int lane, final int end) {
        final String regexp = "s_" + lane + "_" + end + "_(\\d{4})_" + fileType + ".txt(.gz)?";
        return getTiledIlluminaBasecallFiles(directory, regexp);
    }

    /**
     * Find Illumina files that are end-specific.  Finds both .txt and .txt.gz files
     * @param directory Where to look for the files.
     * @param fileType (e.g. "qseq")
     * @param lane
     * @param end Should be either 1 or 2.
     * @param tiles Requested tiles.
     * @return a list of files of the form s_<lane>_<end>_<tile>_<fileType>.txt(.gz)?, in the order specified by tiles parameter,
     * along with the tile number for each file.
     * @throws IlluminaFileNotFoundException If a requested tile does not exist.
     */
    public static IlluminaFileMap getEndedIlluminaBasecallFiles(final File directory, final String fileType,
                                                                    final int lane, final int end, final List<Integer> tiles) {
        if (tiles == null) {
            return getEndedIlluminaBasecallFiles(directory, fileType, lane, end);
        }
        final IlluminaFileMap ret = new IlluminaFileMap();
        for (int i = 0; i < tiles.size(); ++i) {
            final String filename = "s_" + lane + "_" + end + "_" + tileNumberFormatter.get().format(tiles.get(i)) + "_" + fileType + ".txt";
            File f = new File(directory, filename);
            if (!f.exists()) {
                f = new File(directory, filename + ".gz");
                if (!f.exists()) {
                    throw new IlluminaFileNotFoundException(f, "Requested tile file " + f + " not found, with or without .gz.");
                }
            }
            ret.put(tiles.get(i), f);

        }
        return ret;
    }

    /**
     * Search the directory for files matching regexp, and return sorted list of them.
     *
     * @param directory Where to find the files.
     * @param regexp Regular expression for matching the files.  group #1 must encapsulate the tile #.
     * @return Array of TiledIlluminaFile matching regexp, sorted in ascending tile # (numeric, not lexical).
     */
    static IlluminaFileMap getTiledIlluminaBasecallFiles(final File directory, final String regexp) {
        final Pattern pattern = Pattern.compile(regexp);
        final File[] files = getNonEmptyFilesMatchingRegexp(directory, pattern);
        final IlluminaFileMap ret = new IlluminaFileMap();
        for (int i = 0; i < files.length; ++i) {
            final Matcher m = pattern.matcher(files[i].getName());
            m.matches();
            final String tileString = m.group(1);
            final int tile = Integer.parseInt(tileString);
            ret.put(tile, files[i]);
        }

        return ret;
    }

    private static File[] getNonEmptyFilesMatchingRegexp(final File directory, final Pattern pattern) {
        final File[] files = IoUtil.getFilesMatchingRegexp(directory, pattern);
        final ArrayList<File> nonEmptyFiles = new ArrayList<File>(files.length);
        for (final File file : files) {
            if (file.length() > 0) {
                nonEmptyFiles.add(file);
            }
        }
        return nonEmptyFiles.toArray(new File[nonEmptyFiles.size()]);
    }

    /**
     * Check if a particular non-ended file type exists for the given lane.
     * @param directory where to look for files
     * @param fileType e.g. "barcode.txt", "qseq", etc
     * @param lane
     * @return true if any files of the given type exist.
     */
    public static boolean nonEndedIlluminaBasecallFilesExist(final File directory, final String fileType, final int lane) {
        return getNonEndedIlluminaBasecallFiles(directory, fileType, lane).size() > 0;
    }

    /**
     * Get list of basecall files of given type that are not distinguished by end
     * @param directory where to look for files
     * @param fileType e.g. "barcode.txt", "qseq" etc
     * @param lane
     * @return list of files of the given type, sorted in ascending tile # (numeric, not lexical).
     */
    public static IlluminaFileMap getNonEndedIlluminaBasecallFiles(final File directory, final String fileType, final int lane) {
        final String regexp = "s_" + lane + "_(\\d{4})_" + fileType + ".txt(.gz)?";
        return getTiledIlluminaBasecallFiles(directory, regexp);
    }

    /**
     * Get list of basecall files of given type that are not distinguished by end.
     *
     * @param directory Where to look for files.
     * @param fileType e.g. "barcode.txt", "qseq" etc.
     * @param lane
     * @param tiles List of desired tiles.  A file must exist for each requested tile.
     * @return List of TiledIlluminaFiles, in the order of the tiles argument.  Each file may or may not have .gz extension.
     * @throws IlluminaFileNotFoundException if a file is not found for a tile.
     */
    public static IlluminaFileMap getNonEndedIlluminaBasecallFiles(final File directory, final String fileType,
                                                                       final int lane, final List<Integer> tiles) {
        if (tiles == null) {
            return getNonEndedIlluminaBasecallFiles(directory, fileType, lane);
        }
        final IlluminaFileMap ret = new IlluminaFileMap();
        for (int i = 0; i < tiles.size(); ++i) {
            final String filename = makeNonEndedIlluminaBasecallFilename(fileType, lane, tiles.get(i));
            File f = new File(directory, filename);
            if (!f.exists()) {
                f = new File(directory, filename + ".gz");
                if (!f.exists()) {
                    throw new IlluminaFileNotFoundException(f, "Requested tile file " + f + " not found, with or without .gz.");
                }
            }
            ret.put(tiles.get(i), f);

        }
        return ret;
    }

    /**
     * For the given directory/L00<Lane>/C<X>/s_<lane>_<tile>.<fileExt> return a Map Tile -> CycledFilesIterator where there
     * will be one X for each cycle and one tile for every value in tiles.  If tiles == null, use all available tiles.
     * @param directory basecalls directory
     * @param fileExt file extension of the files to look for
     * @param lane Lane to iterate over
     * @param tiles Tiles to to use, if null use all tiles
     * @param expectedCycles Number of files expected per tile
     * @return Map (Tile -> CycledFilesIterator)
     */
    public static CycleIlluminaFileMap getCyledIlluminaFiles(final File directory, final String fileExt, final int lane, final List<Integer> tiles, int expectedCycles) {

        final File[] tempCycleDirs = IoUtil.getFilesMatchingRegexp(directory, CYCLE_SUBDIRECTORY_PATTERN);

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

        final CycleIlluminaFileMap cycledMap = new CycleIlluminaFileMap();
        if (tiles == null) {
            final String regexp = "s_" + lane + "_(\\d+)." + fileExt;
            final IlluminaFileMap tileMap = IlluminaFileUtil.getTiledIlluminaBasecallFiles(tempCycleDirs[0], regexp);
            for(int i = 0; i < tileMap.lastKey(); i++) {
                cycledMap.put(i, new CycleFilesIterator(directory, lane, i, fileExt));
            }
            cycledMap.assertValid(new ArrayList<Integer>(tileMap.keySet()), expectedCycles);
        } else {
            for(int i = 0; i < tiles.size(); i++) {
                cycledMap.put(i, new CycleFilesIterator(directory, lane, i, fileExt));
            }
            cycledMap.assertValid(tiles, expectedCycles);
        }

        return cycledMap;
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
    /**
     * Create a non-end-specifid Illumina-style filename
     * @param fileType E.g. "seq" or "barcode"
     * @param lane
     * @param tile
     * @return String of the form "s_<lane>_<tile>_<fileType>.txt"  It may be necessary to append ".gz" to find a file.
     */
    public static String makeNonEndedIlluminaBasecallFilename(final String fileType, final int lane, final int tile) {
        return "s_" + lane + "_" + tileNumberFormatter.get().format(tile) + "_" + fileType + ".txt";
    }
}
