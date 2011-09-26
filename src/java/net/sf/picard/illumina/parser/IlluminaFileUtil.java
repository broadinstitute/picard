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

import net.sf.picard.util.MathUtil;
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
 * @author alecw@broadinstitute.org
 */
public class IlluminaFileUtil {


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
        return getEndedIlluminaBasecallFiles(directory, fileType, lane, end).length > 0;
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
    public static TiledIlluminaFile[] getEndedIlluminaBasecallFiles(final File directory, final String fileType, final int lane, final int end) {
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
    public static TiledIlluminaFile[] getEndedIlluminaBasecallFiles(final File directory, final String fileType,
                                                                    final int lane, final int end, final List<Integer> tiles) {
        if (tiles == null) {
            return getEndedIlluminaBasecallFiles(directory, fileType, lane, end);
        }
        final TiledIlluminaFile[] ret = new TiledIlluminaFile[tiles.size()];
        for (int i = 0; i < tiles.size(); ++i) {
            final String filename = "s_" + lane + "_" + end + "_" + tileNumberFormatter.get().format(tiles.get(i)) + "_" + fileType + ".txt";
            File f = new File(directory, filename);
            if (!f.exists()) {
                f = new File(directory, filename + ".gz");
                if (!f.exists()) {
                    throw new IlluminaFileNotFoundException(f, "Requested tile file " + f + " not found, with or without .gz.");
                }
            }
            ret[i] = new TiledIlluminaFile(f, tiles.get(i));

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
    static TiledIlluminaFile[] getTiledIlluminaBasecallFiles(final File directory, final String regexp) {
        final Pattern pattern = Pattern.compile(regexp);
        final File[] files = getNonEmptyFilesMatchingRegexp(directory, pattern);
        final TiledIlluminaFile[] ret = new TiledIlluminaFile[files.length];
        for (int i = 0; i < files.length; ++i) {
            final Matcher m = pattern.matcher(files[i].getName());
            if (!m.matches()) {
                throw new PicardException(files[i].getName() + " does not match " + m.pattern());
            }
            final String tileString = m.group(1);
            final int tile = Integer.parseInt(tileString);
            ret[i] = new TiledIlluminaFile(files[i], tile);
        }
        Arrays.sort(ret, new TiledIlluminaFileComparator());
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
        return getNonEndedIlluminaBasecallFiles(directory, fileType, lane).length > 0;
    }

    /**
     * Get list of basecall files of given type that are not distinguished by end
     * @param directory where to look for files
     * @param fileType e.g. "barcode.txt", "qseq" etc
     * @param lane
     * @return list of files of the given type, sorted in ascending tile # (numeric, not lexical).
     */
    public static TiledIlluminaFile[] getNonEndedIlluminaBasecallFiles(final File directory, final String fileType, final int lane) {
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
    public static TiledIlluminaFile[] getNonEndedIlluminaBasecallFiles(final File directory, final String fileType,
                                                                       final int lane, final List<Integer> tiles) {
        if (tiles == null) {
            return getNonEndedIlluminaBasecallFiles(directory, fileType, lane);
        }
        final TiledIlluminaFile[] ret = new TiledIlluminaFile[tiles.size()];
        for (int i = 0; i < tiles.size(); ++i) {
            final String filename = makeNonEndedIlluminaBasecallFilename(fileType, lane, tiles.get(i));
            File f = new File(directory, filename);
            if (!f.exists()) {
                f = new File(directory, filename + ".gz");
                if (!f.exists()) {
                    throw new IlluminaFileNotFoundException(f, "Requested tile file " + f + " not found, with or without .gz.");
                }
            }
            ret[i] = new TiledIlluminaFile(f, tiles.get(i));

        }
        return ret;
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

    /**
     * Sort TiledIlluminaFiles numerically by tile number.
     */
    static class TiledIlluminaFileComparator implements Comparator<TiledIlluminaFile> {
        @Override
        public int compare(final TiledIlluminaFile v1, final TiledIlluminaFile v2) {
            return MathUtil.compare(v1.tile, v2.tile);
        }
    }
}
