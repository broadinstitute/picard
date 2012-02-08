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
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * General utils for dealing with IlluminaFiles as well as utils for specific, support formats.
 * This class contains utils that span across multiple Illumina files but it's primary intent
 * was to provide support for basic file types.  Each supported file type can be accessed
 * via a factory method (make<filetype>Ft).  When IlluminaFileUtil is created it is parameterized
 * by basecallDir and lane and all IlluminaFileTypes created by IlluminaFileUtil will also be
 * parameterized in this fashion.
 *
 * @author jburke@broadinstitute.org
 */
class IlluminaFileUtil {

    public enum SupportedIlluminaFormat {
        Qseq,
        Bcl,
        Cif,
        Cnf,
        Locs,
        Clocs,
        Pos,
        Filter,
        Barcode
    }

    public final File intensityDir;
    public final File intensityLaneDir;
    public final File basecallDir;
    public final File basecallLaneDir;
    public final int lane;


    /** A regex string matching only qseq files */
    private final String QSeqRegex = "s_(\\d+)_(\\d)_(\\d{4}).qseq.txt";   //s_lane_end_tile
    private final QSeqIlluminaFileUtil qseq;
    private final PerTilePerCycleFileUtil bcl;
    private final PerTilePerCycleFileUtil cif;
    private final PerTilePerCycleFileUtil cnf;
    private final PerTileFileUtil pos;
    private final PerTileFileUtil locs;
    private final PerTileFileUtil clocs;
    private final PerTileFileUtil filter;
    private final PerTileFileUtil barcode;
    private final Map<SupportedIlluminaFormat, ParameterizedFileUtil> utils;

    public IlluminaFileUtil(final File basecallDir, final int lane) {
        this.basecallDir  = basecallDir;
        this.intensityDir = basecallDir.getParentFile();
        this.lane = lane;

        this.basecallLaneDir = new File(basecallDir, longLaneStr(lane));
        this.intensityLaneDir = new File(intensityDir, longLaneStr(lane));

        utils = new HashMap<SupportedIlluminaFormat, ParameterizedFileUtil>();

        qseq = new QSeqIlluminaFileUtil();
        utils.put(SupportedIlluminaFormat.Qseq, qseq);

        bcl  = new PerTilePerCycleFileUtil(".bcl", basecallLaneDir);
        utils.put(SupportedIlluminaFormat.Bcl, bcl);

        cif  = new PerTilePerCycleFileUtil(".cif");
        utils.put(SupportedIlluminaFormat.Cif, cif);

        cnf  = new PerTilePerCycleFileUtil(".cnf");
        utils.put(SupportedIlluminaFormat.Cnf, cnf);

        locs    = new PerTileFileUtil(".locs",    false);
        utils.put(SupportedIlluminaFormat.Locs, locs);

        clocs   = new PerTileFileUtil(".clocs",   false);
        utils.put(SupportedIlluminaFormat.Clocs, clocs);

        pos     = new PerTileFileUtil("_pos.txt", false, intensityDir);
        utils.put(SupportedIlluminaFormat.Pos,  pos);

        filter  = new PerTileFileUtil(".filter",  true, basecallLaneDir);
        utils.put(SupportedIlluminaFormat.Filter, filter);

        barcode = new PerTileFileUtil("_barcode.txt", true, basecallDir);
        utils.put(SupportedIlluminaFormat.Barcode, barcode);
    }

    /** Return the lane we're inspecting */
    public int getLane() {
        return lane;
    }

    /** Given a file type, get the Parameterized File Util object associated with it*/
    public ParameterizedFileUtil getUtil(final SupportedIlluminaFormat format) {
        return utils.get(format);
    }

    /** Get the available tiles for the given formats, if the formats have tile lists that differ then
     * throw an exception, if any of the format
     */
    public List<Integer> getTiles(final List<SupportedIlluminaFormat> formats) {
        if(formats == null) {
            throw new PicardException("Format list provided to getTiles was null!");
        }

        if(formats.size() < 0) {
            throw new PicardException("0 Formats were specified.  You need to specify at least SupportedIlluminaFormat to use getTiles");
        }

        final List<Integer> tiles = utils.get(formats.get(0)).getTiles();
        for(int i = 0; i < formats.size(); i++) {
            final List<Integer> fmTiles = utils.get(formats.get(i)).getTiles();
            if(tiles.size() != fmTiles.size() || !tiles.containsAll(fmTiles)) {
                throw new PicardException("Formats do not have the same number of tiles! " + summarizeTileCounts(formats));
            }
        }

        return tiles;
    }

    public QSeqIlluminaFileUtil qseq() {
        return qseq;
    }

    public PerTilePerCycleFileUtil bcl() {
       return bcl;
    }

    public PerTilePerCycleFileUtil cif() {
        return cif;
    }

    public PerTilePerCycleFileUtil cnf() {
        return cnf;
    }

    public PerTileFileUtil locs() {
        return locs;
    }
    public PerTileFileUtil clocs() {
        return clocs;
    }

    public PerTileFileUtil pos() {
        return pos;
    }

    public PerTileFileUtil filter() {
        return filter;
    }

    public PerTileFileUtil barcode() {
        return barcode;
    }

    public static String UNPARAMETERIZED_PER_TILE_PATTERN = "s_(\\d+)_(\\d{1,4})";
    public static String UNPARAMETERIZED_QSEQ_PATTERN     = "s_(\\d+)_(\\d)_(\\d{4})_qseq\\.txt(\\.gz|\\.bz2)?";
    private static final Pattern CYCLE_SUBDIRECTORY_PATTERN = Pattern.compile("^C(\\d+)\\.1$");

    public static String makeParameterizedLaneAndTileRegex(final int lane) {
        if(lane < 0) {
            throw new PicardException("Lane (" + lane + ") cannot be negative");
        }
        return "s_" + lane + "_(\\d{1,4})";
    }

    public static String makeParameterizedQseqRegex(final int lane) {
        if(lane < 0) {
            throw new PicardException("Lane (" + lane + ") cannot be negative");
        }
        return "s_" + lane + "_(\\d)_(\\d{4})_qseq\\.txt(\\.gz|\\.bz2)?";
    }


    /** An object providing utilities for locating Illumina files of specific types */
    abstract class ParameterizedFileUtil {
        /** The file extension for this class, file extension does not have the standard meaning
         * in this instance.  It means, all the characters that come after the identifying portion of
         * the file (after lane, tile, and end that is).  So _qseq.txt and .filter are both file extensions.
         */
        public final String extension;

        /** A pattern that will match files of this type for this lane*/
        public final Pattern pattern;

        /** A pattern that will match files of this type for this lane*/
        public final Pattern unparameterizedPattern;

        /** If you think of the file system as a tree, this is the deepest directory(node) on the tree that
         * still contains all of the files for this given type (e.g. If we're talking about BCLs the directory
         * structure is:
         *
         * BaseCall Dir
         * |
         *      L001
         * |     |        |
         * C1.1 C2.1 ... Cn.1
         * |     |        |
         * bcl Files ... bclFiles
         *
         * L001 is the base because it contains every BCL file in the run (though those files are nested in
         * other folders).
         */
        protected final File base;

        public ParameterizedFileUtil(final String unparameterizedPattern, final String patternStr, final String extension, final File base) {
            this.pattern                  = Pattern.compile(escapePeriods(patternStr));
            this.unparameterizedPattern   = Pattern.compile(escapePeriods(unparameterizedPattern));
            this.extension = extension;
            this.base      = base;
        }

        /** The period separator is expectexd in the file extension, since some do not start with it */
        private String escapePeriods(final String preEscaped) {
            return preEscaped.replaceAll("\\.", "\\."); //In the first one the \\ is inside a regex in the second it's NOT
        }

        /**
         * Determine whether or not files are available
         * @return return true if files are found matching this types pattern, false otherwise
         */
        public abstract boolean filesAvailable();

        /**
         * Illumina file names contain at least lane and tile information and sometimes end info. Return all
         * available lane tile and end information.
         * @param fileName Filename to analyze for data
         * @return A LaneTileEnd object with discovered values or null if that value is not available in the given file name
         */
        public abstract LaneTileEnd fileToLaneTileEnd(final String fileName);

        /**
         * Return a list of all tiles available for this file format and run
         * @return A List of tile integers
         */
        public abstract List<Integer> getTiles();
    }

    /** Represents file types that have one file per tile */
    class PerTileFileUtil extends ParameterizedFileUtil {
        protected final boolean txtBased;
        protected final boolean padTile;
        protected final IlluminaFileMap fileMap;
        protected final List<Integer> tiles;

        public PerTileFileUtil(final String fileNameEndPattern, boolean padTile, final File base) {
            super(makeLTRegex(processTxtExtension(fileNameEndPattern)), makeLTRegex(processTxtExtension(fileNameEndPattern), lane), fileNameEndPattern, base);
            this.txtBased = fileNameEndPattern.endsWith(".txt");
            this.padTile  = padTile;
            this.fileMap  = getTiledFiles(base, pattern, this);

            if(fileMap.size() > 0) {
                this.tiles = Collections.unmodifiableList(new ArrayList<Integer>(this.fileMap.keySet()));
            } else {
                this.tiles = new ArrayList<Integer>();
            }
        }

        public PerTileFileUtil(final String fileNameEndPattern, boolean padTile) {
            this(fileNameEndPattern, padTile, intensityLaneDir);
        }

        @Override
        public boolean filesAvailable() {
            return !fileMap.isEmpty();
        }

        /**
         * Returns only lane and tile information as PerTileFt's do not have End information.
         * @param fileName Filename to analyze for data
         * @return A LaneTileEnd object with the discovered Lane and Tile information and a null end field.
         */
        public LaneTileEnd fileToLaneTileEnd(final String fileName) {
            return laneAndTileFromFirstTwoMatches(fileName, unparameterizedPattern);
        }

        public IlluminaFileMap getFiles() {
            return fileMap;
        }

        public IlluminaFileMap getFiles(final List<Integer> tiles) {
            return fileMap.keep(tiles);
        }

        public List<Integer> getTiles() {
            return tiles;
        }
    }

    /**
     * A base class for file types that occur 1 for each tile/cycle.
     */
    class PerTilePerCycleFileUtil extends ParameterizedFileUtil {
        private final CycleIlluminaFileMap cycleFileMap;
        private int numCycles;
        private List<Integer> tiles;
        public PerTilePerCycleFileUtil(final String fileNameEndPattern, final File base) {
            super(makeLTRegex(fileNameEndPattern), makeLTRegex(fileNameEndPattern, lane), fileNameEndPattern, base);
            this.cycleFileMap = getPerTilePerCycleFiles(); //sideEffect, assigned to numCycles

            if(cycleFileMap.size() > 0) {
                this.tiles = Collections.unmodifiableList(new ArrayList<Integer>(this.cycleFileMap.keySet()));
            } else {
                this.tiles = new ArrayList<Integer>();
            }
        }

        public PerTilePerCycleFileUtil(final String fileNameEndPattern) {
            this(fileNameEndPattern, intensityLaneDir);
        }

        /**
         * Returns only lane and tile information as PerTilePerCycleFt's do not have End information.
         * @param fileName Filename to analyze for data
         * @return A LaneTileEnd object with the discovered Lane and Tile information and a null end field.
         */
        public LaneTileEnd fileToLaneTileEnd(final String fileName) {
            return laneAndTileFromFirstTwoMatches(fileName, unparameterizedPattern);
        }

        /**
         * Given a cycle directory, return a list of tiles in that directory.  If expectedTiles equals null
         * return all files discovered otherwise filter by expectedTiles.
         * @param cycleDir The file object of the cycle directory we are searching
         * @return A list of tile integers describing the tiles available in a cycle directory
         */
        private List<Integer> getTilesInCycleDir(final File cycleDir) {
            final File [] files = IoUtil.getFilesMatchingRegexp(cycleDir, pattern);
            final List<Integer> tiles = new ArrayList<Integer>();
            for(final File file : files) {
                if(file.length() > 0) {
                    tiles.add(fileToLaneTileEnd(file.getName()).tile);
                }
            }

            return tiles;
        }

        /**
         * For the given tiles, populate a CycleIlluminaFileMap that contains all these tiles and will iterate through
         * all the files for these tiles in expectedBase
         * Side Effect: Assigns numCycles
         * @return A CycleIlluminaFileMap with the listed (or all) tiles for at least expectedCycles number of cycles(or total available
         * cycles if expectedCycles is null)
         */
        private CycleIlluminaFileMap getPerTilePerCycleFiles() {
            final CycleIlluminaFileMap cycledMap = new CycleIlluminaFileMap();

            final File laneDir = base;
            final File[] tempCycleDirs = IoUtil.getFilesMatchingRegexp(laneDir, CYCLE_SUBDIRECTORY_PATTERN);

            if (tempCycleDirs == null || tempCycleDirs.length == 0) {
                return cycledMap;
            }

            File cycle1Dir = null;
            final Integer [] cycles = new Integer[tempCycleDirs.length];
            for (int i = 0; i < tempCycleDirs.length; ++i) {
                cycles[i] = getCycleFromDir(tempCycleDirs[i]);
                if(cycles[i] == 1) {
                    cycle1Dir = tempCycleDirs[i];
                }
            }

            Arrays.sort(cycles);
            for (int i = 0; i < tempCycleDirs.length; ++i) {
                if(cycles[i] != i + 1) {
                    StringBuilder presentCycles = new StringBuilder(tempCycleDirs.length);
                    for(int j = 0; j < tempCycleDirs.length; j++) {
                        presentCycles.append(cycles[j]);
                        presentCycles.append(",");
                    }
                    throw new PicardException("Missing cycle dir " + (i+1) + presentCycles.toString());
                }
            }

            numCycles = tempCycleDirs.length;
            final List<Integer> tiles = getTilesInCycleDir(cycle1Dir);
            for(final int tile : tiles) {
                cycledMap.put(tile, new CycleFilesIterator(laneDir, lane, tile, extension)); //Gonna have a problem here if we ever get a (.txt.gz for these types of files)
            }

            return cycledMap;
        }

        public int getNumCycles() {
            return numCycles;
        }

        public CycleIlluminaFileMap getFiles() {
            return cycleFileMap;
        }

        public CycleIlluminaFileMap getFiles(final List<Integer> tiles) {
            return cycleFileMap.keep(tiles);
        }

        /**
         * Discover all files of this type in expectedBase that match pattern and construct a list of tiles
         * available based on these files.  The same number of tiles is expected in each cycle dir.
         * @return A list of tile integers for all tiles available
         */
        public List<Integer> getTiles() {
            return tiles;
        }

        public boolean filesAvailable() {
            return !cycleFileMap.isEmpty();
        }
    }

    /** QSeq files are really tiled and ended so define it's own nested format since no other file types
     * are structured the same. */
    class QSeqIlluminaFileUtil extends ParameterizedFileUtil {
        private final List<Integer> tiles;
        private final List<IlluminaFileMap> fileMaps;
        public QSeqIlluminaFileUtil() {
            super(UNPARAMETERIZED_QSEQ_PATTERN, makeParameterizedQseqRegex(lane), "_qseq.txt", basecallDir);
            fileMaps = getFiles();

            if(fileMaps.size() > 0) {
                tiles = Collections.unmodifiableList(new ArrayList<Integer>(fileMaps.get(0).keySet()));
            } else {
                tiles = new ArrayList<Integer>();
            }
        }

        /** Make a qSeq regex string with the lane and end already filled in */
        private String makeLaneAndEndSpecificRegex(final int lane, final int end) {
            return "^s_" + lane + "_" + end + "_\\d{4}_qseq\\.txt(\\.gz|\\.bz2)?$";
        }

        /**
         * Return the number of ends found in the basecallDir
         * @return The highest end number found among the files in the basecallDir
         */
        public int numberOfEnds() {
            return fileMaps.size();
        }

        /**
         * Given a file name return it's Lane, Tile, and End information
         * @param fileName The name of a file to analyze
         * @return The lane, tile, and end of the file with the given name
         */
        @Override
        public LaneTileEnd fileToLaneTileEnd(final String fileName) {
            final Matcher matcher = unparameterizedPattern.matcher(fileName);
            if(!matcher.matches()) {
                return null;
            }
            return new LaneTileEnd(Integer.parseInt(matcher.group(1)), Integer.parseInt(matcher.group(3)), Integer.parseInt(matcher.group(2)));
        }

        /**
         * For each tile in tiles with the given end find the corresponding QSeq file.  Place that qseq file in an IlluminaFileMap
         * and after all tiles are processed, return that fileMap;
         * @param end A single end integer
         * @return A map of tiles->Files where each file is represents the given tile and end
         */
        private IlluminaFileMap getFiles(final int end) {
            final String regex = makeLaneAndEndSpecificRegex(lane, end);
            return getTiledFiles(basecallDir, Pattern.compile(regex), this);
        }

        /**
         * Return a list of illumina file map, where index 0 contains files for end 1, index 1 contains files for end 2, etc...
         * @return An list of illuminaFileMaps with containing all files for all ends for each given tile
         */
        public List<IlluminaFileMap> getFiles() {
            final List<IlluminaFileMap> readTileMap = new ArrayList<IlluminaFileMap>();

            boolean emptyMap = false;
            for(int i = 1; !emptyMap; i++) {
                final IlluminaFileMap fm = getFiles(i);
                if(fm.isEmpty()) {
                    emptyMap = true;
                } else {
                    readTileMap.add(fm);
                }
            }
            return readTileMap;
        }

        public List<IlluminaFileMap> getFiles(final List<Integer> tiles) {
            final List<IlluminaFileMap> filteredMaps = new ArrayList<IlluminaFileMap>();
            for(final IlluminaFileMap fm : fileMaps) {
                filteredMaps.add(fm.keep(tiles));
            }

            return filteredMaps;
        }

        public List<Integer> getTiles() {
            return tiles;
        }

        public boolean filesAvailable() {
            return !tiles.isEmpty();
        }
    };

    /** A support class for return lane tile and end information for a given file */
    static class LaneTileEnd {
        public final Integer lane;
        public final Integer tile;
        public final Integer end;

        public LaneTileEnd(final Integer lane, final Integer tile, final Integer end) {
            this.lane = lane;
            this.tile = tile;
            this.end  = end;
        }

        public LaneTileEnd(final Integer lane, final Integer tile) {
            this(lane, tile, null);
        }
    }

    /** Return a regex string for finding Lane and Tile given a file extension pattern */
    private static String makeLTRegex(final String fileNameEndPattern) {
        return "^" + UNPARAMETERIZED_PER_TILE_PATTERN + fileNameEndPattern + "$";
    }

    /** Return a regex string for finding Lane and Tile given a file extension pattern */
    private static String makeLTRegex(final String fileNameEndPattern, final int lane) {
        return "^" + makeParameterizedLaneAndTileRegex(lane) + fileNameEndPattern + "$";
    }

    private static int getCycleFromDir(File tempCycleDir) {
        final char [] name = tempCycleDir.getName().toCharArray();
        if(name[0] != 'C') {
            throw new PicardException("Invalid cycle directory name " + tempCycleDir.getName());
        }

        String intStr = "";
        boolean periodFound = false;
        for(int i = 1; i < name.length && !periodFound; i++) {
            if(name[i] == '.') {
                periodFound = true;
            } else if(name[i] == '1' || name[i] == '2' || name[i] == '3' ||
               name[i] == '4' || name[i] == '5' || name[i] == '6' ||
               name[i] == '7' || name[i] == '8' || name[i] == '9' ||
               name[i] == '0') {
                intStr += name[i];
            } else {
               throw new PicardException("Invalid cycle directory name " + tempCycleDir.getAbsolutePath());
            }
        }

        return Integer.parseInt(intStr);
    }

    /** Given a pattern and file name return a LaneTileEnd with the first two matches to the pattern returned
     * as the lane and tile respectively */
    private static LaneTileEnd laneAndTileFromFirstTwoMatches(final String fileName, final Pattern pattern) {
        final Matcher matcher = pattern.matcher(fileName);
        if(!matcher.matches()) {
            return null;
        }
        return new LaneTileEnd(Integer.parseInt(matcher.group(1)), Integer.parseInt(matcher.group(2)));
    }

    /**
     * Return a string representing the Lane in the format "L00<lane>"
     * @param lane The lane to transform
     * @return A long string representation of the name
     */
    private static String longLaneStr(final int lane) {
        String lstr = String.valueOf(lane);
        final int zerosToAdd = 3 - lstr.length();

        for(int i = 0; i < zerosToAdd; i++) {
            lstr = "0" + lstr;
        }
        return "L" + lstr;
    }

    /** Return all files that match pattern of the given file type in the given base directory */
    private static IlluminaFileMap getTiledFiles(final File baseDirectory, final Pattern pattern, final ParameterizedFileUtil ift) {
        final IlluminaFileMap fileMap = new IlluminaFileMap();
        if(baseDirectory.exists()) {
            IoUtil.assertDirectoryIsReadable(baseDirectory);
            final File [] files = IoUtil.getFilesMatchingRegexp(baseDirectory, pattern);
            for(final File file : files) {
                if(file.length() > 0) {
                    final LaneTileEnd lt = ift.fileToLaneTileEnd(file.getName());
                    fileMap.put(lt.tile, file);
                }
            }
        }

        return fileMap;
    }

    /** For filename patterns that end with .txt tack on the option .gz extension */
    private static String processTxtExtension(final String fileNameEndPattern) {
        if(fileNameEndPattern.endsWith(".txt")) {
            return fileNameEndPattern + "(\\.gz|\\.bz2)?";
        } else {
            return fileNameEndPattern;
        }
    }


    private String liToStr(final List<Integer> intList) {
        if(intList.size() == 0)
            return "";

        String summary = String.valueOf(intList.get(0));
        for(int i = 1; i < intList.size(); i++) {
            summary += ", " + String.valueOf(intList.get(i));
        }

        return summary;
    }

    private String summarizeTileCounts(final List<SupportedIlluminaFormat> formats) {
        String summary;
        ParameterizedFileUtil pfu = utils.get(formats.get(0));
        List<Integer> tiles = pfu.getTiles();
        summary = pfu.extension + "(" + liToStr(tiles) + ")";

        for(int i = 0; i < formats.size(); i++) {
            pfu = utils.get(formats.get(i));
            tiles = pfu.getTiles();

            summary += ", " + pfu.extension + "(" + liToStr(tiles) + ")";
        }

        return summary;
    }
}
