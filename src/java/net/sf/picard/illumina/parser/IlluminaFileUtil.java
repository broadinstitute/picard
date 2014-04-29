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
import net.sf.picard.illumina.parser.fakers.BarcodeFileFaker;
import net.sf.picard.illumina.parser.fakers.BciFileFaker;
import net.sf.picard.illumina.parser.fakers.BclFileFaker;
import net.sf.picard.illumina.parser.fakers.CifFileFaker;
import net.sf.picard.illumina.parser.fakers.ClocsFileFaker;
import net.sf.picard.illumina.parser.fakers.CnfFileFaker;
import net.sf.picard.illumina.parser.fakers.FileFaker;
import net.sf.picard.illumina.parser.fakers.FilterFileFaker;
import net.sf.picard.illumina.parser.fakers.LocsFileFaker;
import net.sf.picard.illumina.parser.fakers.MultiTileBclFileFaker;
import net.sf.picard.illumina.parser.fakers.MultiTileLocsFileFaker;
import net.sf.picard.illumina.parser.fakers.PosFileFaker;
import net.sf.picard.illumina.parser.fakers.QSeqFileFaker;
import net.sf.picard.illumina.parser.readers.BclReader;
import net.sf.picard.illumina.parser.readers.TileMetricsOutReader;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;
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
public class IlluminaFileUtil {

    public enum SupportedIlluminaFormat {
        Qseq,
        Bcl,
        Cif,
        Cnf,
        Locs,
        Clocs,
        Pos,
        Filter,
        Barcode,
        MultiTileFilter,
        MultiTileLocs,
        MultiTileBcl
    }

    private final File intensityLaneDir;
    private final File basecallDir;
    private final int lane;


    /**
     * A regex string matching only qseq files
     */
    private final QSeqIlluminaFileUtil qseq;
    private final PerTilePerCycleFileUtil bcl;
    private final PerTilePerCycleFileUtil cif;
    private final PerTilePerCycleFileUtil cnf;
    private final PerTileFileUtil pos;
    private final PerTileFileUtil locs;
    private final PerTileFileUtil clocs;
    private final PerTileFileUtil filter;
    private final PerTileFileUtil barcode;
    private final MultiTileFilterFileUtil multiTileFilter;
    private final MultiTileLocsFileUtil multiTileLocs;
    private final MultiTileBclFileUtil multiTileBcl;
    private final File tileMetricsOut;
    private final Map<SupportedIlluminaFormat, ParameterizedFileUtil> utils;


    public IlluminaFileUtil(final File basecallDir, final int lane) {
        this.basecallDir = basecallDir;
        final File intensityDir = basecallDir.getParentFile();
        final File dataDir = intensityDir.getParentFile();
        final File interopDir = new File(dataDir.getParentFile(), "InterOp");
        this.lane = lane;

        final File basecallLaneDir = new File(basecallDir, longLaneStr(lane));
        this.intensityLaneDir = new File(intensityDir, longLaneStr(lane));

        utils = new HashMap<SupportedIlluminaFormat, ParameterizedFileUtil>();

        qseq = new QSeqIlluminaFileUtil();
        utils.put(SupportedIlluminaFormat.Qseq, qseq);

        bcl = new PerTilePerCycleFileUtil(inferBclExtension(basecallLaneDir), basecallLaneDir, new BclFileFaker());
        utils.put(SupportedIlluminaFormat.Bcl, bcl);

        cif = new PerTilePerCycleFileUtil(".cif", new CifFileFaker());
        utils.put(SupportedIlluminaFormat.Cif, cif);

        cnf = new PerTilePerCycleFileUtil(".cnf", new CnfFileFaker());
        utils.put(SupportedIlluminaFormat.Cnf, cnf);

        locs = new PerTileFileUtil(".locs", false, new LocsFileFaker());
        utils.put(SupportedIlluminaFormat.Locs, locs);

        clocs = new PerTileFileUtil(".clocs", false, new ClocsFileFaker());
        utils.put(SupportedIlluminaFormat.Clocs, clocs);

        pos = new PerTileFileUtil("_pos.txt", false, intensityDir, new PosFileFaker());
        utils.put(SupportedIlluminaFormat.Pos, pos);

        filter = new PerTileFileUtil(".filter", true, basecallLaneDir, new FilterFileFaker());
        utils.put(SupportedIlluminaFormat.Filter, filter);

        barcode = new PerTileFileUtil("_barcode.txt", true, basecallDir, new BarcodeFileFaker());
        utils.put(SupportedIlluminaFormat.Barcode, barcode);

        multiTileFilter = new MultiTileFilterFileUtil(basecallLaneDir);
        utils.put(SupportedIlluminaFormat.MultiTileFilter, multiTileFilter);

        multiTileLocs = new MultiTileLocsFileUtil(new File(intensityDir, basecallLaneDir.getName()), basecallLaneDir);
        utils.put(SupportedIlluminaFormat.MultiTileLocs, multiTileLocs);

        multiTileBcl = new MultiTileBclFileUtil(basecallLaneDir);
        utils.put(SupportedIlluminaFormat.MultiTileBcl, multiTileBcl);

        tileMetricsOut = new File(interopDir, "TileMetricsOut.bin");
    }

    /**
     * Return the lane we're inspecting
     */
    public int getLane() {
        return lane;
    }

    /**
     * Given a file type, get the Parameterized File Util object associated with it
     */
    public ParameterizedFileUtil getUtil(final SupportedIlluminaFormat format) {
        return utils.get(format);
    }

    /**
     * Return the list of tiles we would expect for this lane based on the metrics found in InterOp/TileMetricsOut.bin
     */
    public List<Integer> getExpectedTiles() {
        IoUtil.assertFileIsReadable(tileMetricsOut);
        //Used just to ensure predictable ordering
        final TreeSet<Integer> expectedTiles = new TreeSet<Integer>();

        final Iterator<TileMetricsOutReader.IlluminaTileMetrics> tileMetrics = new TileMetricsOutReader(tileMetricsOut);
        while (tileMetrics.hasNext()) {
            final TileMetricsOutReader.IlluminaTileMetrics tileMetric = tileMetrics.next();

            if (tileMetric.getLaneNumber() == lane) {
                if (!expectedTiles.contains(tileMetric.getTileNumber())) {
                    expectedTiles.add(tileMetric.getTileNumber());
                }
            }
        }

        CloserUtil.close(tileMetrics);
        return new ArrayList<Integer>(expectedTiles);
    }

    /**
     * Get the available tiles for the given formats, if the formats have tile lists that differ then
     * throw an exception, if any of the format
     */
    public List<Integer> getActualTiles(final List<SupportedIlluminaFormat> formats) {
        if (formats == null) {
            throw new PicardException("Format list provided to getTiles was null!");
        }

        if (formats.size() == 0) {
            throw new PicardException(
                    "0 Formats were specified.  You need to specify at least SupportedIlluminaFormat to use getTiles");
        }

        final List<Integer> tiles = utils.get(formats.get(0)).getTiles();
        for (int i = 0; i < formats.size(); i++) {
            final List<Integer> fmTiles = utils.get(formats.get(i)).getTiles();
            if (tiles.size() != fmTiles.size() || !tiles.containsAll(fmTiles)) {
                throw new PicardException(
                        "Formats do not have the same number of tiles! " + summarizeTileCounts(formats));
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

    public MultiTileFilterFileUtil multiTileFilter() {
        return multiTileFilter;
    }

    public MultiTileLocsFileUtil multiTileLocs() {
        return multiTileLocs;
    }

    public MultiTileBclFileUtil multiTileBcl() {
        return multiTileBcl;
    }

    public File tileMetricsOut() {
        return tileMetricsOut;
    }

    public static final String UNPARAMETERIZED_PER_TILE_PATTERN = "s_(\\d+)_(\\d{1,5})";
    public static final String UNPARAMETERIZED_QSEQ_PATTERN = "s_(\\d+)_(\\d)_(\\d{4})_qseq\\.txt(\\.gz|\\.bz2)?";
    private static final Pattern CYCLE_SUBDIRECTORY_PATTERN = Pattern.compile("^C(\\d+)\\.1$");

    public static String makeParameterizedLaneAndTileRegex(final int lane) {
        if (lane < 0) {
            throw new PicardException("Lane (" + lane + ") cannot be negative");
        }
        return "s_" + lane + "_(\\d{1,5})";
    }

    public static String makeParameterizedQseqRegex(final int lane) {
        if (lane < 0) {
            throw new PicardException("Lane (" + lane + ") cannot be negative");
        }
        return "s_" + lane + "_(\\d)_(\\d{4})_qseq\\.txt(\\.gz|\\.bz2)?";
    }

    /**
     * An object providing utilities for locating Illumina files of specific types
     */
    public abstract class ParameterizedFileUtil {
        /**
         * The file extension for this class, file extension does not have the standard meaning
         * in this instance.  It means, all the characters that come after the identifying portion of
         * the file (after lane, tile, and end that is).  So _qseq.txt and .filter are both file extensions
         */
        public final String extension;

        /**
         * A pattern that will match files of this type for this lane
         */
        public final Pattern pattern;

        /**
         * A pattern that will match files of this type for this lane
         */
        public final Pattern unparameterizedPattern;

        /**
         * If you think of the file system as a tree, this is the deepest directory(node) on the tree that
         * still contains all of the files for this given type (e.g. If we're talking about BCLs the directory
         * structure is:
         * <p/>
         * BaseCall Dir
         * |
         * L001
         * |     |        |
         * C1.1 C2.1 ... Cn.1
         * |     |        |
         * bcl Files ... bclFiles
         * <p/>
         * L001 is the base because it contains every BCL file in the run (though those files are nested in
         * other folders).
         */
        protected final File base;

        protected final FileFaker faker;

        public ParameterizedFileUtil(final String unparameterizedPattern, final String patternStr,
                                     final String extension, final File base,
                                     final FileFaker faker) {
            this.pattern = Pattern.compile(escapePeriods(patternStr));
            this.unparameterizedPattern = Pattern.compile(escapePeriods(unparameterizedPattern));
            this.extension = extension;
            this.base = base;
            this.faker = faker;
        }

        /**
         * The period separator is expected in the file extension, since some do not start with it
         */
        private String escapePeriods(final String preEscaped) {
            return preEscaped
                    .replaceAll("\\.", "\\."); //In the first one the \\ is inside a regex in the second it's NOT
        }

        /**
         * Determine whether or not files are available
         *
         * @return return true if files are found matching this types pattern, false otherwise
         */
        public abstract boolean filesAvailable();

        /**
         * Illumina file names contain at least lane and tile information and sometimes end info. Return all
         * available lane tile and end information.
         *
         * @param fileName Filename to analyze for data
         *
         * @return A LaneTileEnd object with discovered values or null if that value is not available in the given file name
         */
        public abstract LaneTileEnd fileToLaneTileEnd(final String fileName);

        /**
         * Return a list of all tiles available for this file format and run
         *
         * @return A List of tile integers
         */
        public abstract List<Integer> getTiles();


        /**
         * Given the expected tiles/expected cycles for this file type, return a list of error messages describing any
         * missing/or malformed files
         *
         * @param expectedTiles  An ordered list of tile numbers
         * @param expectedCycles An ordered list of cycle numbers that may contain gaps
         *
         * @return A list of error messages for this format
         */
        public abstract List<String> verify(List<Integer> expectedTiles, int[] expectedCycles);

        /**
         * Given the expected tiles/expected cycles for this file type create a set of fake files such that the
         * verification criteria are met.
         *
         * @param expectedTiles An ordered list of tile numbers
         * @param cycles        An ordered list of cycle numbers that may contain gaps
         * @param format        The format of the files that are to be faked
         *
         * @return A list of error messages for this format
         */
        public abstract List<String> fakeFiles(List<Integer> expectedTiles, int[] cycles,
                                               SupportedIlluminaFormat format);
    }

    /**
     * Represents file types that have one file per tile
     */
    class PerTileFileUtil extends ParameterizedFileUtil {
        protected final boolean txtBased;
        protected final boolean padTile;
        protected final IlluminaFileMap fileMap;
        protected final List<Integer> tiles;

        public PerTileFileUtil(final String fileNameEndPattern, final boolean padTile, final File base,
                               final FileFaker fileFaker) {
            super(makeLTRegex(processTxtExtension(fileNameEndPattern)),
                    makeLTRegex(processTxtExtension(fileNameEndPattern), lane), fileNameEndPattern, base,
                    fileFaker);
            this.txtBased = fileNameEndPattern.endsWith(".txt");
            this.padTile = padTile;
            this.fileMap = getTiledFiles(base, pattern, this);

            if (fileMap.size() > 0) {
                this.tiles = Collections.unmodifiableList(new ArrayList<Integer>(this.fileMap.keySet()));
            } else {
                this.tiles = new ArrayList<Integer>();
            }
        }

        public PerTileFileUtil(final String fileNameEndPattern, final boolean padTile, final FileFaker fileFaker) {
            this(fileNameEndPattern, padTile, intensityLaneDir, fileFaker);
        }

        @Override
        public boolean filesAvailable() {
            return !fileMap.isEmpty();
        }

        /**
         * Returns only lane and tile information as PerTileFt's do not have End information.
         *
         * @param fileName Filename to analyze for data
         *
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

        @Override
        public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
            final List<String> failures = new LinkedList<String>();

            if (!base.exists()) {
                failures.add("Base directory(" + base.getAbsolutePath() + ") does not exist!");
            } else {
                for (final Integer tile : expectedTiles) {
                    if (!tiles.contains(tile)) {
                        failures.add("Missing tile " + tile + " for file type " + extension + ".");
                    } else if (fileMap.get(tile).length() == 0) {
                        failures.add("Tile " + tile + " is empty for file type " + extension + ".");
                    }
                }
            }

            return failures;
        }

        @Override
        public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] cycles,
                                      final SupportedIlluminaFormat format) {
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

    /**
     * A base class for file types that occur 1 for each tile/cycle.
     */
    class PerTilePerCycleFileUtil extends ParameterizedFileUtil {
        private final CycleIlluminaFileMap cycleFileMap;
        private final List<Integer> tiles;
        private int[] detectedCycles;

        public PerTilePerCycleFileUtil(final String fileNameEndPattern, final File base, final FileFaker fileFaker) {
            super(makeLTRegex(fileNameEndPattern), makeLTRegex(fileNameEndPattern, lane), fileNameEndPattern, base,
                    fileFaker);
            this.cycleFileMap = getPerTilePerCycleFiles(); //sideEffect, assigned to numCycles

            if (cycleFileMap.size() > 0) {
                this.tiles = Collections.unmodifiableList(new ArrayList<Integer>(this.cycleFileMap.keySet()));
            } else {
                this.tiles = new ArrayList<Integer>();
            }
        }

        public PerTilePerCycleFileUtil(final String fileNameEndPattern, final FileFaker fileFaker) {
            this(fileNameEndPattern, intensityLaneDir, fileFaker);
        }

        /**
         * Returns only lane and tile information as PerTilePerCycleFt's do not have End information.
         *
         * @param fileName Filename to analyze for data
         *
         * @return A LaneTileEnd object with the discovered Lane and Tile information and a null end field.
         */
        public LaneTileEnd fileToLaneTileEnd(final String fileName) {
            return laneAndTileFromFirstTwoMatches(fileName, unparameterizedPattern);
        }

        /**
         * Given a cycle directory, return a list of tiles in that directory.  If expectedTiles equals null
         * return all files discovered otherwise filter by expectedTiles.
         *
         * @param cycleDir The file object of the cycle directory we are searching
         *
         * @return A list of tile integers describing the tiles available in a cycle directory
         */
        private List<Integer> getTilesInCycleDir(final File cycleDir) {
            final File[] files = IoUtil.getFilesMatchingRegexp(cycleDir, pattern);
            final List<Integer> tiles = new ArrayList<Integer>();
            for (final File file : files) {
                if (file.length() > 0) {
                    tiles.add(fileToLaneTileEnd(file.getName()).tile);
                }
            }

            return tiles;
        }

        /**
         * For the given tiles, populate a CycleIlluminaFileMap that contains all these tiles and will iterate through
         * all the files for these tiles in expectedBase
         * Side Effect: Assigns numCycles
         *
         * @return A CycleIlluminaFileMap with the listed (or all) tiles for at least expectedCycles number of cycles(or total available
         * cycles if expectedCycles is null)
         */
        private CycleIlluminaFileMap getPerTilePerCycleFiles() {
            final CycleIlluminaFileMap cycledMap = new CycleIlluminaFileMap();

            final File laneDir = base;
            final File[] tempCycleDirs;
            tempCycleDirs = IoUtil.getFilesMatchingRegexp(laneDir, CYCLE_SUBDIRECTORY_PATTERN);
            if (tempCycleDirs == null || tempCycleDirs.length == 0) {
                return cycledMap;
            }

            int lowestCycle = Integer.MAX_VALUE;
            int lowestCycleDirIndex = 0;
            final int[] cycles = new int[tempCycleDirs.length];
            for (int i = 0; i < tempCycleDirs.length; ++i) {
                cycles[i] = getCycleFromDir(tempCycleDirs[i]);
                if (cycles[i] < lowestCycle) {
                    lowestCycle = cycles[i];
                    lowestCycleDirIndex = i;
                }
            }

            final File firstCycleDir = tempCycleDirs[lowestCycleDirIndex];

            Arrays.sort(cycles);
            detectedCycles = cycles;

            final List<Integer> tiles = getTilesInCycleDir(firstCycleDir);
            for (final int tile : tiles) {
                cycledMap.put(tile, new CycleFilesIterator(laneDir, lane, tile, cycles,
                        extension)); //Gonna have a problem here if we ever get a (.txt.gz for these types of files)
            }

            return cycledMap;
        }

        public CycleIlluminaFileMap getFiles() {
            return cycleFileMap;
        }

        public CycleIlluminaFileMap getFiles(final List<Integer> tiles) {
            return cycleFileMap.keep(tiles, null);
        }

        /**
         * Returns a cycleIlluminaFileMap with all available tiles but limited to the cycles passed in.  Any cycles that are missing
         * cycle files or directories will be removed from the cycle list that is kept.
         *
         * @param cycles Cycles that should be present in the output CycleIlluminaFileMap
         *
         * @return A CycleIlluminaFileMap with all available tiles but at most the cycles passed in by the cycles parameter
         */
        public CycleIlluminaFileMap getFiles(final int[] cycles) {
            //Remove any cycles that were discovered to be NON-EXISTENT when this util was instantiated
            final int[] filteredCycles = removeNonExistentCycles(cycles);
            return cycleFileMap.keep(null, filteredCycles);
        }

        /**
         * Returns a cycleIlluminaFileMap that contains only the tiles and cycles specified (and fewer if the orginal CycleIlluminaFileMap, created
         * on util instantiation, doesn't contain any of these tiles/cycles).
         *
         * @param cycles Cycles that should be present in the output CycleIlluminaFileMap
         *
         * @return A CycleIlluminaFileMap with at most the tiles/cycles listed in the parameters
         */
        public CycleIlluminaFileMap getFiles(final List<Integer> tiles, final int[] cycles) {
            //Remove any cycles that were discovered to be NON-EXISTENT when this util was instantiated
            final int[] filteredCycles = removeNonExistentCycles(cycles);
            return cycleFileMap.keep(tiles, filteredCycles);
        }

        private int[] removeNonExistentCycles(final int[] cycles) {
            final TreeSet<Integer> detectedCyclesSet = new TreeSet<Integer>();
            for (final Integer cycle : detectedCycles) {
                detectedCyclesSet.add(cycle);
            }

            final TreeSet<Integer> inputCyclesSet = new TreeSet<Integer>();
            for (final Integer inputCycle : cycles) {
                inputCyclesSet.add(inputCycle);
            }

            //This also sorts outputCycles
            final int[] outputCycles;
            inputCyclesSet.retainAll(detectedCyclesSet);
            outputCycles = new int[inputCyclesSet.size()];
            int i = 0;
            for (final Integer element : inputCyclesSet) {
                outputCycles[i++] = element;
            }

            return outputCycles;
        }

        public int[] getDetectedCycles() {
            return detectedCycles;
        }

        /**
         * Discover all files of this type in expectedBase that match pattern and construct a list of tiles
         * available based on these files.  The same number of tiles is expected in each cycle dir.
         *
         * @return A list of tile integers for all tiles available
         */
        public List<Integer> getTiles() {
            return tiles;
        }

        public boolean filesAvailable() {
            return !cycleFileMap.isEmpty();
        }

        @Override
        public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
            final List<String> failures = new LinkedList<String>();

            if (!base.exists()) {
                failures.add("Base directory(" + base.getAbsolutePath() + ") does not exist!");
            } else {
                final CycleIlluminaFileMap cfm = getFiles(expectedTiles, expectedCycles);

                final Set<Integer> detectedCycleSet = new HashSet<Integer>();
                for (final Integer cycle : detectedCycles) {
                    detectedCycleSet.add(cycle);
                }

                final Set<Integer> missingCycleSet = new TreeSet<Integer>();
                for (final Integer cycle : expectedCycles) {
                    missingCycleSet.add(cycle);
                }

                missingCycleSet.removeAll(detectedCycleSet);

                for (final Integer tile : expectedTiles) {
                    final CycleFilesIterator cfIterator = cfm.get(tile);
                    if (cfIterator == null) {
                        failures.add("File type " + extension + " is missing tile " + tile);
                    } else if (!cfIterator.hasNext()) {
                        failures.add("File type " + extension + " has 0 cycle files for tile " + tile);
                    } else {
                        int expectedCycleIndex = 0;
                        Long cycleSize = null;

                        while (cfIterator.hasNext() && expectedCycleIndex < expectedCycles.length) {
                            final int currentCycle = expectedCycles[expectedCycleIndex];

                            if (cfIterator.getNextCycle() == currentCycle) {
                                final File cycleFile = cfIterator.next();

                                if (!missingCycleSet.contains(currentCycle)) {
                                    if (!cycleFile.exists()) {
                                        failures.add("Missing file(" + cycleFile.getAbsolutePath() + ")");
                                    } else if (cycleFile.length() == 0) {
                                        failures.add("0 Length tile file(" + cycleFile.getAbsolutePath() + ")");
                                    } else if (cycleSize == null) {
                                        cycleSize = cycleFile.length();
                                    } else if (!extension.equals(".bcl.gz") && cycleSize != cycleFile.length()) {
                                        // TODO: The gzip bcl files might not be the same length despite having the same content,
                                        // for now we're punting on this but this should be looked into at some point
                                        failures.add("File type " + extension
                                                     + " has cycles files of different length.  Current cycle ("
                                                     + currentCycle + ") " +
                                                     "Length of first non-empty file (" + cycleSize
                                                     + ") length of current cycle (" + cycleFile.length() + ")"
                                                     + " File(" + cycleFile.getAbsolutePath() + ")");
                                    }
                                } else {
                                    cfIterator.reset();
                                    throw new PicardException(
                                            "Malformed CycleIlluminaFileMap! CycleIlluminaFileMap has cycle "
                                            + currentCycle
                                            + " even though the directory does not exist!  CycleFileIterator("
                                            + CycleIlluminaFileMap.remainingCyclesToString(cfIterator) + ")");
                                }
                            } else if (!missingCycleSet.contains(currentCycle)) {
                                cfIterator.reset();
                                throw new PicardException(
                                        "Malformed CycleIlluminaFileMap! Tile " + tile + "CycleFileIterator("
                                        + CycleIlluminaFileMap.remainingCyclesToString(cfIterator) + ")");
                            }

                            expectedCycleIndex += 1;
                        }
                    }
                }

                for (final Integer cycle : missingCycleSet) {
                    failures.add("Missing cycle directory " + cycle + " in directory " + base.getAbsolutePath()
                                 + " for file type " + extension);
                }
            }

            return failures;
        }

        @Override
        public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] expectedCycles,
                                      final SupportedIlluminaFormat format) {
            final List<String> failures = new LinkedList<String>();

            if (!base.exists()) {
                base.mkdirs();
            }

            final Set<Integer> detectedCycleSet = new HashSet<Integer>();

            if (detectedCycles == null) {
                detectedCycles = new int[0];
            }

            for (final Integer cycle : detectedCycles) {
                detectedCycleSet.add(cycle);
            }

            final Set<Integer> missingCycleSet = new TreeSet<Integer>();
            for (final Integer cycle : expectedCycles) {
                missingCycleSet.add(cycle);
            }

            missingCycleSet.removeAll(detectedCycleSet);
            for (final Integer cycle : missingCycleSet) {
                final File cycleDirectory = new File(base, "C" + cycle + ".1");
                if (cycleDirectory.mkdirs()) {
                    detectedCycleSet.add(cycle);
                }
            }

            final CycleIlluminaFileMap cfm = getPerTilePerCycleFiles();

            for (final Integer tile : expectedTiles) {
                final CycleFilesIterator cfIterator = cfm.get(tile);
                if (cfIterator == null) {
                    for (final Integer cycle : missingCycleSet) {
                        final File cycleDirectory = new File(base, "C" + cycle + ".1");
                        try {
                            faker.fakeFile(cycleDirectory, tile, lane, extension);
                        } catch (final IOException e) {
                            failures.add(String.format("Could not create fake file %s: %s", tile + extension,
                                    e.getMessage()));
                        }
                    }
                } else if (!cfIterator.hasNext()) {
                    failures.add("File type " + extension + " has 0 cycle files for tile " + tile);
                } else {
                    int expectedCycleIndex = 0;
                    Long cycleSize = null;
                    while (cfIterator.hasNext() && expectedCycleIndex < expectedCycles.length) {
                        final int currentCycle = expectedCycles[expectedCycleIndex];

                        if (cfIterator.getNextCycle() == currentCycle) {
                            final File cycleFile = cfIterator.next();

                            if (cycleSize == null) {
                                // TODO: Assumes that the first file encountered is a valid BCL
                                // TODO: file. This isn't always true.
                                cycleSize = BclReader.getNumberOfClusters(cycleFile);
                            }

                            if (!cycleFile.exists() || cycleFile.length() == 0) {
                                try {
                                    faker.fakeFile(cycleFile, cycleSize.intValue());
                                } catch (final IOException e) {
                                    failures.add("Could not create fake file: " + cycleFile);
                                }
                            }
                        }
                        expectedCycleIndex += 1;
                    }
                }
            }

            for (final Integer cycle : missingCycleSet) {
                failures.add("Missing cycle directory " + cycle + " in directory " + base.getAbsolutePath()
                             + " for file type " + extension);
            }
            return failures;
        }
    }

    /**
     * QSeq files are really tiled and ended so define it's own nested format since no other file types
     * are structured the same.
     */
    class QSeqIlluminaFileUtil extends ParameterizedFileUtil {
        private final List<Integer> tiles;
        private final List<IlluminaFileMap> readFileMaps;

        public QSeqIlluminaFileUtil() {
            super(UNPARAMETERIZED_QSEQ_PATTERN, makeParameterizedQseqRegex(lane), "_qseq.txt", basecallDir,
                    new QSeqFileFaker());
            readFileMaps = getFiles();

            if (readFileMaps.size() > 0) {
                tiles = Collections.unmodifiableList(new ArrayList<Integer>(readFileMaps.get(0).keySet()));
            } else {
                tiles = new ArrayList<Integer>();
            }
        }

        /**
         * Make a qSeq regex string with the lane and end already filled in
         */
        private String makeLaneAndEndSpecificRegex(final int lane, final int end) {
            return "^s_" + lane + "_" + end + "_\\d{4}_qseq\\.txt(\\.gz|\\.bz2)?$";
        }

        /**
         * Return the number of ends found in the basecallDir
         *
         * @return The highest end number found among the files in the basecallDir
         */
        public int numberOfEnds() {
            return readFileMaps.size();
        }

        /**
         * Given a file name return it's Lane, Tile, and End information
         *
         * @param fileName The name of a file to analyze
         *
         * @return The lane, tile, and end of the file with the given name
         */
        @Override
        public LaneTileEnd fileToLaneTileEnd(final String fileName) {
            final Matcher matcher = unparameterizedPattern.matcher(fileName);
            if (!matcher.matches()) {
                return null;
            }
            return new LaneTileEnd(Integer.parseInt(matcher.group(1)), Integer.parseInt(matcher.group(3)),
                    Integer.parseInt(matcher.group(2)));
        }

        /**
         * For each tile in tiles with the given end find the corresponding QSeq file.  Place that qseq file in an IlluminaFileMap
         * and after all tiles are processed, return that fileMap;
         *
         * @param end A single end integer
         *
         * @return A map of tiles->Files where each file is represents the given tile and end
         */
        private IlluminaFileMap getFiles(final int end) {
            final String regex = makeLaneAndEndSpecificRegex(lane, end);
            return getTiledFiles(basecallDir, Pattern.compile(regex), this);
        }

        /**
         * Return a list of illumina file map, where index 0 contains files for end 1, index 1 contains files for end 2, etc...
         *
         * @return An list of illuminaFileMaps with containing all files for all ends for each given tile
         */
        public List<IlluminaFileMap> getFiles() {
            final List<IlluminaFileMap> readTileMap = new ArrayList<IlluminaFileMap>();

            boolean emptyMap = false;
            for (int i = 1; !emptyMap; i++) {
                final IlluminaFileMap fm = getFiles(i);
                if (fm.isEmpty()) {
                    emptyMap = true;
                } else {
                    readTileMap.add(fm);
                }
            }
            return readTileMap;
        }

        public List<IlluminaFileMap> getFiles(final List<Integer> tiles) {
            final List<IlluminaFileMap> filteredMaps = new ArrayList<IlluminaFileMap>();

            for (final IlluminaFileMap fm : readFileMaps) {
                filteredMaps.add(fm.keep(tiles));
            }

            return filteredMaps;
        }

        public List<Integer> getTiles() {
            return tiles;
        }

        @Override
        public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
            final List<String> failures = new LinkedList<String>();

            if (!this.base.exists()) {
                failures.add("Base directory( " + this.base.getAbsolutePath() + ") does not exist!");
            } else {
                final List<IlluminaFileMap> fileMapPerRead = getFiles(expectedTiles);
                final int[] qseqReadLengths = new int[numberOfEnds()];
                int lastCycle = 0;
                for (int i = 0; i < qseqReadLengths.length; i++) {
                    final File currentReadForTile = fileMapPerRead.get(i).get(expectedTiles.get(0));
                    qseqReadLengths[i] = QseqReadParser.getReadLength(currentReadForTile);
                    lastCycle += qseqReadLengths[i];
                }

                final Range cycleRange = new Range(1, lastCycle);
                for (final int expectedCycle : expectedCycles) {
                    if (expectedCycle < cycleRange.start || expectedCycle > cycleRange.end) {
                        failures.add("Expected cycle(" + expectedCycle
                                     + ") is not within the range provided by available qseqs.  " +
                                     "Min Available Cycle(" + cycleRange.start + ") Max Available Cycle("
                                     + cycleRange.end + ") Length of Qseqs( " + StringUtil.join(", ", qseqReadLengths));
                    }
                }

                //ensure that those same ends exist for each expectedTile
                for (int i = 1; i < expectedTiles.size(); i++) {
                    final Integer tile = expectedTiles.get(i);
                    for (int j = 0; j < qseqReadLengths.length; j++) {
                        final File currentReadForTile = fileMapPerRead.get(j).get(tile);
                        if (currentReadForTile == null || !currentReadForTile.exists()) {
                            failures.add("Missing file " + "s_" + lane + "_" + (j + 1) + "_" + longTileStr(tile)
                                         + "_qseq.txt");
                        }
                    }
                }
            }

            return failures;
        }

        @Override
        public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] expectedCycles,
                                      final SupportedIlluminaFormat format) {
            final List<String> failures = new LinkedList<String>();

            if (!this.base.exists()) {
                failures.add("Base directory( " + this.base.getAbsolutePath() + ") does not exist!");
            } else {
                final List<IlluminaFileMap> fileMapPerRead = getFiles(expectedTiles);
                final int[] qseqReadLengths = new int[numberOfEnds()];
                int lastCycle = 0;
                for (int i = 0; i < qseqReadLengths.length; i++) {
                    final File currentReadForTile = fileMapPerRead.get(i).get(expectedTiles.get(0));
                    qseqReadLengths[i] = QseqReadParser.getReadLength(currentReadForTile);
                    lastCycle += qseqReadLengths[i];
                }

                final Range cycleRange = new Range(1, lastCycle);
                for (final int expectedCycle : expectedCycles) {
                    if (expectedCycle < cycleRange.start || expectedCycle > cycleRange.end) {
                        failures.add("Expected cycle(" + expectedCycle
                                     + ") is not within the range provided by available qseqs.  " +
                                     "Min Available Cycle(" + cycleRange.start + ") Max Available Cycle("
                                     + cycleRange.end + ") Length of Qseqs( " + StringUtil.join(", ", qseqReadLengths));
                    }
                }

                //ensure that those same ends exist for each expectedTile
                for (int i = 1; i < expectedTiles.size(); i++) {
                    final Integer tile = expectedTiles.get(i);
                    for (int j = 0; j < qseqReadLengths.length; j++) {
                        final File currentReadForTile = fileMapPerRead.get(j).get(tile);
                        if (currentReadForTile == null || !currentReadForTile.exists()) {
                            failures.add("Missing file " + "s_" + lane + "_" + (j + 1) + "_" + longTileStr(tile)
                                         + "_qseq.txt");
                        }
                    }
                }
            }

            return failures;
        }

        public boolean filesAvailable() {
            return !tiles.isEmpty();
        }
    }

    /**
     * For file types for which there is one file per lane, with fixed record size, and all the tiles in it,
     * so the s_<lane>.bci file can be used to figure out where each tile starts and ends.
     */
    abstract class MultiTileFileUtil<OUTPUT_RECORD extends IlluminaData> extends ParameterizedFileUtil {
        protected final File bci;
        protected TileIndex tileIndex;
        protected File dataFile;

        MultiTileFileUtil(final String extension, final File base, final File bciDir, final FileFaker fileFaker) {
            super(makeLaneRegex(extension), makeLaneRegex(extension, lane), extension, base, fileFaker);
            bci = new File(bciDir, "s_" + lane + ".bci");
            if (bci.exists()) {
                tileIndex = new TileIndex(bci);
            } else {
                tileIndex = null;
            }
            final File[] filesMatchingRegexp = IoUtil.getFilesMatchingRegexp(base, pattern);
            if (filesMatchingRegexp == null || filesMatchingRegexp.length == 0) {
                dataFile = null;
            } else if (filesMatchingRegexp.length == 1) {
                dataFile = filesMatchingRegexp[0];
            } else {
                throw new PicardException("More than one filter file found in " + base.getAbsolutePath());
            }
        }

        @Override
        public boolean filesAvailable() {
            return tileIndex != null && dataFile != null && dataFile.exists();
        }

        @Override
        public LaneTileEnd fileToLaneTileEnd(final String fileName) {
            throw new UnsupportedOperationException();
        }

        @Override
        public List<Integer> getTiles() {
            if (tileIndex == null) {
                return Collections.EMPTY_LIST;
            }
            return tileIndex.getTiles();
        }

        /**
         * expectedCycles are not checked in this implementation.
         */
        @Override
        public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
            if (tileIndex == null) {
                return Collections.singletonList("Tile index(" + bci.getAbsolutePath() + ") does not exist!");
            }
            return tileIndex.verify(expectedTiles);
        }

        @Override
        public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] expectedCycles,
                                      final SupportedIlluminaFormat format) {
            //we need to fake a bci file for the tile index
            final BciFileFaker bciFileFaker = new BciFileFaker();
            try {
                bciFileFaker.fakeBciFile(bci, expectedTiles);
                tileIndex = new TileIndex(bci);
                faker.fakeFile(base, expectedTiles, lane, extension);
                final File[] filesMatchingRegexp = IoUtil.getFilesMatchingRegexp(base, pattern);
                if (filesMatchingRegexp == null || filesMatchingRegexp.length == 0) {
                    dataFile = null;
                } else if (filesMatchingRegexp.length == 1) {
                    dataFile = filesMatchingRegexp[0];
                } else {
                    throw new PicardException("More than one filter file found in " + base.getAbsolutePath());
                }
            } catch (final IOException e) {
                return Collections.singletonList("Could not create tile index file: " + bci.getAbsolutePath());
            }
            return tileIndex.verify(expectedTiles);
        }

        abstract IlluminaParser<OUTPUT_RECORD> makeParser(List<Integer> requestedTiles);
    }

    class MultiTileFilterFileUtil extends MultiTileFileUtil<PfData> {

        /**
         * @param basecallLaneDir location of .filter file and also .bci file
         */
        MultiTileFilterFileUtil(final File basecallLaneDir) {
            super(".filter", basecallLaneDir, basecallLaneDir, new FilterFileFaker());
        }

        @Override
        IlluminaParser<PfData> makeParser(final List<Integer> requestedTiles) {
            return new MultiTileFilterParser(tileIndex, requestedTiles, dataFile);
        }
    }

    class MultiTileLocsFileUtil extends MultiTileFileUtil<PositionalData> {

        MultiTileLocsFileUtil(final File basecallLaneDir, final File bciDir) {
            super(".locs", basecallLaneDir, bciDir, new MultiTileLocsFileFaker());
        }

        @Override
        IlluminaParser<PositionalData> makeParser(final List<Integer> requestedTiles) {
            return new MultiTileLocsParser(tileIndex, requestedTiles, dataFile, lane);
        }
    }

    /**
     * NextSeq-style bcl's have all tiles for a cycle in a single file.
     */
    class MultiTileBclFileUtil extends ParameterizedFileUtil {
        final File basecallLaneDir;
        final File bci;
        final TileIndex tileIndex;
        final SortedMap<Integer, File> cycleFileMap = new TreeMap<Integer, File>();

        MultiTileBclFileUtil(final File basecallLaneDir) {
            // Since these file names do not contain lane number, first two args to ctor are the same.
            super("^(\\d{4}).bcl.bgzf$", "^(\\d{4}).bcl.bgzf$", ".bcl.bgzf", basecallLaneDir,
                    new MultiTileBclFileFaker());
            this.basecallLaneDir = basecallLaneDir;
            bci = new File(basecallLaneDir, "s_" + lane + ".bci");
            // Do this once rather than when deciding if these files exist and again later.
            final File[] cycleFiles = IoUtil.getFilesMatchingRegexp(base, pattern);
            if (cycleFiles != null) {
                for (final File file : cycleFiles) {
                    final String fileName = file.getName();
                    final String cycleNum = fileName.substring(0, fileName.indexOf('.'));
                    cycleFileMap.put(Integer.valueOf(cycleNum), file);
                }
            }
            if (bci.exists()) {
                tileIndex = new TileIndex(bci);
            } else {
                tileIndex = null;
            }

        }

        public CycleIlluminaFileMap getFiles(final List<Integer> tiles, final int[] cycles) {
            // Filter input list of cycles according to which actually exist
            final ArrayList<Integer> goodCycleList = new ArrayList<Integer>(cycles.length);
            for (final int cycle : cycles) {
                if (cycleFileMap.containsKey(cycle)) {
                    goodCycleList.add(cycle);
                }
            }
            // Ensure cycles are sorted.
            Collections.sort(goodCycleList);
            final int[] goodCycles = new int[goodCycleList.size()];
            for (int i = 0; i < goodCycles.length; ++i) {
                goodCycles[i] = goodCycleList.get(i);
            }

            // Create the map.
            final CycleIlluminaFileMap cycledMap = new CycleIlluminaFileMap();
            if (goodCycles.length > 0) {
                for (final int tile : tiles) {
                    cycledMap.put(tile,
                            new MultiTileBclCycleFilesIterator(basecallLaneDir, lane, tile, goodCycles, extension));
                }
            }
            return cycledMap;
        }

        @Override
        public boolean filesAvailable() {
            return bci.exists() && cycleFileMap.size() > 0;
        }

        @Override
        public LaneTileEnd fileToLaneTileEnd(final String fileName) {
            throw new UnsupportedOperationException();
        }


        @Override
        public List<Integer> getTiles() {
            if (tileIndex == null) {
                return Collections.EMPTY_LIST;
            }
            return tileIndex.getTiles();
        }

        @Override
        public List<String> verify(final List<Integer> expectedTiles, final int[] expectedCycles) {
            if (tileIndex == null) {
                return Collections.singletonList("Tile index(" + bci.getAbsolutePath() + ") does not exist!");
            }
            final List<String> ret = tileIndex.verify(expectedTiles);
            for (final int expectedCycle : expectedCycles) {
                if (!cycleFileMap.containsKey(expectedCycle)) {
                    ret.add(expectedCycle + ".bcl.bgzf not found in " + base);
                }
            }
            return ret;
        }

        @Override
        public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] expectedCycles,
                                      final SupportedIlluminaFormat format) {
            if (tileIndex == null) {
                return Collections.singletonList("Tile index(" + bci.getAbsolutePath() + ") does not exist!");
            }
            final List<String> ret = tileIndex.verify(expectedTiles);
            for (final int expectedCycle : expectedCycles) {
                if (!cycleFileMap.containsKey(expectedCycle)) {
                    ret.add(expectedCycle + ".bcl.bgzf not found in " + base);
                }
            }
            return ret;
        }
    }

    /**
     * A support class for return lane tile and end information for a given file
     */
    static class LaneTileEnd {
        public final Integer lane;
        public final Integer tile;
        public final Integer end;

        public LaneTileEnd(final Integer lane, final Integer tile, final Integer end) {
            this.lane = lane;
            this.tile = tile;
            this.end = end;
        }

        public LaneTileEnd(final Integer lane, final Integer tile) {
            this(lane, tile, null);
        }
    }

    /**
     * Return a regex string for finding Lane and Tile given a file extension pattern
     */
    public static String makeLTRegex(final String fileNameEndPattern) {
        return "^" + UNPARAMETERIZED_PER_TILE_PATTERN + fileNameEndPattern + "$";
    }

    /**
     * Return a regex string for finding Lane and Tile given a file extension pattern
     */
    private static String makeLTRegex(final String fileNameEndPattern, final int lane) {
        return "^" + makeParameterizedLaneAndTileRegex(lane) + fileNameEndPattern + "$";
    }

    private static String makeLaneRegex(final String fileNameEndPattern) {
        return "^s_(\\d+)" + fileNameEndPattern + "$";
    }

    private static String makeLaneRegex(final String fileNameEndPattern, final int lane) {
        return "^s_" + lane + fileNameEndPattern + "$";
    }

    private static int getCycleFromDir(final File tempCycleDir) {
        final char[] name = tempCycleDir.getName().toCharArray();
        if (name[0] != 'C') {
            throw new PicardException("Invalid cycle directory name " + tempCycleDir.getName());
        }

        String intStr = "";
        boolean periodFound = false;
        for (int i = 1; i < name.length && !periodFound; i++) {
            if (name[i] == '.') {
                periodFound = true;
            } else if (name[i] == '1' || name[i] == '2' || name[i] == '3' ||
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

    /**
     * Given a pattern and file name return a LaneTileEnd with the first two matches to the pattern returned
     * as the lane and tile respectively
     */
    private static LaneTileEnd laneAndTileFromFirstTwoMatches(final String fileName, final Pattern pattern) {
        final Matcher matcher = pattern.matcher(fileName);
        if (!matcher.matches()) {
            return null;
        }
        return new LaneTileEnd(Integer.parseInt(matcher.group(1)), Integer.parseInt(matcher.group(2)));
    }

    /**
     * Return a string representing the Lane in the format "L00<lane>"
     *
     * @param lane The lane to transform
     *
     * @return A long string representation of the name
     */
    public static String longLaneStr(final int lane) {
        String lstr = String.valueOf(lane);
        final int zerosToAdd = 3 - lstr.length();

        for (int i = 0; i < zerosToAdd; i++) {
            lstr = "0" + lstr;
        }
        return "L" + lstr;
    }

    /**
     * Return a string representing the Lane in the format "000<tile>"
     *
     * @param tile The tile to transform
     *
     * @return A long string representation of the name
     */
    private static String longTileStr(final int tile) {
        String tstr = String.valueOf(tile);
        final int zerosToAdd = 4 - tstr.length();

        for (int i = 0; i < zerosToAdd; i++) {
            tstr = "0" + tstr;
        }
        return tstr;
    }

    /**
     * Return all files that match pattern of the given file type in the given base directory
     */
    private static IlluminaFileMap getTiledFiles(final File baseDirectory, final Pattern pattern,
                                                 final ParameterizedFileUtil ift) {
        final IlluminaFileMap fileMap = new IlluminaFileMap();
        if (baseDirectory.exists()) {
            IoUtil.assertDirectoryIsReadable(baseDirectory);
            final File[] files = IoUtil.getFilesMatchingRegexp(baseDirectory, pattern);
            for (final File file : files) {
                if ( file.length() > 0) {
                    final LaneTileEnd lt = ift.fileToLaneTileEnd(file.getName());
                    fileMap.put(lt.tile, file);
                }
            }
        }

        return fileMap;
    }

    /**
     * For filename patterns that end with .txt tack on the option .gz extension
     */
    private static String processTxtExtension(final String fileNameEndPattern) {
        if (fileNameEndPattern.endsWith(".txt")) {
            return fileNameEndPattern + "(\\.gz|\\.bz2)?";
        } else {
            return fileNameEndPattern;
        }
    }


    private String liToStr(final List<Integer> intList) {
        if (intList.size() == 0) {
            return "";
        }

        String summary = String.valueOf(intList.get(0));
        for (int i = 1; i < intList.size(); i++) {
            summary += ", " + String.valueOf(intList.get(i));
        }

        return summary;
    }

    private String summarizeTileCounts(final List<SupportedIlluminaFormat> formats) {
        String summary;
        ParameterizedFileUtil pfu = utils.get(formats.get(0));
        List<Integer> tiles = pfu.getTiles();
        summary = pfu.extension + "(" + liToStr(tiles) + ")";

        for (final SupportedIlluminaFormat format : formats) {
            pfu = utils.get(format);
            tiles = pfu.getTiles();

            summary += ", " + pfu.extension + "(" + liToStr(tiles) + ")";
        }

        return summary;
    }

    /**
     * We want to be able to predetermine if the BCL files are gzipped or not and we also want to verify
     * that all of the files are the same. Look through all of the cycle dirs in this lane and grab all
     * BCL (gzipped or not) files in the tree. Determine the exension and then verify that they're all the same.
     * <p/>
     * If there are no BCL files, return the standard extension (i.e. ".bcl") to conserve backwards compatibility
     */
    private String inferBclExtension(final File laneDir) {
        final Pattern bclExtensionPattern = Pattern.compile(".*.bcl(\\.gz)?$");
        final String bclGzipExtension = ".bcl.gz";
        String bclExtension = ".bcl";

        final File[] cycleDirs = IoUtil.getFilesMatchingRegexp(laneDir, CYCLE_SUBDIRECTORY_PATTERN);
        final List<File> allBclFiles = new ArrayList<File>();
        if (cycleDirs != null && cycleDirs.length > 0) {
            // Get all of the BCL files in the various cycle dirs
            for (final File cycleDir : cycleDirs) {
                allBclFiles.addAll(Arrays.asList(IoUtil.getFilesMatchingRegexp(cycleDir, bclExtensionPattern)));
            }

            if (allBclFiles.size() > 0) {
                // Define the extension to be the one the first file has. After that, verify that all files have the
                // same extension
                if (allBclFiles.get(0).getPath().endsWith(bclGzipExtension)) {
                    bclExtension = bclGzipExtension;
                }

                for (final File bclFile : allBclFiles) {
                    if (!bclFile.getPath().endsWith(bclExtension)) {
                        throw new PicardException(
                                "Not all BCL files in " + laneDir.getAbsolutePath() + " have the same extension!");
                    }
                }
            }
        }

        return bclExtension;
    }
}
