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
package picard.illumina.parser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.illumina.parser.fakers.BarcodeFileFaker;
import picard.illumina.parser.fakers.BclFileFaker;
import picard.illumina.parser.fakers.ClocsFileFaker;
import picard.illumina.parser.fakers.FilterFileFaker;
import picard.illumina.parser.fakers.LocsFileFaker;
import picard.illumina.parser.fakers.PosFileFaker;
import picard.illumina.parser.readers.TileMetricsOutReader;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
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
    public static final Pattern CYCLE_SUBDIRECTORY_PATTERN = Pattern.compile("^C(\\d+)\\.1$");

    public enum SupportedIlluminaFormat {
        Bcl,
        Locs,
        Clocs,
        Pos,
        Filter,
        Barcode,
        MultiTileFilter,
        MultiTileLocs,
        MultiTileBcl
    }

    private final File basecallLaneDir;
    private final File intensityLaneDir;
    private final File basecallDir;
    private final File barcodeDir;
    private final File intensityDir;
    private final int lane;

    private final File tileMetricsOut;
    private final Map<SupportedIlluminaFormat, ParameterizedFileUtil> utils = new HashMap<SupportedIlluminaFormat, ParameterizedFileUtil>();

    public IlluminaFileUtil(final File basecallDir, final int lane) {
		this(basecallDir, null, lane);
	}


	public IlluminaFileUtil(final File basecallDir, File barcodeDir, final int lane) {
        this.lane = lane;
        this.basecallDir = basecallDir;
        this.barcodeDir = barcodeDir;
        this.intensityDir = basecallDir.getParentFile();
        final File dataDir = intensityDir.getParentFile();
        this.basecallLaneDir = new File(basecallDir, longLaneStr(lane));
        this.intensityLaneDir = new File(intensityDir, longLaneStr(lane));
        final File interopDir = new File(dataDir.getParentFile(), "InterOp");
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
        ParameterizedFileUtil parameterizedFileUtil = utils.get(format);
        if (parameterizedFileUtil == null) {
            switch (format) {
                case Bcl:
                    final ParameterizedFileUtil bclFileUtil = new PerTilePerCycleFileUtil(".bcl", basecallLaneDir, new BclFileFaker(), lane);
                    final ParameterizedFileUtil gzBclFileUtil = new PerTilePerCycleFileUtil(".bcl.gz", basecallLaneDir, new BclFileFaker(), lane);
                    if (bclFileUtil.filesAvailable() && !gzBclFileUtil.filesAvailable()) {
                        parameterizedFileUtil = bclFileUtil;
                    } else if (!bclFileUtil.filesAvailable() && gzBclFileUtil.filesAvailable()) {
                        parameterizedFileUtil = gzBclFileUtil;
                    } else if (!bclFileUtil.filesAvailable() && !gzBclFileUtil.filesAvailable()) {
                        parameterizedFileUtil = bclFileUtil;
                    } else {
                        throw new PicardException(
                                "Not all BCL files in " + basecallLaneDir.getAbsolutePath() + " have the same extension!");
                    }
                    utils.put(SupportedIlluminaFormat.Bcl, parameterizedFileUtil);
                    break;
                case Locs:
                    parameterizedFileUtil = new PerTileFileUtil(".locs", intensityLaneDir, new LocsFileFaker(), lane);
                    utils.put(SupportedIlluminaFormat.Locs, parameterizedFileUtil);
                    break;
                case Clocs:
                    parameterizedFileUtil = new PerTileFileUtil(".clocs", intensityLaneDir, new ClocsFileFaker(), lane);
                    utils.put(SupportedIlluminaFormat.Clocs, parameterizedFileUtil);
                    break;
                case Pos:
                    parameterizedFileUtil = new PerTileFileUtil("_pos.txt", intensityDir, new PosFileFaker(), lane);
                    utils.put(SupportedIlluminaFormat.Pos, parameterizedFileUtil);
                    break;
                case Filter:
                    parameterizedFileUtil = new PerTileFileUtil(".filter", basecallLaneDir, new FilterFileFaker(), lane);
                    utils.put(SupportedIlluminaFormat.Filter, parameterizedFileUtil);
                    break;
                case Barcode:
                    parameterizedFileUtil = new PerTileFileUtil("_barcode.txt", barcodeDir != null ? barcodeDir : basecallDir, new BarcodeFileFaker(), lane, false);
                    utils.put(SupportedIlluminaFormat.Barcode, parameterizedFileUtil);
                    break;
                case MultiTileFilter:
                    parameterizedFileUtil = new MultiTileFilterFileUtil(basecallLaneDir, lane);
                    utils.put(SupportedIlluminaFormat.MultiTileFilter, parameterizedFileUtil);
                    break;
                case MultiTileLocs:
                    parameterizedFileUtil = new MultiTileLocsFileUtil(new File(intensityDir, basecallLaneDir.getName()), basecallLaneDir, lane);
                    utils.put(SupportedIlluminaFormat.MultiTileLocs, parameterizedFileUtil);
                    break;
                case MultiTileBcl:
                    parameterizedFileUtil = new MultiTileBclFileUtil(basecallLaneDir, lane);
                    utils.put(SupportedIlluminaFormat.MultiTileBcl, parameterizedFileUtil);
                    break;
            }
        }
        return parameterizedFileUtil;
    }

    /**
     * Return the list of tiles we would expect for this lane based on the metrics found in InterOp/TileMetricsOut.bin
     */
    public List<Integer> getExpectedTiles() {
        IOUtil.assertFileIsReadable(tileMetricsOut);
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

        if (formats.isEmpty()) {
            throw new PicardException(
                    "0 Formats were specified.  You need to specify at least SupportedIlluminaFormat to use getTiles");
        }

        final List<Integer> tiles = getUtil(formats.get(0)).getTiles();
        for (int i = 1; i < formats.size(); i++) {
            final List<Integer> fmTiles = getUtil(formats.get(i)).getTiles();
            if (tiles.size() != fmTiles.size() || !tiles.containsAll(fmTiles)) {
                throw new PicardException(
                        "Formats do not have the same number of tiles! " + summarizeTileCounts(formats));
            }
        }

        return tiles;
    }

    public File tileMetricsOut() {
        return tileMetricsOut;
    }

    /*
     * Return a string representing the Lane in the format "L00<lane>"
     *
     * @param lane The lane to transform
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


    private String liToStr(final List<Integer> intList) {
        if (intList.isEmpty()) {
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
        ParameterizedFileUtil pfu = getUtil(formats.get(0));
        List<Integer> tiles = pfu.getTiles();
        summary = pfu.extension + "(" + liToStr(tiles) + ")";

        for (final SupportedIlluminaFormat format : formats) {
            pfu = getUtil(format);
            tiles = pfu.getTiles();

            summary += ", " + pfu.extension + "(" + liToStr(tiles) + ")";
        }

        return summary;
    }
}
