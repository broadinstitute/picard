/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IterableAdapter;
import picard.PicardException;
import picard.illumina.parser.readers.TileMetricsOutReader;
import picard.illumina.parser.readers.TileMetricsOutReader.IlluminaTileMetrics;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Utility for reading the tile data from an Illumina run directory's TileMetricsOut.bin file
 *
 * @author mccowan
 */
public class TileMetricsUtil {
    /** The path to the directory containing the tile metrics file relative to the basecalling directory. */
    public static String INTEROP_SUBDIRECTORY_NAME = "InterOp";
    
    /** The expected name of the tile metrics output file. */
    public static String TILE_METRICS_OUT_FILE_NAME = "TileMetricsOut.bin";

    /** Returns the path to the TileMetrics file given the basecalling directory. */
    public static File renderTileMetricsFileFromBasecallingDirectory(final File illuminaRunDirectory) {
        return new File(new File(illuminaRunDirectory, INTEROP_SUBDIRECTORY_NAME), TILE_METRICS_OUT_FILE_NAME);
    }
    
    /**
     * Returns an unmodifiable collection of tile data read from the provided file. For each tile we will extract:
     *     - lane number
     *     - tile number
     *     - density
     *     - cluster ID
     *     - Phasing & Prephasing for first template read (if available)
     *     - Phasing & Prephasing for second template read (if available)
     */
    public static Collection<Tile> parseTileMetrics(final File tileMetricsOutFile, final ReadStructure readStructure) throws FileNotFoundException {
        // Get the tile metrics lines from TileMetricsOut, keeping only the last value for any Lane/Tile/Code combination
        final Collection<IlluminaTileMetrics> tileMetrics = determineLastValueForLaneTileMetricsCode(new TileMetricsOutReader
                (tileMetricsOutFile));

        // Collect the tiles by lane & tile, and then collect the metrics by lane
        final Map<String, ? extends Collection<IlluminaTileMetrics>> locationToMetricsMap = partitionTileMetricsByLocation(tileMetrics);
        final Collection<Tile> tiles = new LinkedList<>();
        for (final Map.Entry<String, ? extends Collection<IlluminaTileMetrics>> entry : locationToMetricsMap.entrySet()) {
            final Collection<IlluminaTileMetrics> tileRecords = entry.getValue();

            // Get a mapping from metric code number to the corresponding IlluminaTileMetrics
            final Map<Integer, ? extends Collection<IlluminaTileMetrics>> codeMetricsMap = partitionTileMetricsByCode(tileRecords);

            final Set<Integer> observedCodes = codeMetricsMap.keySet();
            if (!(observedCodes.contains(IlluminaMetricsCode.DENSITY_ID.getMetricsCode()) && observedCodes.contains(IlluminaMetricsCode.CLUSTER_ID.getMetricsCode())))
                throw new PicardException(String.format("Expected to find cluster and density record codes (%s and %s) in records read for tile location %s (lane:tile), but found only %s.",
                        IlluminaMetricsCode.CLUSTER_ID.getMetricsCode(), IlluminaMetricsCode.DENSITY_ID.getMetricsCode(), entry.getKey(), observedCodes));

            final IlluminaTileMetrics densityRecord = CollectionUtil.getSoleElement(codeMetricsMap.get(IlluminaMetricsCode.DENSITY_ID.getMetricsCode()));
            final IlluminaTileMetrics clusterRecord = CollectionUtil.getSoleElement(codeMetricsMap.get(IlluminaMetricsCode.CLUSTER_ID.getMetricsCode()));

            // Snag the phasing data for each read in the read structure. For both types of phasing values, this is the median of all of the individual values seen
            final Collection<TilePhasingValue> tilePhasingValues = getTilePhasingValues(codeMetricsMap, readStructure);

            tiles.add(new Tile(densityRecord.getLaneNumber(), densityRecord.getTileNumber(), densityRecord.getMetricValue(), clusterRecord.getMetricValue(),
                tilePhasingValues.toArray(new TilePhasingValue[tilePhasingValues.size()])));
        }

        return Collections.unmodifiableCollection(tiles);
    }

    /** Pulls out the phasing & prephasing value for the template reads and returns a collection of TilePhasingValues representing these */
    private static Collection<TilePhasingValue> getTilePhasingValues(final Map<Integer, ? extends Collection<IlluminaTileMetrics>> codeMetricsMap, final ReadStructure readStructure) {
        boolean isFirstRead = true;
        final Collection<TilePhasingValue> tilePhasingValues = new ArrayList<>();
        for (int descriptorIndex = 0; descriptorIndex < readStructure.descriptors.size(); descriptorIndex++) {
            if (readStructure.descriptors.get(descriptorIndex).type == ReadType.Template) {
                final TileTemplateRead tileTemplateRead = isFirstRead ? TileTemplateRead.FIRST : TileTemplateRead.SECOND;
                // For both phasing & prephasing, pull out the value and create a TilePhasingValue for further processing
                final int phasingCode = IlluminaMetricsCode.getPhasingCode(descriptorIndex, IlluminaMetricsCode.PHASING_BASE);
                final int prePhasingCode = IlluminaMetricsCode.getPhasingCode(descriptorIndex, IlluminaMetricsCode.PREPHASING_BASE);

                if (!(codeMetricsMap.containsKey(phasingCode) && codeMetricsMap.containsKey(prePhasingCode))) {
                    throw new PicardException("Don't have both phasing and prephasing values for tile");
                }

                tilePhasingValues.add(new TilePhasingValue(tileTemplateRead,
                        CollectionUtil.getSoleElement(codeMetricsMap.get(phasingCode)).getMetricValue(),
                        CollectionUtil.getSoleElement(codeMetricsMap.get(prePhasingCode)).getMetricValue()));
                isFirstRead = false;
            }
        }

        return tilePhasingValues;
    }

    /** According to Illumina, for every lane/tile/code combination they will only use the last value. Filter out the previous values */
    private static Collection<IlluminaTileMetrics> determineLastValueForLaneTileMetricsCode(final Iterator<IlluminaTileMetrics>
                                                                                                    tileMetricsIterator) {
        final Map<TileMetricsOutReader.IlluminaLaneTileCode, IlluminaTileMetrics> filteredTileMetrics = new HashMap<>();
        for (final IlluminaTileMetrics illuminaTileMetrics : new IterableAdapter<>(tileMetricsIterator)) {
            filteredTileMetrics.put(illuminaTileMetrics.getLaneTileCode(), illuminaTileMetrics);
        }

        return filteredTileMetrics.values();
    }

    private static String renderMetricLocationKey(final IlluminaTileMetrics metric) {
        return String.format("%s:%s", metric.getLaneNumber(), metric.getTileNumber());
    }

    // Wrapper around CollectionUtil.Partitioner, purely to de-bulk the actual methods
    private static Map<Integer, ? extends Collection<IlluminaTileMetrics>> partitionTileMetricsByCode(final Collection<IlluminaTileMetrics> tileMetrics) {
        return tileMetrics.stream().collect(Collectors.groupingBy(IlluminaTileMetrics::getMetricCode));
    }

    // Wrapper around CollectionUtil.Partitioner, purely to de-bulk the actual methods
    private static Map<String, ? extends Collection<IlluminaTileMetrics>> partitionTileMetricsByLocation(final Collection<IlluminaTileMetrics> tileMetrics) {
        return tileMetrics.stream().collect(Collectors.groupingBy(TileMetricsUtil::renderMetricLocationKey));
    }

}
