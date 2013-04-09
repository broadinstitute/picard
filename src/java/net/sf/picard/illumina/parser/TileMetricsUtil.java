package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;
import net.sf.picard.illumina.parser.readers.TileMetricsOutReader;
import net.sf.picard.illumina.parser.readers.TileMetricsOutReader.IlluminaTileMetrics;
import net.sf.picard.util.CollectionUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Utility for reading the tile data from an Illumina run directory's TileMetricsOut.bin file
 *
 * @author mccowan
 */
public class TileMetricsUtil {
    private final static Integer DENSITY_ID_CODE = 100;
    private final static Integer CLUSTER_ID_CODE = 102;

    /** The path to the directory containing the tile metrics file relative to the basecalling directory. */
    public static String INTEROP_SUBDIRECTORY_NAME = "InterOp";
    
    /** The expected name of the tile metrics output file. */
    public static String TILE_METRICS_OUT_FILE_NAME = "TileMetricsOut.bin";

    /** Returns the path to the TileMetrics file given the basecalling directory. */
    public static File renderTileMetricsFileFromBasecallingDirectory(final File illuminaRunDirectory) {
        return new File(new File(illuminaRunDirectory, INTEROP_SUBDIRECTORY_NAME), TILE_METRICS_OUT_FILE_NAME);
    }
    
    /**
     * Returns an unmodifiable collection of tile data read from the provided file.
     */
    public static Collection<Tile> parseTileMetrics(final File tileMetricsOutFile) throws FileNotFoundException {

        final Collection<IlluminaTileMetrics> metrics = CollectionUtil.makeCollection(new TileMetricsOutReader(tileMetricsOutFile));
        final Map<String, Collection<IlluminaTileMetrics>> locationToMetricsMap = CollectionUtil.partition(metrics, new CollectionUtil.Partitioner<IlluminaTileMetrics, String>() {
            @Override
            public String getPartition(final IlluminaTileMetrics metric) {
                return renderMetricLocationKey(metric);
            }
        });

        final Collection<Tile> tiles = new LinkedList<Tile>();
        for (final Map.Entry<String, Collection<IlluminaTileMetrics>> entry : locationToMetricsMap.entrySet()) {
            final Collection<IlluminaTileMetrics> tileRecords = entry.getValue();
            final Map<Integer, Collection<IlluminaTileMetrics>> codeMetricsMap = CollectionUtil.partition(tileRecords, new CollectionUtil.Partitioner<IlluminaTileMetrics, Integer>() {
                @Override
                public Integer getPartition(final IlluminaTileMetrics metric) {
                    return metric.getMetricCode();
                }
            });
            final Set<Integer> observedCodes = codeMetricsMap.keySet();
            if (!(observedCodes.contains(DENSITY_ID_CODE) && observedCodes.contains(CLUSTER_ID_CODE)))
                throw new PicardException(String.format("Expected to find cluster and density record codes (%s and %s) in records read for tile location %s (lane:tile), but found only %s.", CLUSTER_ID_CODE, DENSITY_ID_CODE, entry.getKey(), observedCodes));

            final IlluminaTileMetrics densityRecord = CollectionUtil.getSoleElement(codeMetricsMap.get(DENSITY_ID_CODE));
            final IlluminaTileMetrics clusterRecord = CollectionUtil.getSoleElement(codeMetricsMap.get(CLUSTER_ID_CODE));
            tiles.add(new Tile(densityRecord.getLaneNumber(), densityRecord.getTileNumber(), densityRecord.getMetricValue(), clusterRecord.getMetricValue()));
        }

        return Collections.unmodifiableCollection(tiles);
    }

    
    
    private static String renderMetricLocationKey(final IlluminaTileMetrics metric) {
        return String.format("%s:%s", metric.getLaneNumber(), metric.getTileNumber());
    }

    /**
     * Describes a tile.
     */
    public static class Tile {
        private final int lane, tile;
        private final float density, clusters;

        /** Returns the number of this tile's parent lane. */
        public int getLaneNumber() {
            return lane;
        }

        protected Tile(final int lane, final int tile, final float density, final float clusters) {
            this.lane = lane;
            this.tile = tile;
            this.density = density;
            this.clusters = clusters;
        }

        /** Returns the number/name of this tile. */
        public int getTileNumber() {
            return tile;
        }

        /** Returns the cluster density of this tile, in units of [cluster/mm^2]. */
        public float getClusterDensity() {
            return density;
        }

        /** Returns the number of on this tile. */
        public float getClusterCount() {
            return clusters;
        }
    }
}
