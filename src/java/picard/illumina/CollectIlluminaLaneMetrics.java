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

package picard.illumina;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Illumina;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.Tile;
import picard.illumina.parser.TileMetricsUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Collection;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Command-line wrapper around {@link IlluminaLaneMetricsCollector}.
 * @author mccowan
 */

@CommandLineProgramProperties(
        usage = CollectIlluminaLaneMetrics.USAGE,
        usageShort = CollectIlluminaLaneMetrics.USAGE,
        programGroup = Illumina.class
)
public class CollectIlluminaLaneMetrics extends CommandLineProgram {
    static final String USAGE = "Collects Illumina lane metrics for the given basecalling analysis directory";

    @Option(doc = "The Illumina run directory of the run for which the lane metrics are to be generated")
    public File RUN_DIRECTORY;

    @Option(doc = "The directory to which the output file will be written")
    public File OUTPUT_DIRECTORY;

    @Option(doc = "The prefix to be prepended to the file name of the output file; an appropriate suffix will be applied", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public String OUTPUT_PREFIX;

    @Option(doc = ReadStructure.PARAMETER_DOC, shortName = "RS")
    public ReadStructure READ_STRUCTURE;

    @Override
    protected int doWork() {
        final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile = this.getMetricsFile();
        final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile = this.getMetricsFile();
        IlluminaLaneMetricsCollector.collectLaneMetrics(RUN_DIRECTORY, OUTPUT_DIRECTORY, OUTPUT_PREFIX, laneMetricsFile, phasingMetricsFile, READ_STRUCTURE);
        return 0;
    }

    public static void main(final String[] args) {
        new CollectIlluminaLaneMetrics().instanceMainWithExit(args);
    }

    /**
     * Utility for collating Tile records from the Illumina TileMetrics file into lane-level and phasing-level metrics.
     */
    public static class IlluminaLaneMetricsCollector {

        private final static Log LOG = Log.getInstance(IlluminaLaneMetricsCollector.class);

        /** Returns a partitioned collection of lane number to Tile objects from the provided basecall directory. */
        public static Map<Integer, ? extends Collection<Tile>> readLaneTiles(final File illuminaRunDirectory, final ReadStructure readStructure) {
            final Collection<Tile> tiles;
            try {
                tiles = TileMetricsUtil.parseTileMetrics(TileMetricsUtil.renderTileMetricsFileFromBasecallingDirectory(illuminaRunDirectory), readStructure);
            } catch (final FileNotFoundException e) {
                throw new PicardException("Unable to open laneMetrics file.", e);
            }

            return tiles.stream().collect(Collectors.groupingBy(Tile::getLaneNumber));
        }

        /** Parses the tile data from the basecall directory and writes to both the lane and phasing metrics files */
        public static void collectLaneMetrics(final File runDirectory, final File outputDirectory, final String outputPrefix,
                                              final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile,
                                              final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile,
                                              final ReadStructure readStructure) {
            final Map<Integer, ? extends Collection<Tile>> laneTiles = readLaneTiles(runDirectory, readStructure);
            writeLaneMetrics(laneTiles, outputDirectory, outputPrefix, laneMetricsFile);
            writePhasingMetrics(laneTiles, outputDirectory, outputPrefix, phasingMetricsFile);
        }

        public static File writePhasingMetrics(final Map<Integer, ? extends Collection<Tile>> laneTiles, final File outputDirectory,
                                               final String outputPrefix, final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile) {
            laneTiles.entrySet().stream().forEach(entry -> IlluminaPhasingMetrics.getPhasingMetricsForTiles(entry.getKey().longValue(),
                    entry.getValue()).forEach(phasingMetricsFile::addMetric));

            return writeMetrics(phasingMetricsFile, outputDirectory, outputPrefix, IlluminaPhasingMetrics.getExtension());
        }

        public static File writeLaneMetrics(final Map<Integer, ? extends Collection<Tile>> laneTiles, final File outputDirectory,
                                            final String outputPrefix, final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile) {
            laneTiles.entrySet().stream().forEach(entry -> {
                final IlluminaLaneMetrics laneMetric = new IlluminaLaneMetrics();
                laneMetric.LANE = entry.getKey().longValue();
                laneMetric.CLUSTER_DENSITY = calculateLaneDensityFromTiles(entry.getValue());
                laneMetricsFile.addMetric(laneMetric);
            });

            return writeMetrics(laneMetricsFile, outputDirectory, outputPrefix, IlluminaLaneMetrics.getExtension());
        }

        private static File writeMetrics(final MetricsFile<MetricBase, Comparable<?>> metricsFile, final File outputDirectory,
                                         final String outputPrefix, final String outputExtension) {
            final File outputFile = new File(outputDirectory, String.format("%s.%s", outputPrefix, outputExtension));
            LOG.info(String.format("Writing %s lane metrics to %s ...", metricsFile.getMetrics().size(), outputFile));
            metricsFile.write(outputFile);
            return outputFile;
        }

        private static double calculateLaneDensityFromTiles(final Collection<Tile> tiles) {
            double area = 0;
            double clusters = 0;
            for (final Tile tile : tiles) {
                area += (tile.getClusterCount() / tile.getClusterDensity());
                clusters += tile.getClusterCount();
            }
            return clusters / area;
        }
    }
}
