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

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.BaseCallingProgramGroup;
import picard.illumina.parser.*;
import picard.illumina.parser.readers.TileMetricsOutReader;

import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Command-line wrapper around {@link IlluminaLaneMetricsCollector}.
 * @author mccowan
 */

@CommandLineProgramProperties(
        summary = CollectIlluminaLaneMetrics.USAGE_SUMMARY + CollectIlluminaLaneMetrics.USAGE_DETAILS,
        oneLineSummary = CollectIlluminaLaneMetrics.USAGE_SUMMARY,
        programGroup = BaseCallingProgramGroup.class
)
@DocumentedFeature
public class CollectIlluminaLaneMetrics extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Collects Illumina lane metrics for the given BaseCalling analysis directory.";
    static final String USAGE_DETAILS = "This tool produces quality control metrics on cluster density for each lane of an Illumina flowcell. " +
            "This tool takes Illumina TileMetrics data and places them into directories containing lane- and phasing-level metrics. " +
            "In this context, phasing refers to the fraction of molecules that fall behind or jump ahead (prephasing) during a read cycle." +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectIlluminaLaneMetrics \\<br />" +
            "      RUN_DIR=test_run \\<br />" +
            "      OUTPUT_DIRECTORY=Lane_output_metrics \\<br />" +
            "      OUTPUT_PREFIX=experiment1 \\<br />" +
            "      READ_STRUCTURE=25T8B25T " +
            "</pre>" +
            "<p>Please see the CollectIlluminaLaneMetrics " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectIlluminaLaneMetrics'>definitions</a> " +
            "for a complete description of the metrics produced by this tool.</p>" +
            "<hr />"
    ;
    @Argument(doc = "The Illumina run directory of the run for which the lane metrics are to be generated")
    public File RUN_DIRECTORY;

    @Argument(doc = "The directory to which the output file will be written")
    public File OUTPUT_DIRECTORY;

    @Argument(doc = "The prefix to be prepended to the file name of the output file; an appropriate suffix will be applied", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public String OUTPUT_PREFIX;

    @Argument(doc = ReadStructure.PARAMETER_DOC + "\nIf not given, will use the RunInfo.xml in the run directory.", shortName = "RS", optional = true)
    public ReadStructure READ_STRUCTURE;

    @Argument(doc="Append the given file extension to all metric file names (ex. OUTPUT.illumina_lane_metrics.EXT). None if null", shortName = "EXT", optional = true)
    public String FILE_EXTENSION = null;

    @Override
    protected int doWork() {
        final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile = this.getMetricsFile();
        final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile = this.getMetricsFile();

        if (READ_STRUCTURE == null) {
            final File runInfo = new File(RUN_DIRECTORY + "/" + "RunInfo.xml");
            IOUtil.assertFileIsReadable(runInfo);
            try {
                final Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(runInfo);
                final NodeList reads = document.getElementsByTagName("Read");
                final List<ReadDescriptor> descriptors = new ArrayList<>(reads.getLength());
                for (int i = 0; i < reads.getLength(); i++) {
                    final Node read = reads.item(i);
                    final NamedNodeMap attributes = read.getAttributes();
                    final int readNumber = Integer.parseInt(attributes.getNamedItem("Number").getNodeValue());
                    final int numCycles = Integer.parseInt(attributes.getNamedItem("NumCycles").getNodeValue());
                    final boolean isIndexedRead = attributes.getNamedItem("IsIndexedRead").getNodeValue().toUpperCase().equals("Y");
                    if (readNumber != i + 1) throw new PicardException("Read number in RunInfo.xml was out of order: " + (i+1) + " != " + readNumber);
                    descriptors.add(new ReadDescriptor(numCycles, isIndexedRead ? ReadType.Barcode: ReadType.Template));
                }
                READ_STRUCTURE = new ReadStructure(descriptors);
            } catch (final Exception e) {
                throw new PicardException(e.getMessage());
            }
        }

        IlluminaLaneMetricsCollector.collectLaneMetrics(RUN_DIRECTORY, OUTPUT_DIRECTORY, OUTPUT_PREFIX,
                laneMetricsFile, phasingMetricsFile,
                READ_STRUCTURE, FILE_EXTENSION == null ? "" : FILE_EXTENSION, VALIDATION_STRINGENCY);
        return 0;
    }

    /**
     * Utility for collating Tile records from the Illumina TileMetrics file into lane-level and phasing-level metrics.
     */
    public static class IlluminaLaneMetricsCollector {

        private final static Log LOG = Log.getInstance(IlluminaLaneMetricsCollector.class);

        /**
         * Returns a partitioned collection of lane number to Tile objects from the provided basecall directory.
         */
        public static Map<Integer, ? extends Collection<Tile>> readLaneTiles(final File illuminaRunDirectory,
                                                                             final ReadStructure readStructure,
                                                                             final ValidationStringency validationStringency,
                                                                             final int tileMetricsVersion) {
            final Collection<Tile> tiles;

            final List<File> tileMetricsOutFiles = TileMetricsUtil.findTileMetricsFiles(illuminaRunDirectory, readStructure.totalCycles);
            if (tileMetricsVersion == TileMetricsOutReader.TileMetricsVersion.THREE.version) {
                tiles = TileMetricsUtil.parseClusterRecordsFromTileMetrics(
                        tileMetricsOutFiles,
                        TileMetricsUtil.renderPhasingMetricsFilesFromBasecallingDirectory(illuminaRunDirectory),
                        readStructure
                );
            } else {
                tiles = TileMetricsUtil.parseTileMetrics(
                        tileMetricsOutFiles.get(0),
                        readStructure,
                        validationStringency
                );
            }
            return tiles.stream().filter(tile -> tile.getLaneNumber() > 0).collect(Collectors.groupingBy(Tile::getLaneNumber));
        }

        /** Parses the tile data from the basecall directory and writes to both the lane and phasing metrics files */
        public static void collectLaneMetrics(final File runDirectory, final File outputDirectory, final String outputPrefix,
                                              final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile,
                                              final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile,
                                              final ReadStructure readStructure, final String fileExtension,
                                              final ValidationStringency validationStringency) {
            int tileMetricsVersion = determineTileMetricsVersion(runDirectory, readStructure);
            final Map<Integer, ? extends Collection<Tile>> laneTiles = readLaneTiles(runDirectory, readStructure, validationStringency, tileMetricsVersion);
            writeLaneMetrics(laneTiles, outputDirectory, outputPrefix, laneMetricsFile, fileExtension);
            writePhasingMetrics(laneTiles, outputDirectory, outputPrefix, phasingMetricsFile, fileExtension, tileMetricsVersion);
        }

        private static int determineTileMetricsVersion(File illuminaRunDirectory, ReadStructure readStructure) {
            final List<File> tileMetricsOutFiles = TileMetricsUtil.findTileMetricsFiles(illuminaRunDirectory, readStructure.totalCycles);
            int version = new TileMetricsOutReader(tileMetricsOutFiles.get(0)).getVersion();
            if (!tileMetricsOutFiles.stream().allMatch(metricFile -> {
                int fileVersion = new TileMetricsOutReader(metricFile).getVersion();
                boolean matches = fileVersion == version;
                if (!matches) {
                    LOG.error(String.format("File %s version %d does not match expected version %d.", metricFile.getAbsolutePath(), fileVersion, version));
                }
                return matches;
            })) {
                throw new PicardException("Not all tile metrics files match expected version: " + version);
            }
            return version;
        }

        private static void writePhasingMetrics(final Map<Integer, ? extends Collection<Tile>> laneTiles, final File outputDirectory,
                                               final String outputPrefix, final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile,
                                               final String fileExtension, final int tileMetricsVersion) {
            laneTiles.forEach((key, value) -> IlluminaPhasingMetrics.getPhasingMetricsForTiles(key.longValue(),
                    value, tileMetricsVersion == TileMetricsOutReader.TileMetricsVersion.TWO.version).forEach(phasingMetricsFile::addMetric));

            writeMetrics(phasingMetricsFile, outputDirectory, outputPrefix, IlluminaPhasingMetrics.getExtension() + fileExtension);
        }

        private static void writeLaneMetrics(final Map<Integer, ? extends Collection<Tile>> laneTiles, final File outputDirectory,
                                            final String outputPrefix, final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile,
                                            final String fileExtension) {
            laneTiles.forEach((key, value) -> {
                final IlluminaLaneMetrics laneMetric = new IlluminaLaneMetrics();
                laneMetric.LANE = key.longValue();
                laneMetric.CLUSTER_DENSITY = calculateLaneDensityFromTiles(value);
                laneMetricsFile.addMetric(laneMetric);
            });

            writeMetrics(laneMetricsFile, outputDirectory, outputPrefix, IlluminaLaneMetrics.getExtension() + fileExtension);
        }

        private static void writeMetrics(final MetricsFile<MetricBase, Comparable<?>> metricsFile, final File outputDirectory,
                                         final String outputPrefix, final String outputExtension) {
            final File outputFile = new File(outputDirectory, String.format("%s.%s", outputPrefix, outputExtension));
            LOG.info(String.format("Writing %s lane metrics to %s ...", metricsFile.getMetrics().size(), outputFile));
            metricsFile.write(outputFile);
        }

        private static double calculateLaneDensityFromTiles(final Collection<Tile> tiles) {
            double area = 0;
            double clusters = 0;
            for (final Tile tile : tiles) {
                if (tile.getClusterDensity() > 0) area += (tile.getClusterCount() / tile.getClusterDensity());
                clusters += tile.getClusterCount();
            }
            return (area > 0) ? clusters / area : 0.0;
        }
    }
}
