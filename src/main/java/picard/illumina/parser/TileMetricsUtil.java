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

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.illumina.parser.readers.EmpiricalPhasingMetricsOutReader;
import picard.illumina.parser.readers.TileMetricsOutReader;
import picard.illumina.parser.readers.TileMetricsOutReader.IlluminaTileMetrics;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Utility for reading the tile data from an Illumina run directory's TileMetricsOut.bin file
 *
 * @author mccowan
 */
public class TileMetricsUtil {

    private static final Log log = Log.getInstance(TileMetricsUtil.class);

    /**
     * The path to the directory containing the tile metrics file relative to the basecalling directory.
     */
    public static String INTEROP_SUBDIRECTORY_NAME = "InterOp";

    /**
     * The expected name of the tile metrics output file.
     */
    public static String TILE_METRICS_OUT_FILE_NAME = "TileMetricsOut.bin";

    private final static Log LOG = Log.getInstance(TileMetricsUtil.class);

    /**
     * Finds all of the tile metrics files for a given run directory and cycle count.
     * @param illuminaRunDirectory The run directory
     * @param numCycles The number of cycle directories to check.
     * @return A list of all tile metrics files.
     */
    public static List<File> findTileMetricsFiles(File illuminaRunDirectory, int numCycles) {
        Path interOpDir = illuminaRunDirectory.toPath().resolve(INTEROP_SUBDIRECTORY_NAME);
        final List<Path> pathsToTest = new ArrayList<>();

        pathsToTest.add(interOpDir.resolve(TILE_METRICS_OUT_FILE_NAME));

        // check cycles in reverse order.
        for (int i = numCycles; i > 0; i--) {
            pathsToTest.add(interOpDir.resolve(String.format("C%d.1/%s", i, TILE_METRICS_OUT_FILE_NAME)));
        }

        final List<File> files = pathsToTest.stream()
                .filter(Files::exists)
                .map(Path::toFile)
                .collect(Collectors.toList());
        if (files.isEmpty()) {
            throw new IllegalStateException(String.format("No %s file found in %s or any of its cycle directories.\"", INTEROP_SUBDIRECTORY_NAME, interOpDir));
        }
        return files;
    }

    private static Collection<Tile> getTileClusterRecords(
            final Map<String, ? extends Collection<IlluminaTileMetrics>> locationToMetricsMap,
            final Map<Integer, Map<Integer, Collection<TilePhasingValue>>> phasingValues,
            final float density) {

        final Collection<Tile> tiles = new LinkedList<>();
        for (final Map.Entry<String, ? extends Collection<IlluminaTileMetrics>> entry : locationToMetricsMap.entrySet()) {
            final Collection<IlluminaTileMetrics> tileRecords = entry.getValue();

            final IlluminaTileMetrics record = CollectionUtil.getSoleElement(tileRecords);

            //only create for cluster records
            if (record.isClusterRecord()) {
                final Collection<TilePhasingValue> tilePhasingValues = phasingValues.get(record.getLaneNumber()).get(record.getTileNumber());
                tiles.add(new Tile(record.getLaneNumber(), record.getTileNumber(), density, record.getMetricValue(),
                        tilePhasingValues.toArray(new TilePhasingValue[tilePhasingValues.size()])));
            }
        }
        return Collections.unmodifiableCollection(tiles);
    }

    public static Collection<Tile> parseClusterRecordsFromTileMetrics(
            final Collection<File> tileMetricsOutFiles,
            final Map<Integer, File> phasingMetricsFiles,
            final ReadStructure readStructure
    ) {
        final Map<Integer, Map<Integer, Collection<TilePhasingValue>>> phasingValues = getTilePhasingValues(phasingMetricsFiles, readStructure);
        for (File tileMetricsOutFile : tileMetricsOutFiles) {
            final TileMetricsOutReader tileMetricsIterator = new TileMetricsOutReader(tileMetricsOutFile);
            final float density = tileMetricsIterator.getDensity();
            final Collection<IlluminaTileMetrics> tileMetrics = determineLastValueForLaneTileMetricsCode(tileMetricsIterator);
            final Map<String, ? extends Collection<IlluminaTileMetrics>> locationToMetricsMap = partitionTileMetricsByLocation(tileMetrics);
            final Collection<Tile> tiles = getTileClusterRecords(locationToMetricsMap, phasingValues, density);
            if (!tiles.isEmpty()) {
                return tiles;
            }
        }
        final String pathsString = tileMetricsOutFiles
                .stream()
                .map(File::getAbsolutePath)
                .collect(Collectors.joining(", "));
        throw new RuntimeException("None of the following input files contained cluster records: " + pathsString);
    }

    /**
     * Returns an unmodifiable collection of tile data read from the provided file. For each tile we will extract:
     * - lane number
     * - tile number
     * - density
     * - cluster ID
     * - Phasing & Prephasing for first template read (if available)
     * - Phasing & Prephasing for second template read (if available)
     */
    public static Collection<Tile> parseTileMetrics(final File tileMetricsOutFile, final ReadStructure readStructure,
                                                    final ValidationStringency validationStringency) {
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
            final Collection<TilePhasingValue> tilePhasingValues = getTilePhasingValues(codeMetricsMap, readStructure, validationStringency);

            tiles.add(new Tile(densityRecord.getLaneNumber(), densityRecord.getTileNumber(), densityRecord.getMetricValue(), clusterRecord.getMetricValue(),
                    tilePhasingValues.toArray(new TilePhasingValue[tilePhasingValues.size()])));
        }

        return Collections.unmodifiableCollection(tiles);
    }

    /**
     * Parse phasing metrics files for phasing data for each read in the read structure for each tile.
     * For both types of phasing values, this is the median of all of the individual values seen
     */
    private static Map<Integer, Map<Integer, Collection<TilePhasingValue>>> getTilePhasingValues(
                                                                     Map<Integer, File> phasingMetricFiles,
                                                                     final ReadStructure readStructure) {
        final Map<Integer, Map<Integer, Collection<TilePhasingValue>>> phasingValues = new HashMap<>();
        int totalCycleCount = 0;

        boolean isFirstRead = true;

        int readNum = 0;
        for (int outputLength : readStructure.readLengths) {
            ReadDescriptor descriptor = readStructure.descriptors.get(readNum++);
            if (descriptor.type == ReadType.Template) {
                final TileTemplateRead tileTemplateRead = isFirstRead ? TileTemplateRead.FIRST : TileTemplateRead.SECOND;

                List<Float> cycleNumWithData = new ArrayList<>();
                Map<Integer, Map<Integer, List<Float>>> phasing = new HashMap<>();
                Map<Integer, Map<Integer, List<Float>>> prePhasing = new HashMap<>();
                for (int cycle = 0; cycle < outputLength; cycle++) {
                    File phasingData = phasingMetricFiles.get(totalCycleCount + 1);
                    if (phasingData != null) {
                        cycleNumWithData.add((float) (cycle + 1));
                        EmpiricalPhasingMetricsOutReader reader = new EmpiricalPhasingMetricsOutReader(phasingData);

                        while (reader.hasNext()) {
                            EmpiricalPhasingMetricsOutReader.IlluminaPhasingMetrics phasingMetrics = reader.next();

                            TileMetricsOutReader.IlluminaLaneTileCode laneTileCode = phasingMetrics.laneTileCode;
                            int tileNumber = laneTileCode.getTileNumber();
                            int laneNumber = laneTileCode.getLaneNumber();
                            phasing
                                .computeIfAbsent(tileNumber, k -> new HashMap<>())
                                .computeIfAbsent(laneNumber, k -> new ArrayList<>())
                                .add(phasingMetrics.phasingWeight);
                            prePhasing
                                .computeIfAbsent(tileNumber, k -> new HashMap<>())
                                .computeIfAbsent(laneNumber, k -> new ArrayList<>())
                                .add(phasingMetrics.prephasingWeight);
                        }
                    }

                    totalCycleCount++;
                }

                phasing.forEach((tileNum, lanePhasingMetrics) -> {
                    Map<Integer, List<Float>> lanePrePhasingMetrics = prePhasing.get(tileNum);
                    lanePhasingMetrics.forEach((laneNum, phasingMetrics) -> {
                        List<Float> prephasingMetrics = lanePrePhasingMetrics.get(laneNum);
                        float[] phasingSlopeAndOffset = computeLinearFit(cycleNumWithData.toArray(new Float[0]), phasingMetrics.toArray(new Float[0]), phasingMetrics.size());
                        float[] prePhasingSlopeAndOffset = computeLinearFit(cycleNumWithData.toArray(new Float[0]), prephasingMetrics.toArray(new Float[0]), phasingMetrics.size());
                        TilePhasingValue value = new TilePhasingValue(tileTemplateRead, phasingSlopeAndOffset[0], prePhasingSlopeAndOffset[0]);
                        phasingValues
                            .computeIfAbsent(laneNum, k -> new HashMap<>())
                            .computeIfAbsent(tileNum, k -> new ArrayList<>())
                            .add(value);
                    });
                });

                isFirstRead = false;
            } else {
                totalCycleCount += outputLength;
            }
        }

        return phasingValues;
    }

    /**
     * Pulls out the phasing & prephasing value for the template reads and returns a collection of TilePhasingValues representing these
     */
    private static Collection<TilePhasingValue> getTilePhasingValues(final Map<Integer, ? extends Collection<IlluminaTileMetrics>> codeMetricsMap, final ReadStructure readStructure, final ValidationStringency validationStringency) {
        boolean isFirstRead = true;
        final Collection<TilePhasingValue> tilePhasingValues = new ArrayList<>();
        for (int descriptorIndex = 0; descriptorIndex < readStructure.descriptors.size(); descriptorIndex++) {
            if (readStructure.descriptors.get(descriptorIndex).type == ReadType.Template) {
                final TileTemplateRead tileTemplateRead = isFirstRead ? TileTemplateRead.FIRST : TileTemplateRead.SECOND;
                // For both phasing & prephasing, pull out the value and create a TilePhasingValue for further processing
                final int phasingCode = IlluminaMetricsCode.getPhasingCode(descriptorIndex, IlluminaMetricsCode.PHASING_BASE);
                final int prePhasingCode = IlluminaMetricsCode.getPhasingCode(descriptorIndex, IlluminaMetricsCode.PREPHASING_BASE);

                final float phasingValue, prePhasingValue;

                // If both the phasing and pre-phasing data are missing, then likely something went wrong when imaging
                // this tile, for example a grain of sand disrupting the path of light to the sensor.  If only one of them
                // is missing, then likely the data is corrupt.
                if (codeMetricsMap.containsKey(phasingCode) && codeMetricsMap.containsKey(prePhasingCode)) {
                    phasingValue = CollectionUtil.getSoleElement(codeMetricsMap.get(phasingCode)).getMetricValue();
                    prePhasingValue = CollectionUtil.getSoleElement(codeMetricsMap.get(prePhasingCode)).getMetricValue();
                } else {
                    final String message = String.format(
                            "Don't have both phasing and prephasing values for %s read cycle %s.  Phasing code was %d and prephasing code was %d.",
                            tileTemplateRead.toString(), descriptorIndex + 1, phasingCode, prePhasingCode
                    );
                    if (!codeMetricsMap.containsKey(phasingCode) && !codeMetricsMap.containsKey(prePhasingCode) && validationStringency != ValidationStringency.STRICT) {
                        // Ignore the error, and use the default (zero) for the phasing values
                        if (validationStringency == ValidationStringency.LENIENT) {
                            LOG.warn(message);
                        }
                    } else {
                        throw new PicardException(message);
                    }
                    phasingValue = 0;
                    prePhasingValue = 0;
                }

                tilePhasingValues.add(new TilePhasingValue(tileTemplateRead, phasingValue, prePhasingValue));
                isFirstRead = false;
            }
        }

        return tilePhasingValues;
    }

    /**
     * According to Illumina, for every lane/tile/code combination they will only use the last value. Filter out the previous values
     */
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

    public static Map<Integer, File> renderPhasingMetricsFilesFromBasecallingDirectory(File illuminaRunDirectory) {
        File[] cycleDirs = IOUtil.getFilesMatchingRegexp(new File(illuminaRunDirectory, INTEROP_SUBDIRECTORY_NAME),
                IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN);

        Map<Integer, File> phasingMetrics = new HashMap<>();
        Arrays.asList(cycleDirs)
                .forEach(cycleDir -> {
                    File[] filesMatchingRegexp = IOUtil.getFilesMatchingRegexp(
                            cycleDir, "EmpiricalPhasingMetricsOut.bin");
                    if (filesMatchingRegexp.length > 0) {
                        phasingMetrics.put(PerTilePerCycleFileUtil.getCycleFromDir(cycleDir),
                                filesMatchingRegexp[0]);
                    }
                });
        return phasingMetrics;
    }

    /*line fit for phasing weights
    /** Compute the best fit line (slope + offset) given the points defined by (x_values[i], y_values[i])
     *
     * @param x_values X-values of the list of points
     * @param y_values Y-values of the list of points
     * @param slope slope of the best fit line (calculated by this function)
     * @param offset offset of the best fit line (calculated by this function)
     * @param sample_count number of weights to use for fitting
     * @return true if linear fit could be computed, false if not

    bool compute_linear_fit(const std::vector<float>& x_values, const std::vector<float>& y_values, float& slope, float& offset, size_t sample_count=0)
    {
        if(sample_count==0 || sample_count > x_values.size()) sample_count = x_values.size();
        if (x_values.size() <= 1 || x_values.size() != y_values.size())
        {
            //TODO: throw error?
            return false;
        }
        float sx = 0; // Sum(x)
        float sy = 0; // Sum(y)
        float sxx = 0; // Sum(x*x)
        float sxy = 0; // Sum(x*y)
        for (size_t i = 0; i < sample_count; i++)
        {
            sx += x_values[i];
            sy += y_values[i];
            sxy += x_values[i] * y_values[i];
            sxx += x_values[i] * x_values[i];
        }
        const float denominator = sample_count * sxx - sx * sx;
        if (denominator > std::numeric_limits<float>::epsilon())
        {
            slope = (sample_count * sxy - sx * sy) / denominator;
            offset = (sy * sxx - sx * sxy) / denominator;
            return true;
        }
        return false;
    }*/

    private static float[] computeLinearFit(Float[] xValues, Float[] yValues, int sampleCount) {
        if (sampleCount == 0 || sampleCount > xValues.length) sampleCount = xValues.length;
        if (xValues.length <= 1 || xValues.length != yValues.length) {
            throw new PicardException("Can not compute linear fit.");
        }

        float sumX = 0;
        float sumY = 0;
        float sumXX = 0;
        float sumXY = 0;

        //returns
        float slope = 0f;
        float offset = 0f;

        for (int i = 0; i < sampleCount; i++) {
            sumX += xValues[i];
            sumY += yValues[i];
            sumXY += xValues[i] * yValues[i];
            sumXX += xValues[i] * xValues[i];
        }

        float denominator = sampleCount * sumXX - sumX * sumX;

        if (denominator > Math.ulp(denominator)) {
            slope = (sampleCount * sumXY - sumX * sumY) / denominator;
            offset = (sumY * sumXX - sumX * sumXY) / denominator;
        }
        return new float[]{slope, offset};
    }
}
