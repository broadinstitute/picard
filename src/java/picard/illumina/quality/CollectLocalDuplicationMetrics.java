/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.illumina.quality;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProvider;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * In Illumina flowcells, a "cluster" is a contiguous region containing
 * bridge-amplified clones of a DNA fragment.  In the following, groups of
 * duplicate clusters will be called "bunches".
 *
 * @author David Benjamin
 */
@CommandLineProgramProperties(
        usage = "Measure duplication between clusters within a specified maximum separation." +
                "  Input is an Illumina basecalls directory and this tool may be run before sequencing is complete." +
                "  One important application is pad-hopping in HiSeqX flowcells.",
        usageShort = "Measure local duplication in Illumina sequencing from a basecalls directory.",
        programGroup = Metrics.class
)
public class CollectLocalDuplicationMetrics extends CommandLineProgram {
    //Command line options in addition to those inherited from CommandLineProgram
    @Option(doc = "The Illumina basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Basename for metrics file. Resulting file will be" +
            " <OUTPUT>" + SUMMARY_METRICS_EXTENSION, optional = false)
    public File OUTPUT;

    @Option(doc = "The fraction of local duplicate bunches to output in detailed metrics.", optional = true)
    public double PROB_EXPLICIT_OUTPUT = 0;

    @Option(doc = "Lane number.", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Option(doc = "Run this many PerTileLocalDuplicationMetricsExtractor in parallel.  If NUM_PROCESSORS = 0, use all available cores. " +
            "If NUM_PROCESSORS < 0 use all but |NUM_PROCESSORS| cores.", optional = true)
    public int NUM_PROCESSORS = 1;

    @Option(doc = "Number of tiles on which to calculate metrics.  Default of 8 gives a good lane average.", optional = true)
    public int N_TILES = 8;

    @Option(doc = "Index of first tile (0 to 95).  Default -1 is an evenly-spaced sample over the lane", optional = true)
    public int TILE_INDEX = -1;

    @Option(doc = "Number of bases to look at.  Due to sequencing error comparing too many bases" +
            " may yield an underestimate of duplication.", optional = true)
    public int NUM_BASES = 24;

    @Option(doc = "Max distance in pixels between duplicate clusters in the same bunch." +
            "  If two clusters are farther than MAX_SEPARATION from one another but within MAX_SEPRATION " +
            "from a third cluster, all three belong to the same bunch.  Pads in the HiSeqX are " +
            "separated by 20 pixels and for the purpose of measuring pad-hopping an accurate and fast " +
            "approximation is to consider all duplicates on the same tile to be pad-hopping, hence the infinite default.", optional = true)
    public double MAX_SEPARATION = Double.POSITIVE_INFINITY;

    private static final Log LOG = Log.getInstance(CollectLocalDuplicationMetrics.class);

    //Set up a LocalDuplicationSummaryMetrics and a List of LocalDuplicationDetailMetrics for each tile
    private final Map<Integer, LocalDuplicationSummaryMetrics> tileToSummaryMetrics = new HashMap<Integer, LocalDuplicationSummaryMetrics>();
    private final Map<Integer, List<LocalDuplicationDetailMetrics>> tileToDetailedMetrics = new HashMap<Integer, List<LocalDuplicationDetailMetrics>>();

    //Add "T" to the number of cycles to create a "TemplateRead" of the desired length.
    private final ReadStructure READ_STRUCTURE = new ReadStructure(NUM_BASES + "T");

    public static final String DETAILED_METRICS_EXTENSION = "local_duplication_detailed_metrics";
    public static final String SUMMARY_METRICS_EXTENSION = "local_duplication_summary_metrics";

    public static final int TILES_PER_LANE = 96;

    //Add error-checking for the command line arguments specific to this program
    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<String>();

        if (NUM_BASES < 1) errors.add("Must consider at least one base (and really fewer than 6 is nonsensical)");

        if (N_TILES < 1 || N_TILES > TILES_PER_LANE) {
            errors.add(String.format("Must process at least 1 and at most %d tiles", TILES_PER_LANE));
        }

        if (TILE_INDEX < -1) errors.add("Must choose a non-negative tile index or -1 to take a lane average");
        if (TILE_INDEX >= TILES_PER_LANE ) {
            errors.add(String.format("Tile index may be at most %d", TILES_PER_LANE - 1));
        }

        if ( PROB_EXPLICIT_OUTPUT < 0 || PROB_EXPLICIT_OUTPUT > 1 ) {
            errors.add("PROB_EXPLICIT_OUTPUT must be a probability (from 0 to 1).");
        }

        if (errors.size() > 0) return errors.toArray(new String[errors.size()]);
        else return super.customCommandLineValidation();
    }

    /** Stock main method for any CommandLineProgram. */
    public static void main(final String[] args) { new CollectLocalDuplicationMetrics().instanceMainWithExit(args); }

    @Override
    protected int doWork() {
        /**
         * Each tile is processed on a single thread by a PerTileLocalDuplicationMetricsExtractor, which asks
         * the IlluminaDataProviderFactory for an IlluminaDataProvider, which is an iterator for all the
         * ClusterData on a single tile.  ClusterData contains the raw data of a read and its x-y coordinates.
         */
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, READ_STRUCTURE,
                new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY),
                IlluminaDataType.BaseCalls, IlluminaDataType.PF, IlluminaDataType.Position);

        final File summaryMetricsFileName = new File(OUTPUT + "." + SUMMARY_METRICS_EXTENSION);
        final File detailedMetricsFileName = new File(OUTPUT + "." + DETAILED_METRICS_EXTENSION);

        IOUtil.assertFileIsWritable(summaryMetricsFileName);
        if (PROB_EXPLICIT_OUTPUT > 0) IOUtil.assertFileIsWritable(detailedMetricsFileName);

        final int numProcessors = NUM_PROCESSORS + ((NUM_PROCESSORS > 0) ? 0 : Runtime.getRuntime().availableProcessors());
        final ExecutorService pool = Executors.newFixedThreadPool(numProcessors);
        LOG.info(String.format("Processing with %d threads.", numProcessors));

        final List<Integer> allTiles = new ArrayList<Integer>(factory.getAvailableTiles());
        Collections.sort(allTiles);

        List<Integer> tilesToProcess;
        //default case of evenly-spaced tiles to average over the lane
        if (TILE_INDEX == -1) {
            int offset = TILES_PER_LANE/(N_TILES * 2);
            tilesToProcess = new ArrayList<Integer>();
            for (int i = 0; i < N_TILES; i++) {
                tilesToProcess.add(allTiles.get(offset + ((TILES_PER_LANE * i) / N_TILES)));
            }
        }
        else {
            final int firstTile = TILE_INDEX;
            final int lastTile = Math.min(allTiles.size(), firstTile + N_TILES);
            tilesToProcess = allTiles.subList(firstTile, lastTile);
        }

        LOG.info(String.format("Computing metrics for %d tiles.", + tilesToProcess.size()));

        final List<PerTileLocalDuplicationMetricsExtractor> extractors = new ArrayList<PerTileLocalDuplicationMetricsExtractor>(tilesToProcess.size());
        for (final int tile : tilesToProcess) {
            tileToSummaryMetrics.put(tile, new LocalDuplicationSummaryMetrics(Integer.toString(tile)));
            tileToDetailedMetrics.put(tile, new ArrayList<LocalDuplicationDetailMetrics>());

            extractors.add(new PerTileLocalDuplicationMetricsExtractor(tile, tileToSummaryMetrics.get(tile),
                    tileToDetailedMetrics.get(tile), factory, PROB_EXPLICIT_OUTPUT, MAX_SEPARATION, NUM_BASES));
        }
        try {
            for (final PerTileLocalDuplicationMetricsExtractor extractor : extractors) pool.submit(extractor);
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch (final Throwable e) {
            // Cancel if current thread also interrupted
            LOG.error(e, "Problem submitting extractors to thread pool or awaiting shutdown of thread pool.  Attempting to kill thread pool.");
            pool.shutdownNow();
            return 2;
        }
        LOG.info(String.format("Processed all %d tiles.", extractors.size()));

        // Check for exceptions from extractors
        for (final PerTileLocalDuplicationMetricsExtractor extractor : extractors) {
            if (extractor.getException() != null) {
                LOG.error("Abandoning calculation because one or more PerTileLocalDuplicationMetricsExtractors failed.");
                return 4;
            }
        }

        final MetricsFile<LocalDuplicationDetailMetrics, ?> detailedMetrics = getMetricsFile();
        for (final Collection<LocalDuplicationDetailMetrics> detailedMetricCollection : tileToDetailedMetrics.values()) {
            for (final LocalDuplicationDetailMetrics metrics : detailedMetricCollection) {
                detailedMetrics.addMetric(metrics);
            }
        }

        if (PROB_EXPLICIT_OUTPUT > 0) detailedMetrics.write(detailedMetricsFileName);

        final LocalDuplicationSummaryMetrics totalMetrics = new LocalDuplicationSummaryMetrics("All"); // a "fake" tile that will contain the total tally
        for (final LocalDuplicationSummaryMetrics summaryMetrics : tileToSummaryMetrics.values()) {
            totalMetrics.merge(summaryMetrics);
        }
        totalMetrics.calculateDerivedFields();
        final MetricsFile<LocalDuplicationSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
        summaryMetricsFile.addMetric(totalMetrics);

        for (final LocalDuplicationSummaryMetrics summaryMetrics : tileToSummaryMetrics.values()) {
            summaryMetrics.calculateDerivedFields();
            summaryMetricsFile.addMetric(summaryMetrics);
        }
        summaryMetricsFile.write(summaryMetricsFileName);

        return 0;
    }

    /** Extracts metrics from a single tile on its own thread
     */
    private class PerTileLocalDuplicationMetricsExtractor implements Runnable {

        private final int tile;
        private final LocalDuplicationSummaryMetrics summaryMetrics;
        final Collection<LocalDuplicationDetailMetrics> detailedMetrics;
        private Exception exception = null;
        private final IlluminaDataProvider provider;
        private final double probWriteDetailed;
        private final double maxSeparation;
        private final int nBases;
        private final Random random = new Random(1);

        public PerTileLocalDuplicationMetricsExtractor(final int tile, final LocalDuplicationSummaryMetrics summaryMetrics,
                final Collection<LocalDuplicationDetailMetrics> detailedMetrics, final IlluminaDataProviderFactory factory,
                final double probWriteDetailed, final double maxSeparation, final int nBases) {
            this.tile = tile;
            this.summaryMetrics = summaryMetrics;
            this.detailedMetrics = detailedMetrics;
            this.probWriteDetailed = probWriteDetailed;
            this.maxSeparation = maxSeparation;
            this.nBases = nBases;
            this.provider = factory.makeDataProvider(Arrays.asList(tile));
        }

        public Exception getException() { return this.exception; }

        /** run method which extracts accumulates metrics for a tile */
        public void run() {
            try {
                LOG.info(String.format("Extracting metrics for tile %d.", tile));

                Map<String, List<Point>> duplicateSets = new HashMap<String, List<Point>>();

                for (final ClusterData cluster : provider) {
                    if (! cluster.isPf() ) continue;
                    summaryMetrics.READS++;

                    //getBases() returns byte[]. Converting to String costs some speed but makes hashing easier
                    final String allBases = new String(cluster.getRead(0).getBases());
                    final String bases = allBases.substring(0, nBases);

                    List<Point> list = duplicateSets.get(bases);
                    if (list == null) duplicateSets.put(bases, list = new ArrayList<Point>());
                    list.add(new Point(cluster.getX(), cluster.getY()));
                }

                for (final Map.Entry<String, List<Point>> entry : duplicateSets.entrySet()) {
                    final List<Point> points = entry.getValue();
                    if (points.size() == 1) continue; //if there is no duplication
                    final String bases = entry.getKey();

                    if (maxSeparation == Double.POSITIVE_INFINITY) {
                        summaryMetrics.LOCAL_DUPLICATES += points.size() - 1;
                        if (random.nextDouble() < probWriteDetailed) {
                            detailedMetrics.add(new LocalDuplicationDetailMetrics(tile, bases, points));
                        }
                    }
                    else {
                        BunchFinder bunchFinder = new BunchFinder(points, maxSeparation);
                        for (final Bunch bunch : bunchFinder.getBunches()) {
                            if (bunch.size() == 1) continue;
                            summaryMetrics.LOCAL_DUPLICATES += bunch.numDuplicates();
                            if (random.nextDouble() < probWriteDetailed) {
                                detailedMetrics.add(new LocalDuplicationDetailMetrics(tile, bases, bunch));
                            }
                        }
                    }
                }

            } catch (final Exception e) {
                LOG.error(e, String.format("Error processing tile %d.", tile));
                this.exception = e;
            } finally {
                provider.close();
            }
        }
    }

    private class Bunch extends ArrayList<Point> {
        public int numDuplicates() { return size() - 1; }
    }

    private class BunchFinder {
        private ArrayList<Bunch> bunches;
        private int numPoints;

        public BunchFinder(List<Point> points, double maxSeparation) {
            bunches = new ArrayList<Bunch>();
            numPoints = points.size();
            boolean[] visited = new boolean[numPoints];

            for (int root = 0; root < numPoints; root++) {
                if (visited[root]) continue;   //point belongs to a previously-counted component
                Bunch bunch = new Bunch();

                //depth-first search for all points in same Bunch as root
                Stack<Integer> pointStack = new Stack<Integer>();
                pointStack.push(root);
                while (!pointStack.isEmpty()) {
                    int bud = pointStack.pop();
                    bunch.add(points.get(bud));
                    for (int shoot = bud + 1; shoot < numPoints; shoot++) {
                        if (!visited[shoot] && points.get(bud).distance(points.get(shoot)) <= maxSeparation) {
                            pointStack.push(shoot);
                            visited[shoot] = true;
                        }
                    }
                }
                bunches.add(bunch);
            }
        }

        public List<Bunch> getBunches() { return bunches; }
    }

}


