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

import htsjdk.samtools.metrics.MetricBase;
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

import java.awt.Point;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Collect metrics regarding pad-hopping in the Illumina Hi-Seq X.
 * Terminology: in all Illumina Hi-Seq machines, a "cluster" is a discrete area on a flow cell that contains
 * bridge-amplified DNA fragments, hopefully all clones of a single original insert  In the Hi-Seq X, each cluster
 * is constrained within a hexagonal "pad".  In "pad-hopping", contiguous groups of pads have duplicate DNA fragments.
 * In the following code, such a contiguous group will be called a "bunch".
 *
 * @author David Benjamin
 */
@CommandLineProgramProperties(
        usage = "Measure pad-hopping duplication in HiSeqX.",
        usageShort = "Measure pad-hopping duplication in HiSeqX.",
        programGroup = Metrics.class
)

public class CollectPadHoppingMetrics extends CommandLineProgram {
    //Command line options in addition to those inherited from CommandLineProgram
    @Option(doc = "The Illumina basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Basename for metrics file. Resulting file will be" +
            " <OUTPUT>" + summaryMetricsExtension, optional = false)
    public File OUTPUT;

    @Option(shortName = "P", doc = "The fraction of pad-hopping events to output explicitly. Output file will be <OUTPUT>" + detailedMetricsExtension + " (if PROB_EXPLICIT_OUTPUT != 0)", optional = true)
    public double PROB_EXPLICIT_OUTPUT = 0;

    @Option(doc = "Lane number.", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Option(shortName = "NP", doc = "Run this many PerTilePadHoppingMetricsExtractor in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.", optional = true)
    public int NUM_PROCESSORS = 1;

    @Option(shortName = "TP", doc = "Number of tiles on which to calculate pad-hopping metrics.", optional = true)
    public int TILES_TO_PROCESS = 1;

    @Option(shortName = "TI", doc = "Index of first tile (0 to 95)", optional = true)
    public int TILE_INDEX = 0;

    @Option(doc = "Number of cycles to look at. At time of writing PF status gets determined at cycle 24 so numbers greater than this will yield strange results. " +
            "In addition, PF status is currently determined at cycle 24, so running this with any other value is neither tested nor recommended.", optional = true)
    public int N_CYCLES = 24;

    @Option(shortName = "NB", doc = "Number of bases to use for comparing reads.  We might make this less" +
            " than N_CYCLES to avoid masking duplicates by sequencing errors.", optional = true)
    public int N_BASES = 24;

    @Option(shortName = "ND", doc = "Max distance (in Illumina's internal cluster coordinate units) for two custers " +
            "to be considered adjacent.  The distance is 20 +/- 1 for all tiles of all Hi Seq X flowcells.", optional = true)
    public double MAX_NEIGHBOR_DISTANCE = 22.0;

    private static final Log LOG = Log.getInstance(CollectPadHoppingMetrics.class);

    //Set up a PadHoppingSummaryMetric and a List of PadHoppingDetailMetrics for each tile
    private final Map<Integer, PadHoppingSummaryMetric> tileToSummaryMetrics = new LinkedHashMap<Integer, PadHoppingSummaryMetric>();
    private final Map<Integer, List<PadHoppingDetailMetric>> tileToDetailedMetrics = new LinkedHashMap<Integer, List<PadHoppingDetailMetric>>();

    //Add "T" to the number of cycles to create a "TemplateRead" of the desired length.
    private final ReadStructure READ_STRUCTURE = new ReadStructure(N_CYCLES + "T");

    public final static String detailedMetricsExtension = ".pad_hopping_detailed_metrics";
    public final static String summaryMetricsExtension = ".pad_hopping_summary_metrics";

    //Add error-checking for the command line arguments specific to this program
    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<String>();

        if (N_CYCLES < 0)
            errors.add("Number of Cycles to look at must be greater than 0");

        if (TILES_TO_PROCESS < 1)
            errors.add("Must process at least one tile");

        if (TILES_TO_PROCESS > 96)
            errors.add("There are only 96 tiles per lane of the Illumina HiSeqX");

        if (TILE_INDEX < 0)
            errors.add("Must choose a non-negative tile index");

        if (TILE_INDEX > 95)
            errors.add("Tile index may be at most 95 (there are 96 tiles on the HiSeqX)");

        if (PROB_EXPLICIT_OUTPUT > 1 || PROB_EXPLICIT_OUTPUT < 0)
            errors.add("PROB_EXPLICIT_OUTPUT must be a probability, i.e., 0 <= PROB_EXPLICIT_OUTPUT <= 1");

        if (errors.size() > 0) {
            return errors.toArray(new String[errors.size()]);
        } else {
            return super.customCommandLineValidation();
        }
    }

    /** Stock main method for any CommandLineProgram. */
    public static void main(final String[] args) { new CollectPadHoppingMetrics().instanceMainWithExit(args); }

    @Override
    protected int doWork() {
        /**
         * Each tile is processed on a single thread by a PerTilePadHoppingMetricsExtractor, which asks
         * the IlluminaDataProviderFactory for an IlluminaDataProvider, which is an iterator for all the
         * ClusterData on a single tile.  ClusterData contains the raw data of a read and its x-y coordinates.
         */
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, READ_STRUCTURE,
                new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY),
                IlluminaDataType.BaseCalls, IlluminaDataType.PF, IlluminaDataType.QualityScores, IlluminaDataType.Position);

        final File summaryMetricsFileName = new File(OUTPUT + summaryMetricsExtension);
        final File detailedMetricsFileName = new File(OUTPUT + detailedMetricsExtension);

        IOUtil.assertFileIsWritable(summaryMetricsFileName);
        if (PROB_EXPLICIT_OUTPUT > 0) IOUtil.assertFileIsWritable(detailedMetricsFileName);

        final int numProcessors = NUM_PROCESSORS + ((NUM_PROCESSORS > 0) ? 0 : Runtime.getRuntime().availableProcessors());
        final ExecutorService pool = Executors.newFixedThreadPool(numProcessors);
        LOG.info("Processing with " + numProcessors + " PerTilePadHoppingMetricsExtractor(s).");

        List<Integer> allTiles = new ArrayList<Integer>(factory.getAvailableTiles());
        Collections.sort(allTiles);

        int firstTile = TILE_INDEX;
        int lastTile = Math.min(allTiles.size(), firstTile + TILES_TO_PROCESS);
        final List<Integer> tilesToProcess = allTiles.subList(firstTile, lastTile);
        LOG.info("Computing pad hopping metrics for " + tilesToProcess.size() + " tiles.");

        final List<PerTilePadHoppingMetricsExtractor> extractors = new ArrayList<PerTilePadHoppingMetricsExtractor>(tilesToProcess.size());
        for (final int tile : tilesToProcess) {
            tileToSummaryMetrics.put(tile, new PadHoppingSummaryMetric(Integer.toString(tile)));
            tileToDetailedMetrics.put(tile, new ArrayList<PadHoppingDetailMetric>());

            extractors.add(new PerTilePadHoppingMetricsExtractor(tile, tileToSummaryMetrics.get(tile),
                    tileToDetailedMetrics.get(tile), factory, PROB_EXPLICIT_OUTPUT, MAX_NEIGHBOR_DISTANCE));
        }
        try {
            for (final PerTilePadHoppingMetricsExtractor extractor : extractors)
                pool.submit(extractor);
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch (final Throwable e) {
            // Cancel if current thread also interrupted
            LOG.error(e, "Problem submitting extractors to thread pool or awaiting shutdown of thread pool.  Attempting to kill thread pool.");
            pool.shutdownNow();
            return 2;
        }
        LOG.info("Processed all " + extractors.size() + " tiles.");

        // Check for exceptions from extractors
        for (final PerTilePadHoppingMetricsExtractor extractor : extractors) {
            if (extractor.getException() != null) {
                LOG.error("Abandoning calculation because one or more PerTilePadHoppingMetricsExtractors failed.");
                return 4;
            }
        }

        final MetricsFile<PadHoppingDetailMetric, ?> detailedMetrics = getMetricsFile();
        for (final Collection<PadHoppingDetailMetric> detailedMetricCollection : tileToDetailedMetrics.values())
            for (final PadHoppingDetailMetric metric : detailedMetricCollection)
                detailedMetrics.addMetric(metric);

        if (PROB_EXPLICIT_OUTPUT > 0)
            detailedMetrics.write(detailedMetricsFileName);

        final PadHoppingSummaryMetric totalMetric = new PadHoppingSummaryMetric("All"); // a "fake" tile that will contain the total tally
        for (final PadHoppingSummaryMetric summaryMetric : tileToSummaryMetrics.values())
            totalMetric.merge(summaryMetric);
        totalMetric.calculateDerivedFields();
        final MetricsFile<PadHoppingSummaryMetric, ?> summaryMetricsFile = getMetricsFile();
        summaryMetricsFile.addMetric(totalMetric);

        for (final PadHoppingSummaryMetric summaryMetric : tileToSummaryMetrics.values()) {
            summaryMetric.calculateDerivedFields();
            summaryMetricsFile.addMetric(summaryMetric);
        }
        summaryMetricsFile.write(summaryMetricsFileName);

        return 0;
    }

    /** Extracts metrics from a HiSeqX tile.
     * Different tiles use different files so each Extractor has its own thread to avoid waiting for
     * each other's file I/O.
     */
    private static class PerTilePadHoppingMetricsExtractor implements Runnable {

        private final int tile;
        private final PadHoppingSummaryMetric summaryMetric;
        final Collection<PadHoppingDetailMetric> detailedMetrics;
        private Exception exception = null;
        private final IlluminaDataProvider provider;
        final private double pWriteDetailed;
        final private double cutoffDistance;
        final private Random random = new Random();

        public PerTilePadHoppingMetricsExtractor(final int tile, final PadHoppingSummaryMetric summaryMetric,
                final Collection<PadHoppingDetailMetric> detailedMetrics, final IlluminaDataProviderFactory factory,
                final double pWriteDetailed, final double cutoffDistance) {
            this.tile = tile;
            this.summaryMetric = summaryMetric;
            this.detailedMetrics = detailedMetrics;
            this.pWriteDetailed = pWriteDetailed;
            this.cutoffDistance = cutoffDistance;
            this.provider = factory.makeDataProvider(Arrays.asList(tile));
        }

        public Exception getException() { return this.exception; }

        /** run method which extracts accumulates metrics for a tile */
        public void run() {
            try {
                LOG.info("Extracting pad-hopping metrics for tile " + tile);

                //a possible source of improved performance is that this sets up a List<Point> for
                //EVERY read, even ones that are not duplicated
                Map<String, List<Point>> duplicateSets = new HashMap<String, List<Point>>();
                for (final ClusterData cluster : provider) {
                    if (! cluster.isPf() ) continue; //only deal with PF reads
                    summaryMetric.READS++;

                    //getBases() returns byte[].  Converting to String loses performance but is more convenient for hashing
                    //Someone who knows Java better could probably advise me on a better method
                    final String bases = new String(cluster.getRead(0).getBases());

                    List<Point> list = duplicateSets.get(bases);
                    if (list == null)
                        duplicateSets.put(bases, list = new ArrayList<Point>());
                    list.add(new Point(cluster.getX(), cluster.getY()));
                }
                for (Map.Entry<String, List<Point>> entry : duplicateSets.entrySet()) {
                    List<Point> points = entry.getValue();
                    String bases = entry.getKey();
                    if (points.size() > 1) {    //if there is duplication
                        BunchFinder bunchFinder = new BunchFinder(bases, points, cutoffDistance);
                        for (Bunch bunch : bunchFinder.getBunches()) {
                            if (bunch.size() == 1) continue;
                            summaryMetric.PAD_HOPPING_DUPLICATES += bunch.numDuplicates();
                            //randomly add pad-hopping events to detailed metrics
                            if (random.nextDouble() < pWriteDetailed) {
                                Point center = bunch.center();
                                detailedMetrics.add(new PadHoppingDetailMetric(tile, bases, center.getX(), center.getY(), bunch.size()));
                            }
                        }
                    }
                }

            } catch (final Exception e) {
                LOG.error(e, "Error processing tile ", this.tile);
                this.exception = e;
            } finally {
                provider.close();
            }
        }
    }

    /**
     * cluster.getRead(0).getBases() returns a byte[], not a String.  A problem with this is that byte[] hashing
     * is by object identity, not by value.  The simplest fix is just to use the String(byte[]) constructor, but
     * that is really wasteful
     */
    private class BasesWrapper
    {
        private final byte[] bases;

        public BasesWrapper(byte[] data) { this.bases = data; }

        @Override
        public boolean equals(Object other) {
            if ( !(other instanceof BasesWrapper) ) return false;
            for (int i = 0; i < N_BASES; i++) {
                if ( ((BasesWrapper)other).bases[i] != bases[i]) return false;
            }
            return true;
        }

        @Override
        public int hashCode() {
            int hash = 173; // arbitrary seed value
            int multiplier = 37; // arbitrary multiplier value
            for (int i = 0; i < N_BASES; i++) {
                hash = hash * multiplier + bases[i];
            }
            return hash;
        }
    }


    /**
     * A Bunch is little more than a typedef for a list of Points.  It contains a few extra methods for characterizing
     * pad-hopping and stores the read common to all clusters in the bunch.
     *
     * Depending on how much we deeply we wish to study pad-hopping, we could add more methods.
     */
    private static class Bunch extends ArrayList<Point> {
        private final String bases;

        public Bunch(String s) { bases = s; }

        public String getBases() { return bases;}

        public int numDuplicates() { return size() - 1; }

        public Point center() {
            int totalX = 0;
            int totalY = 0;
            for (Point p : this) {
                totalX += p.getX();
                totalY += p.getY();
            }
            return new Point(totalX / size(), totalY / size());
        }
    }

    private static class BunchFinder {
        private ArrayList<Bunch> bunches;
        private int N;  //total number of points

        public BunchFinder(String bases, List<Point> points, double cutoffDistance) {
            bunches = new ArrayList<Bunch>();
            N = points.size();
            boolean[] visited = new boolean[N];

            for (int root = 0; root < N; root++) {
                if (visited[root]) continue;   //point belongs to a previously-counted component
                Bunch bunch = new Bunch(bases);

                //depth-first search for all points in same Bunch as root
                Stack<Integer> DFS = new Stack<Integer>();
                DFS.push(root);
                while (!DFS.isEmpty()) {
                    int bud = DFS.pop();
                    bunch.add(points.get(bud));
                    for (int shoot = bud + 1; shoot < N; shoot++) {
                        if (!visited[shoot] && points.get(bud).distance(points.get(shoot)) <= cutoffDistance) {
                            DFS.push(shoot);
                            visited[shoot] = true;
                        }
                    }
                }
                bunches.add(bunch);
            }
        }

        public List<Bunch> getBunches() { return bunches; }
    }

    /** a metric class for describing pad-hopping bunches **/
    public static class PadHoppingDetailMetric extends MetricBase {
        public Integer TILE;

        public String BASES;

        //The X and Y coordinates of the center of the bunch
        public double X;
        public double Y;

        //The number of clusters in this pad-hopping bunch
        public int SIZE;

        public PadHoppingDetailMetric(final Integer tile, final String bases, final double x, final double y, final int size) {
            TILE = tile;
            BASES = bases;
            X = x;
            Y = y;
            SIZE = size;
        }

        /** This constructor is necessary for reading metrics from file */
        public PadHoppingDetailMetric() { }
    }

    public static class PadHoppingSummaryMetric extends MetricBase {
        /** The Tile that is described by this metric. Can be a string (like "All") to mean some marginal over tiles. * */
        public String TILE = null;

        /** The total number of reads examined */
        public int READS = 0;

        /** Duplicates due to pad-hopping in this tile.  In a bunch of N clusters, N - 1 are duplicates */
        public int PAD_HOPPING_DUPLICATES = 0;

        public double PCT_PAD_HOPPING_DUPLICATES = 0.0;

        // constructor takes a String for tile since we want to have one instance with tile="All". This tile will contain the summary of all the tiles
        public PadHoppingSummaryMetric(final String tile) {
            TILE = tile;
        }

        /** This constructor is necessary for when reading metrics from file */
        public PadHoppingSummaryMetric() { }

        public void merge(final PadHoppingSummaryMetric metric) {
            READS += metric.READS;
            PAD_HOPPING_DUPLICATES += metric.PAD_HOPPING_DUPLICATES;
        }

        public void calculateDerivedFields() {
            if (READS != 0)
                PCT_PAD_HOPPING_DUPLICATES = (double) PAD_HOPPING_DUPLICATES / this.READS;
        }
    }

}
