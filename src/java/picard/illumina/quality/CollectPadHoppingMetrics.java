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
import picard.illumina.parser.ReadData;
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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Collect metrics regarding the reason for reads (sequenced by HiSeqX) not passing the Illumina PF Filter. (BETA)
 *
 * @author David Benjamin
 */
@CommandLineProgramProperties(
        usage = "Detect pad-hopping duplication in HiSeqX.",
        usageShort = "Detect pad-hopping duplication in HiSeqX.",
        programGroup = Metrics.class
)

public class CollectPadHoppingMetrics extends CommandLineProgram {
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
    public int TILES_TO_PROCESS = Integer.MAX_VALUE;

    @Option(doc = "Number of cycles to look at. At time of writing PF status gets determined at cycle 24 so numbers greater than this will yield strange results. " +
            "In addition, PF status is currently determined at cycle 24, so running this with any other value is neither tested nor recommended.", optional = true)
    public int N_CYCLES = 24;

    private static final Log LOG = Log.getInstance(CollectPadHoppingMetrics.class);

    private final Map<Integer, PadHoppingSummaryMetric> tileToSummaryMetrics = new LinkedHashMap<Integer, PadHoppingSummaryMetric>();
    private final Map<Integer, List<PadHoppingDetailMetric>> tileToDetailedMetrics = new LinkedHashMap<Integer, List<PadHoppingDetailMetric>>();

    //Add "T" to the number of cycles to create a "TemplateRead" of the desired length.
    private final ReadStructure READ_STRUCTURE = new ReadStructure(N_CYCLES + "T");

    public final static String detailedMetricsExtension = ".pad_hopping_detailed_metrics";
    public final static String summaryMetricsExtension = ".pad_hopping_summary_metrics";

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<String>();

        if (N_CYCLES < 0) {
            errors.add("Number of Cycles to look at must be greater than 0");
        }

        if (PROB_EXPLICIT_OUTPUT > 1 || PROB_EXPLICIT_OUTPUT < 0) {
            errors.add("PROB_EXPLICIT_OUTPUT must be a probability, i.e., 0 <= PROB_EXPLICIT_OUTPUT <= 1");
        }

        if (errors.size() > 0) {
            return errors.toArray(new String[errors.size()]);
        } else {
            return super.customCommandLineValidation();
        }
    }

    /** Stock main method. */
    public static void main(final String[] args) {
        new CollectPadHoppingMetrics().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {

        //A factory gives a provider for every tile in the flowcell lane
        //A provider is an iterator for ClusterData
        //An extractor is a Runnable wrapper for a provider
        //ClusterData corresponds to a single BAM record (read)
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, READ_STRUCTURE,
                new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY),
                IlluminaDataType.BaseCalls,
                IlluminaDataType.PF,
                IlluminaDataType.QualityScores,
                IlluminaDataType.Position);

        final File summaryMetricsFileName = new File(OUTPUT + summaryMetricsExtension);
        final File detailedMetricsFileName = new File(OUTPUT + detailedMetricsExtension);

        IOUtil.assertFileIsWritable(summaryMetricsFileName);
        if (PROB_EXPLICIT_OUTPUT != 0) {
            IOUtil.assertFileIsWritable(detailedMetricsFileName);
        }

        final int numProcessors;
        if (NUM_PROCESSORS == 0) {
            numProcessors = Runtime.getRuntime().availableProcessors();
        } else if (NUM_PROCESSORS < 0) {
            numProcessors = Runtime.getRuntime().availableProcessors() + NUM_PROCESSORS;
        } else {
            numProcessors = NUM_PROCESSORS;
        }

        // Create thread-pool submit jobs and what for their completion
        LOG.info("Processing with " + numProcessors + " PerTilePadHoppingMetricsExtractor(s).");
        final ExecutorService pool = Executors.newFixedThreadPool(numProcessors);

        //get a random subset of all tiles
        final List<Integer> allTiles = factory.getAvailableTiles();
        Collections.shuffle(allTiles);
        List<Integer> tilesToProcess = new ArrayList<Integer>();
        for (int n = 0; n < allTiles.size() && n < TILES_TO_PROCESS; n++)
            tilesToProcess.add(allTiles.get(n));

        final List<PerTilePadHoppingMetricsExtractor> extractors = new ArrayList<PerTilePadHoppingMetricsExtractor>(tilesToProcess.size());
        for (final int tile : tilesToProcess) {
            tileToSummaryMetrics.put(tile, new PadHoppingSummaryMetric(Integer.toString(tile)));
            tileToDetailedMetrics.put(tile, new ArrayList<PadHoppingDetailMetric>());

            final PerTilePadHoppingMetricsExtractor extractor = new PerTilePadHoppingMetricsExtractor(
                    tile,
                    tileToSummaryMetrics.get(tile),
                    tileToDetailedMetrics.get(tile),
                    factory,
                    PROB_EXPLICIT_OUTPUT
            );
            extractors.add(extractor);
        }
        try {
            for (final PerTilePadHoppingMetricsExtractor extractor : extractors) {
                pool.submit(extractor);
            }
            pool.shutdown();
            // Wait forever for tasks to terminate
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch (final Throwable e) {
            // Cancel if current thread also interrupted
            LOG.error(e, "Parent thread encountered problem submitting extractors to thread pool or awaiting shutdown of threadpool.  Attempting to kill threadpool.");
            pool.shutdownNow();
            return 2;
        }

        LOG.info("Processed " + extractors.size() + " tiles.");

        // Check for exceptions from extractors
        for (final PerTilePadHoppingMetricsExtractor extractor : extractors) {
            if (extractor.getException() != null) {
                LOG.error("Abandoning metrics calculation because one or more PerTilePFMetricsExtractors failed.");
                return 4;
            }
        }

        // Add detailed metrics to file
        final MetricsFile<PadHoppingDetailMetric, ?> detailedMetrics = getMetricsFile();
        for (final Collection<PadHoppingDetailMetric> detailedMetricCollection : tileToDetailedMetrics.values()) {
            for (final PadHoppingDetailMetric metric : detailedMetricCollection) {
                detailedMetrics.addMetric(metric);
            }
        }

        // If detailed metrics were requested, write them now.
        if (PROB_EXPLICIT_OUTPUT > 0) {
            detailedMetrics.write(detailedMetricsFileName);
        }

        // Finish metrics tallying. Looping twice so that the "All" metrics will come out on top.
        final PadHoppingSummaryMetric totalMetric = new PadHoppingSummaryMetric("All"); // a "fake" tile that will contain the total tally
        for (final PadHoppingSummaryMetric summaryMetric : tileToSummaryMetrics.values()) {
            totalMetric.merge(summaryMetric);
        }

        // Derive fields for total metric and add to file
        totalMetric.calculateDerivedFields();
        final MetricsFile<PadHoppingSummaryMetric, ?> summaryMetricsFile = getMetricsFile();
        summaryMetricsFile.addMetric(totalMetric);

        // Prepare each tile's derived fields and add it to the file
        for (final PadHoppingSummaryMetric summaryMetric : tileToSummaryMetrics.values()) {
            summaryMetric.calculateDerivedFields();
            summaryMetricsFile.addMetric(summaryMetric);
        }

        // Actually write the summary metrics to their file.
        summaryMetricsFile.write(summaryMetricsFileName);

        return 0;
    }

    /** Extracts metrics from a HiSeqX tile. */
    private static class PerTilePadHoppingMetricsExtractor implements Runnable {

        private final int tile;
        private final PadHoppingSummaryMetric summaryMetric;
        final Collection<PadHoppingDetailMetric> detailedMetrics;
        private Exception exception = null;
        private final IlluminaDataProvider provider;
        final private double pWriteDetailed;
        final private Random random = new Random();

        /**
         * Constructor
         *
         * @param tile The number of the tile being processed.
         * @param summaryMetric A summaryMetric for collecting the tile data in.
         * @param detailedMetrics A set of metrics for collecting the classification data in.
         * @param factory A dataprovider for IlluminaData
         */
        public PerTilePadHoppingMetricsExtractor(
                final int tile,
                final PadHoppingSummaryMetric summaryMetric,
                final Collection<PadHoppingDetailMetric> detailedMetrics,
                final IlluminaDataProviderFactory factory,
                final double pWriteDetailed
        ) {
            this.tile = tile;
            this.summaryMetric = summaryMetric;
            this.detailedMetrics = detailedMetrics;
            this.pWriteDetailed = pWriteDetailed;
            this.provider = factory.makeDataProvider(Arrays.asList(tile));
        }

        public Exception getException() { return this.exception; }

        /** run method which extracts accumulates metrics for a tile */
        public void run() {
            try {
                LOG.info("Extracting pad-hopping metrics for tile " + tile);

                /**
                 *   Sometimes makeDataProvider takes a while waiting for slow file IO, for each tile the needed set of files
                 *   is non-overlapping sets of files so make the data providers in the individual threads for Extractors
                 *   so they are not all waiting for each others file operations
                 */
                Map<String, List<Point>> clusters = new HashMap<String, List<Point>>();
                while (provider.hasNext()) {
                    final ClusterData cluster = provider.next();
                    // if (cluster.isPf()) . . . only deal with Pf reads??
                    this.summaryMetric.READS++;

                    //getBases() returns byte[].  Converting to String loses performance but
                    //java arrays' hashing and comparison is by object identity, not value
                    //Casting to String is simpler than writing a byte[] wrapper with overridden hashCode()
                    final String bases = new String(cluster.getRead(0).getBases());

                    List<Point> list = clusters.get(bases);
                    if (list == null)
                        clusters.put(bases, list = new ArrayList<Point>());
                    list.add(new Point(cluster.getX(), cluster.getY()));
                }

                for (List<Point> points : clusters.values())
                {
                    if (points.size() > 1) {    //if there is duplication
                        //this.summaryMetric.DUPLICATE_READS += points.size();

                        //partition points into pad-hopping blobs
                        //I should have a nested class that takes a list of points in its constructor
                        //along with a cutoff distance or distance squared if I need to optimize
                    }
                }

                //MORE PAD-HOPPING-SPECIFIC CODE
                //do union-find processing

                //randomly add pad-hopping events to detailed metrics
                //if (random.nextDouble() < pWriteDetailed) {
                //    detailedMetrics.add(new PadHoppingDetailMetric(. . .));
                //}


            } catch (final Exception e) {
                LOG.error(e, "Error processing tile ", this.tile);
                this.exception = e;
            } finally {
                provider.close();
            }
        }
    }

    private static class NeighboringPoints {

        private ArrayList<ArrayList<Point>> groups;

        public NeighboringPoints(List<Point> points) {
            groups = new ArrayList<ArrayList<Point>>();

            int N = points.size();
            boolean[] visited = new boolean[N];


            for (int n = 0; n < N; n++)
                if (!visited[n]) {
                    visited[n] = true;
                    ArrayList<Point> group = new ArrayList<Point>();
                    group.add
                }

        }
    }

    /** a metric class for describing pad-hopping clusters **/
    public static class PadHoppingDetailMetric extends MetricBase {
        // The Tile that is described by this metric.
        public Integer TILE;

        //The X coordinate of the first read in the cluster within the tile
        public int X;

        //The Y coordinate of the read in the cluster within the tile
        public int Y;

        //The number of reads in this pad-hopping cluster
        public int SIZE;

        public PadHoppingDetailMetric(final Integer TILE, final int x, final int y, final int size) {
            this.TILE = TILE;
            X = x;
            Y = y;
            SIZE = size;
        }

        /** This ctor is necessary for when reading metrics from file */
        public PadHoppingDetailMetric() { }
    }

    public static class PadHoppingSummaryMetric extends MetricBase {
        /** The Tile that is described by this metric. Can be a string (like "All") to mean some marginal over tiles. * */
        public String TILE = null;

        /** The total number of reads examined */
        public int READS = 0;

        /** The number of reads contained in pad-hopping clusters in this tile. */
        public int PAD_HOPPING_READS = 0;

        public double PCT_PAD_HOPPING_READS = 0.0;

        // constructor takes a String for tile since we want to have one instance with tile="All". This tile will contain the summary of all the tiles
        public PadHoppingSummaryMetric(final String tile) {
            TILE = tile;
        }

        /** This ctor is necessary for when reading metrics from file */
        public PadHoppingSummaryMetric() { }

        /**
         * Adds the non-calculated fields from the other metric to this one.
         *
         * @param metric
         */
        public void merge(final PadHoppingSummaryMetric metric) {
            this.READS += metric.READS;
            this.PAD_HOPPING_READS += metric.PAD_HOPPING_READS;
        }

        public void calculateDerivedFields() {
            //protect against divide by zero
            if (this.READS != 0) {
                this.PCT_PAD_HOPPING_READS = (double) this.PAD_HOPPING_READS / this.READS;
            }
        }
    }


}
