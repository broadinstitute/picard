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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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

        final List<PerTilePadHoppingMetricsExtractor> extractors = new ArrayList<PerTilePadHoppingMetricsExtractor>(factory.getAvailableTiles().size());
        for (final int tile : factory.getAvailableTiles()) {
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

        //PAD-HOPPING-SPECIFIC FIELDS
        //Need to initialize an array of all x/y/bases/hash data and the union-find data structures

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
                while (provider.hasNext()) {
                    final ClusterData cluster = provider.next();

                    //BEGINNING OF PAD-HOPPING-SPECIFIC CODE
                    //fill arrays of x/y/bases/hash data from cluster data
                    //union-find processing occurs after this while loop
                    // if (cluster.isPf()) . . . only deal with Pf reads??
                    //this.summaryMetric.READS++;

                }   //end of while -- all clusters have been processed

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


    /** a metric class for describing FP failing reads from an Illumina HiSeqX lane * */
    //PAD-HOPPING-SPECIFIC: ENLARGE THIS CLASS!!!!
    public static class PadHoppingDetailMetric extends MetricBase {
        // The Tile that is described by this metric.
        public Integer TILE;

        //The X coordinate of the read within the tile
        public int X;

        //The Y coordinate of the read within the tile
        public int Y;

        public PadHoppingDetailMetric(final Integer TILE, final int x, final int y) {
            this.TILE = TILE;
            X = x;
            Y = y;
        }

        /** This ctor is necessary for when reading metrics from file */
        public PadHoppingDetailMetric() {
        }
    }

    /**
     * Metrics produced by the GetHiSeqXPFFailMetrics program. Used to diagnose lanes from HiSeqX
     * Sequencing, providing the number and fraction of each of the reasons that reads could have not passed PF.
     * Possible reasons are EMPTY (reads from empty wells with no template strand), POLYCLONAL (reads from wells that had more than one strand
     * cloned in them), MISALIGNED (reads from wells that are near the edge of the tile), UNKNOWN (reads that didn't pass PF but couldn't be diagnosed)
     */
    //PAD-HOPPING-SPECIFIC: CHANGE THIS CLASS!!!!
    public static class PadHoppingSummaryMetric extends MetricBase {
        /** The Tile that is described by this metric. Can be a string (like "All") to mean some marginal over tiles. * */
        public String TILE = null;

        /** The total number of reads examined */
        public int READS = 0;

        /** The number of non-PF reads in this tile. */
        public int PF_FAIL_READS = 0;

        /** The fraction of PF_READS */
        public double PCT_PF_FAIL_READS = 0.0;

        /** The number of non-PF reads in this tile that are deemed empty. */
        public int PF_FAIL_EMPTY = 0;

        /** The fraction of non-PF reads in this tile that are deemed empty (as fraction of all non-PF reads). */
        public double PCT_PF_FAIL_EMPTY = 0.0;

        /** The number of non-PF reads in this tile that are deemed multiclonal. */
        public int PF_FAIL_POLYCLONAL = 0;

        /** The fraction of non-PF reads in this tile that are deemed multiclonal (as fraction of all non-PF reads). */
        public double PCT_PF_FAIL_POLYCLONAL = 0.0;

        /** The number of non-PF reads in this tile that are deemed "misaligned". */
        public int PF_FAIL_MISALIGNED = 0;

        /** The fraction of non-PF reads in this tile that are deemed "misaligned" (as fraction of all non-PF reads). */
        public double PCT_PF_FAIL_MISALIGNED = 0.0;

        /** The number of non-PF reads in this tile that have not been classified. */
        public int PF_FAIL_UNKNOWN = 0;

        /** The fraction of non-PF reads in this tile that have not been classified (as fraction of all non-PF reads). */
        public double PCT_PF_FAIL_UNKNOWN = 0.0;

        // constructor takes a String for tile since we want to have one instance with tile="All". This tile will contain the summary of all the tiles
        public PadHoppingSummaryMetric(final String tile) {
            TILE = tile;
        }

        /** This ctor is necessary for when reading metrics from file */
        public PadHoppingSummaryMetric() {
        }

        /**
         * Adds the non-calculated fields from the other metric to this one.
         *
         * @param metric
         */
        public void merge(final PadHoppingSummaryMetric metric) {
            this.READS += metric.READS;
            this.PF_FAIL_READS += metric.PF_FAIL_READS;
            this.PF_FAIL_EMPTY += metric.PF_FAIL_EMPTY;
            this.PF_FAIL_MISALIGNED += metric.PF_FAIL_MISALIGNED;
            this.PF_FAIL_POLYCLONAL += metric.PF_FAIL_POLYCLONAL;
            this.PF_FAIL_UNKNOWN += metric.PF_FAIL_UNKNOWN;
        }

        public void calculateDerivedFields() {
            //protect against divide by zero
            if (this.READS != 0) {
                this.PCT_PF_FAIL_READS = (double) this.PF_FAIL_READS / this.READS;
                this.PCT_PF_FAIL_EMPTY = (double) this.PF_FAIL_EMPTY / this.READS;
                this.PCT_PF_FAIL_MISALIGNED = (double) this.PF_FAIL_MISALIGNED / this.READS;
                this.PCT_PF_FAIL_POLYCLONAL = (double) this.PF_FAIL_POLYCLONAL / this.READS;
                this.PCT_PF_FAIL_UNKNOWN = (double) this.PF_FAIL_UNKNOWN / this.READS;
            }
        }
    }


}
