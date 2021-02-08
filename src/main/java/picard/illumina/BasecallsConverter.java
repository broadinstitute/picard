package picard.illumina;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.PicardException;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.ThreadPoolExecutorUtil;
import picard.util.ThreadPoolExecutorWithExceptions;

import java.io.File;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

/**
 * BasecallsConverter utilizes an underlying IlluminaDataProvider to convert parsed and decoded sequencing data
 * from standard Illumina formats to specific output records (FASTA records/SAM records).
 * <p>
 * The underlying IlluminaDataProvider apply several optional transformations that can include EAMSS filtering,
 * non-PF read filtering and quality score recoding using a BclQualityEvaluationStrategy.
 * <p>
 * The converter can also limit the scope of data that is converted from the data provider by setting the
 * tile to start on (firstTile) and the total number of tiles to process (tileLimit).
 * <p>
 * Additionally, BasecallsConverter can optionally demultiplex reads by outputting barcode specific reads to
 * their associated writers..
 *
 * @param <CLUSTER_OUTPUT_RECORD> The type of record that this converter will convert to.
 */
public abstract class BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    public static final Set<IlluminaDataType> DATA_TYPES_WITH_BARCODE = new HashSet<>(Arrays.asList(
            IlluminaDataType.BaseCalls,
            IlluminaDataType.QualityScores,
            IlluminaDataType.Position,
            IlluminaDataType.PF,
            IlluminaDataType.Barcodes));

    public static final Set<IlluminaDataType> DATA_TYPES_WITHOUT_BARCODE = new HashSet<>(Arrays.asList(
            IlluminaDataType.BaseCalls,
            IlluminaDataType.QualityScores,
            IlluminaDataType.Position,
            IlluminaDataType.PF));
    protected static final Log log = Log.getInstance(UnsortedBasecallsConverter.class);

    protected final IlluminaDataProviderFactory factory;
    protected final boolean demultiplex;
    protected final boolean ignoreUnexpectedBarcodes;
    protected final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    protected final boolean includeNonPfReads;
    protected final int numThreads;
    protected final ProgressLogger readProgressLogger = new ProgressLogger(log, 1000000, "Read");
    protected final ProgressLogger writeProgressLogger = new ProgressLogger(log, 1000000, "Write");
    protected final Map<Integer, List<? extends Runnable>> completedWork = new HashMap<>();
    protected final ThreadPoolExecutorWithExceptions tileWriteExecutor;
    protected final ThreadPoolExecutorWithExceptions tileReadExecutor;
    protected final ThreadPoolExecutorWithExceptions completedWorkExecutor = new ThreadPoolExecutorWithExceptions(1);
    protected ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter = null;
    protected List<Integer> tiles;
    protected boolean tileProcessingComplete = false;

    /**
     * Constructs a new BasecallsConverter object.
     *
     * @param basecallsDir                 Where to read basecalls from.
     * @param barcodesDir                  Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lane                         What lane to process.
     * @param readStructure                How to interpret each cluster.
     * @param barcodeRecordWriterMap       Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                                     one writer stored with key=null.
     * @param demultiplex                  If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param numThreads                   Controls number of threads.
     * @param firstTile                    (For debugging) If non-null, start processing at this tile.
     * @param tileLimit                    (For debugging) If non-null, process no more than this many tiles.
     * @param bclQualityEvaluationStrategy The basecall quality evaluation strategy that is applyed to decoded base calls.
     * @param ignoreUnexpectedBarcodes     If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap.
     * @param applyEamssFiltering          If true, apply EAMSS filtering if parsing BCLs for bases and quality scores.
     * @param includeNonPfReads            If true, will include ALL reads (including those which do not have PF set).
     *                                     This option does nothing for instruments that output cbcls (Novaseqs)
     */
    public BasecallsConverter(
            final File basecallsDir,
            final File barcodesDir,
            final int lane,
            final ReadStructure readStructure,
            final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
            final boolean demultiplex,
            final int numThreads,
            final Integer firstTile,
            final Integer tileLimit,
            final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
            final boolean ignoreUnexpectedBarcodes,
            final boolean applyEamssFiltering,
            final boolean includeNonPfReads,
            final int numWriteThreads
    ) {
        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
        this.ignoreUnexpectedBarcodes = ignoreUnexpectedBarcodes;
        this.demultiplex = demultiplex;
        this.numThreads = numThreads;

        this.factory = new IlluminaDataProviderFactory(basecallsDir,
                barcodesDir, lane, readStructure, bclQualityEvaluationStrategy, getDataTypesFromReadStructure(readStructure, demultiplex));
        this.factory.setApplyEamssFiltering(applyEamssFiltering);
        this.includeNonPfReads = includeNonPfReads;
        this.tiles = factory.getAvailableTiles();
        tiles.sort(TILE_NUMBER_COMPARATOR);
        setTileLimits(firstTile, tileLimit);
        tileWriteExecutor = new ThreadPoolExecutorWithExceptions(numWriteThreads);
        tileWriteExecutor.setKeepAliveTime(500, TimeUnit.MILLISECONDS);
        tileReadExecutor = new ThreadPoolExecutorWithExceptions(numThreads);
        final CompletedWorkChecker workChecker = new CompletedWorkChecker();
        completedWorkExecutor.submit(workChecker);
        completedWorkExecutor.shutdown();
    }

    /**
     * Abstract method for processing tiles of data and outputting records for each barcode.
     *
     * @param barcodes The barcodes used optionally for demultiplexing. Must contain at least a single null value if
     *                 no demultiplexing is being done.
     */
    public abstract void processTilesAndWritePerSampleOutputs(final Set<String> barcodes);

    /**
     * Interface that defines a converter that takes ClusterData and returns OUTPUT_RECORD type objects.
     *
     * @param <OUTPUT_RECORD> The recode type to convert to.
     */
    protected interface ClusterDataConverter<OUTPUT_RECORD> {
        /**
         * Creates the OUTPUT_RECORDs from the cluster
         */
        OUTPUT_RECORD convertClusterToOutputRecord(final ClusterData cluster);
    }

    /**
     * Interface that defines a writer that will write out OUTPUT_RECORD type objects.
     *
     * @param <OUTPUT_RECORD> The recode type to convert to.
     */
    protected interface ConvertedClusterDataWriter<OUTPUT_RECORD> {
        /**
         * Write out a single record of type OUTPUT_RECORD.
         *
         * @param rec The record to write.
         */
        void write(final OUTPUT_RECORD rec);

        /**
         * Closes the writer.
         */
        void close();
    }

    protected void awaitTileProcessingCompletion() {
        tileReadExecutor.shutdown();
        // Wait for all the read threads to complete before checking for errors
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Reading executor", tileReadExecutor, Duration.ofMinutes(5));
        tileProcessingComplete = true;

        try {
            // Check for reading errors
            if (tileReadExecutor.hasError()) {
                interruptAndShutdownExecutors(tileReadExecutor, completedWorkExecutor, tileWriteExecutor);
            }

            synchronized (completedWork) {
                log.debug("Final notification of work complete.");
                completedWork.notifyAll();
            }

            // Wait for tile processing synchronization to complete
            ThreadPoolExecutorUtil.awaitThreadPoolTermination("Tile completion executor", completedWorkExecutor, Duration.ofMinutes(5));

            // Check for tile work synchronization errors
            if (completedWorkExecutor.hasError()) {
                interruptAndShutdownExecutors(tileReadExecutor, completedWorkExecutor, tileWriteExecutor);
            }
        } finally {
            // We are all done scheduling work. Now close the writers.
            barcodeRecordWriterMap.values().forEach(ConvertedClusterDataWriter::close);
        }
    }

    protected void notifyWorkComplete(int tileNum, List<? extends Runnable> pumpList) {
        synchronized (completedWork) {
            log.debug("Notifying completed work. Tile: " + tileNum);
            completedWork.put(tileNum, pumpList);
            completedWork.notifyAll();
        }
    }

    /**
     * CompletedWorkChecker is notified by the TileProcessor threads as work on a tile is complete and the
     * records are ready for writing. It also ensures that tiles are written out in the proper order according
     * by keep track of the current tile index in the sorted list of all tiles to be processed.
     * <p>
     * If a tile is finished and it is next in line to be written the CompletedWorkChecker thread will call
     * writeRecords on the SortedRecordToWriterPump.
     */
    protected class CompletedWorkChecker implements Runnable {
        private int currentTileIndex = 0;

        @Override
        public void run() {
            try {
                checkCompletedWork();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        private void checkCompletedWork() throws InterruptedException {
            synchronized (completedWork) {
                while (currentTileIndex < tiles.size()) {
                    // Wait only if tile processing is still occurring
                    if (!tileProcessingComplete) {
                        log.debug("Waiting for completed work.");
                        completedWork.wait();
                    }
                    final Integer currentTile = tiles.get(currentTileIndex);
                    if (completedWork.containsKey(currentTile)) {
                        if (tileWriteExecutor.getQueue().size() == 0
                                && tileWriteExecutor.getActiveCount() == 0
                                && tileWriteExecutor.getTaskCount() == tileWriteExecutor.getCompletedTaskCount()) {
                            // tileWriteExecutor will report 0 active workers even though the worker is still tidying up
                            // so we add a small sleep to ensure it is finished before moving on to the next tile
                            Thread.sleep(100);
                            log.debug("Writing out tile. Tile: " + currentTile);
                            completedWork.get(currentTile).forEach(tileWriteExecutor::submit);
                            currentTileIndex++;
                        }
                    }
                }
                tileWriteExecutor.shutdown();
                ThreadPoolExecutorUtil.awaitThreadPoolTermination("Tile completion executor", tileWriteExecutor, Duration.ofMinutes(5));

                // Check for writing errors
                if (tileWriteExecutor.hasError()) {
                    interruptAndShutdownExecutors(tileReadExecutor, completedWorkExecutor, tileWriteExecutor);
                }
            }
        }
    }

    /**
     * A comparator used to sort Illumina tiles in their proper order.
     * Because the tile number is followed by a colon, a tile number that is a prefix of another tile number
     * should sort after. (e.g. 10 sorts after 100). Tile numbers with the same number of digits are sorted numerically.
     */
    public static final Comparator<Integer> TILE_NUMBER_COMPARATOR = (integer1, integer2) -> {
        final String s1 = integer1.toString();
        final String s2 = integer2.toString();

        if (s1.length() < s2.length()) {
            if (s2.startsWith(s1)) {
                return 1;
            }
        } else if (s2.length() < s1.length() && s1.startsWith(s2)) {
            return -1;
        }
        return s1.compareTo(s2);
    };

    /**
     * Applies an lane and tile based regex to return all files matching that regex for each tile.
     *
     * @param baseDirectory The directory to search for tiled files.
     * @param pattern       The pattern used to match files.
     * @return A file array of all of the tile based files that match the regex pattern.
     */
    public static File[] getTiledFiles(final File baseDirectory, final Pattern pattern) {
        return IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
    }

    /**
     * Given a read structure return the data types that need to be parsed for this run
     *
     * @param readStructure The read structure that defines how the read is set up.
     * @param demultiplex   If true, output is split by barcode, otherwise all are written to the same output stream.
     * @return A data type array for each piece of data needed to satisfy the read structure.
     */
    protected static Set<IlluminaDataType> getDataTypesFromReadStructure(final ReadStructure readStructure,
                                                                         final boolean demultiplex) {
        if (!readStructure.hasSampleBarcode() || !demultiplex) {
            return DATA_TYPES_WITHOUT_BARCODE;
        } else {
            return DATA_TYPES_WITH_BARCODE;
        }
    }

    /**
     * Gets the data provider factory used to create the underlying data provider.
     *
     * @return A factory used for create the underlying data provider.
     */
    protected IlluminaDataProviderFactory getFactory() {
        return factory;
    }

    /**
     * Must be called before doTileProcessing.  This is not passed in the ctor because often the
     * IlluminaDataProviderFactory is needed in order to construct the converter.
     *
     * @param converter Converts ClusterData to CLUSTER_OUTPUT_RECORD
     */
    protected void setConverter(final ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter) {
        this.converter = converter;
    }

    /**
     * Uses the firstTile and tileLimit parameters to set which tiles will be processed. The processor will
     * start with firstTile and continue to process tiles in order until it has processed at most tileLimit tiles.
     *
     * @param firstTile The tile to begin processing at.
     * @param tileLimit The maximum number of tiles to process.
     */
    protected void setTileLimits(final Integer firstTile, final Integer tileLimit) {
        if (firstTile != null) {
            int i;
            for (i = 0; i < tiles.size(); ++i) {
                if (tiles.get(i).intValue() == firstTile.intValue()) {
                    tiles = tiles.subList(i, tiles.size());
                    break;
                }
            }
            if (tiles.get(0).intValue() != firstTile.intValue()) {
                throw new PicardException("firstTile=" + firstTile + ", but that tile was not found.");
            }
        }
        if (tileLimit != null && tiles.size() > tileLimit) {
            tiles = tiles.subList(0, tileLimit);
        }
    }

    protected void interruptAndShutdownExecutors(ThreadPoolExecutorWithExceptions... executors) {
        int tasksRunning = Arrays.stream(executors).mapToInt(test -> test.shutdownNow().size()).sum();
        throw new PicardException("Exceptions in tile processing. There were " + tasksRunning
                + " tasks were still running or queued and have been cancelled.");
    }
}
