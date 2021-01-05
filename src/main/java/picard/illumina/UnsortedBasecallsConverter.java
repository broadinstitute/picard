package picard.illumina;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.PicardException;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.ThreadPoolExecutorUtil;
import picard.util.ThreadPoolExecutorWithExceptions;

import java.io.File;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * UnortedBasecallsConverter utilizes an underlying IlluminaDataProvider to convert parsed and decoded sequencing data
 * from standard Illumina formats to specific output records (FASTA records/SAM records). This data is processed
 * on a tile by tile basis.
 * <p>
 * The underlying IlluminaDataProvider applies several optional transformations that can include EAMSS filtering,
 * non-PF read filtering and quality score recoding using a BclQualityEvaluationStrategy.
 * <p>
 * The converter can also limit the scope of data that is converted from the data provider by setting the
 * tile to start on (firstTile) and the total number of tiles to process (tileLimit).
 * <p>
 * Additionally, BasecallsConverter can optionally demultiplex reads by outputting barcode specific reads to
 * their associated writers.
 */
public class UnsortedBasecallsConverter<CLUSTER_OUTPUT_RECORD> extends BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(UnsortedBasecallsConverter.class);
    private final ProgressLogger readProgressLogger = new ProgressLogger(log, 1000000, "Read");
    private final ProgressLogger writeProgressLogger = new ProgressLogger(log, 1000000, "Write");
    private final Map<Integer, Queue<ClusterData>> tileReadCache = new HashMap<>();
    private boolean tileProcessingComplete = false;
    private boolean tileProcessingError = false;
    private int tilesProcessing = 0;

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
     * @param numThreads                   Controls number of threads.  If <= 0, the number of threads allocated is
     *                                     available cores - numProcessors.
     * @param firstTile                    (For debugging) If non-null, start processing at this tile.
     * @param tileLimit                    (For debugging) If non-null, process no more than this many tiles.
     * @param bclQualityEvaluationStrategy The basecall quality evaluation strategy that is applyed to decoded base calls.
     * @param ignoreUnexpectedBarcodes     If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap.
     * @param applyEamssFiltering          If true, apply EAMSS filtering if parsing BCLs for bases and quality scores.
     * @param includeNonPfReads            If true, will include ALL reads (including those which do not have PF set).
     */
    protected UnsortedBasecallsConverter(
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
            final boolean includeNonPfReads
    ) {
        super(basecallsDir, barcodesDir, lane, readStructure, barcodeRecordWriterMap, demultiplex,
                numThreads, firstTile, tileLimit, bclQualityEvaluationStrategy,
                ignoreUnexpectedBarcodes, applyEamssFiltering, includeNonPfReads);
    }

    /**
     * Set up tile processing and record writing threads for this converter.  This creates a tile reading thread
     * pool of size 4. The tile processing threads notify the completed work checking thread when they are
     * done processing a thread. The completed work checking thread will then dispatch the record writing for tiles
     * in order.
     *
     * @param barcodes The barcodes used for demultiplexing. When there is no demultiplexing done this should be a Set
     *                 containing a single null value.
     */
    @Override
    public void processTilesAndWritePerSampleOutputs(final Set<String> barcodes) {
        final ThreadPoolExecutorWithExceptions completedWorkExecutor = new ThreadPoolExecutorWithExceptions(1);
        final CompletedWorkChecker workChecker = new CompletedWorkChecker();
        completedWorkExecutor.submit(workChecker);
        completedWorkExecutor.shutdown();

        final ThreadPoolExecutorWithExceptions tileReadExecutor = new ThreadPoolExecutorWithExceptions(numThreads);
        int MAX_TILES_IN_CACHE = 4;

        int tilesSubmitted = 0;
        while( tilesSubmitted < tiles.size()){
            if(tilesProcessing < MAX_TILES_IN_CACHE) {
                int tile = tiles.get(tilesSubmitted);
                tileReadExecutor.submit(new TileReadProcessor(tile));
                tilesProcessing++;
                tilesSubmitted++;
            } else {
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    throw new PicardException("Tile processing thread interrupted: " + e.getMessage());
                }
            }
        }
        tileReadExecutor.shutdown();

        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Reading executor", tileReadExecutor, Duration.ofMinutes(5));
        tileProcessingComplete = true;
        tileProcessingError = tileReadExecutor.hasError();

        synchronized (tileReadCache) {
            log.debug("Final notification of work complete.");
            tileReadCache.notifyAll();
        }
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Tile completion executor", completedWorkExecutor, Duration.ofMinutes(5));

        if (tileProcessingError ||
                completedWorkExecutor.hasError()) {
            int tasksStillRunning = completedWorkExecutor.shutdownNow().size();
            throw new PicardException("Exceptions in tile processing. There were " + tasksStillRunning
                    + " tasks were still running or queued and have been cancelled.");
        }

        barcodeRecordWriterMap.values().forEach(ConvertedClusterDataWriter::close);
    }

    /**
     * TileProcessor is a Runnable that process all records for a given tile. It uses the underlying
     * IlluminaDataProvider to iterate over cluster data for a specific tile. Records are added to a
     * cache as they are read.
     * <p>
     * After the tile processing is complete it notifies the CompletedWorkChecker that data is ready
     * for writing.
     */
    private class TileReadProcessor implements Runnable {
        private final int tileNum;

        TileReadProcessor(final int tileNum) {
            this.tileNum = tileNum;
        }
        @Override
        public void run() {
            Queue<ClusterData> queue = new ArrayDeque<>();
            final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(tileNum);

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                queue.add(cluster);
            }

            dataProvider.close();
            synchronized (tileReadCache) {
                log.debug("Notifying completed work. Tile: " + tileNum);
                tileReadCache.put(tileNum, queue);
                tileReadCache.notifyAll();
            }

        }
    }

    /**
     * CompletedWorkChecker is notified by the TileProcessor threads as work on a tile is complete and the
     * records are ready for writing. It also ensures that tiles are written out in the proper order according
     * by keep track of the current tile index in the sorted list of all tiles to be processed.
     * <p>
     * If a tile is finished and it is next in line to be written the CompletedWorkChecker thread will submit
     * work to the tileWriteExecutor which uses RecordToWriterPump to write records.
     */
    private class CompletedWorkChecker implements Runnable {
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
            synchronized (tileReadCache) {
                while (currentTileIndex < tiles.size()) {
                    if (tileProcessingError) {
                        throw new InterruptedException("Tile processing error, shutting down completed work checker.");
                    }
                    // Wait only if tile processing is still occurring
                    if (!tileProcessingComplete) {
                        log.debug("Waiting for completed work.");
                        tileReadCache.wait();
                    }
                    final Integer currentTile = tiles.get(currentTileIndex);
                    if (tileReadCache.containsKey(currentTile)) {
                        log.debug("Writing out tile. Tile: " + currentTile);
                        Queue<ClusterData> clusterData = tileReadCache.get(currentTile);
                        ClusterData cluster;
                        while( (cluster = clusterData.poll()) != null){
                            if (includeNonPfReads || cluster.isPf()) {
                                barcodeRecordWriterMap.get(cluster.getMatchedBarcode()).write(converter.convertClusterToOutputRecord(cluster));
                                writeProgressLogger.record(null, 0);
                            }
                        }
                        tileReadCache.remove(currentTile);
                        currentTileIndex++;
                        tilesProcessing--;
                    }
                }
            }
        }
    }
}
