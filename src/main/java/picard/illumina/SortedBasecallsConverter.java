package picard.illumina;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
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

/**
 * SortedBasecallsConverter utilizes an underlying IlluminaDataProvider to convert parsed and decoded sequencing data
 * from standard Illumina formats to specific output records (FASTA records/SAM records). This data is processed
 * on a tile by tile basis and sorted based on a output record comparator.
 * <p>
 * The underlying IlluminaDataProvider apply several optional transformations that can include EAMSS filtering,
 * non-PF read filtering and quality score recoding using a BclQualityEvaluationStrategy.
 * <p>
 * The converter can also limit the scope of data that is converted from the data provider by setting the
 * tile to start on (firstTile) and the total number of tiles to process (tileLimit).
 * <p>
 * Additionally, BasecallsConverter can optionally demultiplex reads by outputting barcode specific reads to
 * their associated writers.
 */
public class SortedBasecallsConverter<CLUSTER_OUTPUT_RECORD> extends BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(SortedBasecallsConverter.class);

    private final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator;
    private final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype;
    private final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass;
    private final int maxReadsInRamPerTile;
    private final List<File> tmpDirs;
    private final ProgressLogger readProgressLogger = new ProgressLogger(log, 1000000, "Read");
    private final ProgressLogger writeProgressLogger = new ProgressLogger(log, 1000000, "Write");
    private final Map<Integer, List<SortedRecordToWriterPump>> completedWork = new HashMap<>();
    private boolean tileProcessingComplete = false;
    private boolean tileProcessingError = false;

    /**
     * Constructs a new SortedBaseCallsConverter.
     *
     * @param basecallsDir                 Where to read basecalls from.
     * @param barcodesDir                  Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lane                         What lane to process.
     * @param readStructure                How to interpret each cluster.
     * @param barcodeRecordWriterMap       Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                                     one writer stored with key=null.
     * @param demultiplex                  If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param maxReadsInRamPerTile         Configures number of reads each tile will store in RAM before spilling to disk.
     * @param tmpDirs                      For SortingCollection spilling.
     * @param numThreads                   Controls number of threads.
     * @param firstTile                    (For debugging) If non-null, start processing at this tile.
     * @param tileLimit                    (For debugging) If non-null, process no more than this many tiles.
     * @param outputRecordComparator       For sorting output records within a single tile.
     * @param codecPrototype               For spilling output records to disk.
     * @param outputRecordClass            Class needed to create SortingCollections.
     * @param bclQualityEvaluationStrategy The basecall quality evaluation strategy that is applyed to decoded base calls.
     * @param ignoreUnexpectedBarcodes     If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap.
     * @param applyEamssFiltering          If true, apply EAMSS filtering if parsing BCLs for bases and quality scores.
     * @param includeNonPfReads            If true, will include ALL reads (including those which do not have PF set).
     *                                     This option does nothing for instruments that output cbcls (Novaseqs)
     */
    protected SortedBasecallsConverter(
            final File basecallsDir,
            final File barcodesDir,
            final int lane,
            final ReadStructure readStructure,
            final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
            final boolean demultiplex,
            final int maxReadsInRamPerTile,
            final List<File> tmpDirs,
            final int numThreads,
            final Integer firstTile,
            final Integer tileLimit,
            final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator,
            final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype,
            final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass,
            final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
            final boolean ignoreUnexpectedBarcodes,
            final boolean applyEamssFiltering,
            final boolean includeNonPfReads
    ) {
        super(basecallsDir, barcodesDir, lane, readStructure, barcodeRecordWriterMap, demultiplex,
                numThreads, firstTile, tileLimit, bclQualityEvaluationStrategy,
                ignoreUnexpectedBarcodes, applyEamssFiltering, includeNonPfReads);

        this.tmpDirs = tmpDirs;
        this.maxReadsInRamPerTile = maxReadsInRamPerTile;
        this.codecPrototype = codecPrototype;
        this.outputRecordComparator = outputRecordComparator;
        this.outputRecordClass = outputRecordClass;
    }

    /**
     * Set up tile processing and record writing threads for this converter.  This creates a tile processing thread
     * pool of size `numThreads`. The tile processing threads notify the completed work checking thread when they are
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

        //  Thread by surface tile
        final ThreadPoolExecutorWithExceptions tileProcessingExecutor = new ThreadPoolExecutorWithExceptions(numThreads);
        for (final Integer tile : tiles) {
            tileProcessingExecutor.submit(new TileProcessor(tile, barcodes));
        }
        tileProcessingExecutor.shutdown();

        // Wait for all the threads to complete before checking for errors
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Reading executor", tileProcessingExecutor, Duration.ofMinutes(5));
        tileProcessingComplete = true;
        tileProcessingError = tileProcessingExecutor.hasError();
        synchronized (completedWork) {
            log.debug("Final notification of work complete.");
            completedWork.notifyAll();
        }
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Tile completion executor", completedWorkExecutor, Duration.ofMinutes(5));

        // We are all done scheduling work. Now schedule the closers.
        barcodeRecordWriterMap.values().forEach(ConvertedClusterDataWriter::close);

        if (tileProcessingError ||
                completedWorkExecutor.hasError()) {
            int tasksStillRunning = completedWorkExecutor.shutdownNow().size();
            throw new PicardException("Exceptions in tile processing. There were " + tasksStillRunning
                    + " tasks were still running or queued and have been cancelled.");
        }
    }

    /**
     * SortedRecordToWriterPump takes a collection of output records and writes them using a
     * ConvertedClusterDataWriter.
     */
    private class SortedRecordToWriterPump {
        private final SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection;
        private final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer;

        SortedRecordToWriterPump(final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer,
                                 final SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection) {
            this.writer = writer;
            this.recordCollection = recordCollection;
        }

        public void writeRecords() {
            for (final CLUSTER_OUTPUT_RECORD record : recordCollection) {
                writer.write(record);
                writeProgressLogger.record(null, 0);
            }
        }
    }

    /**
     * TileProcessor is a Runnable that process all records for a given tile. It uses the underlying
     * IlluminaDataProvider to iterate over cluster data for a specific tile. Records are added to a
     * SortingCollection as they are read and decoded. This processor also optionally filters non-PF reads.
     * In addition, it will optionally demultiplex by barcode.
     * <p>
     * After the tile processing is complete it notifies the SortedRecordToWriterPump that data is ready
     * for writing.
     */
    private class TileProcessor implements Runnable {
        private final int tileNum;
        private final Map<String, SortingCollection<CLUSTER_OUTPUT_RECORD>> barcodeToRecordCollection;

        TileProcessor(final int tileNum, final Set<String> barcodes) {
            this.tileNum = tileNum;
            this.barcodeToRecordCollection = new HashMap<>(barcodes.size(), 1.0f);
            for (String barcode : barcodes) {
                SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection = newSortingCollection();
                this.barcodeToRecordCollection.put(barcode, recordCollection);
            }
        }

        @Override
        public void run() {
            final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(tileNum);

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                if (cluster.isPf() || includeNonPfReads) {
                    final String barcode = (demultiplex ? cluster.getMatchedBarcode() : null);
                    addRecord(barcode, converter.convertClusterToOutputRecord(cluster));
                }
            }

            dataProvider.close();

            final List<SortedRecordToWriterPump> writerList = new ArrayList<>();
            barcodeToRecordCollection.forEach((barcode, value) -> {
                value.doneAdding();
                final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer = barcodeRecordWriterMap.get(barcode);
                log.debug("Writing out barcode " + barcode);
                writerList.add(new SortedRecordToWriterPump(writer, value));
            });

            notifyWorkComplete(writerList);

            log.debug("Finished processing tile " + tileNum);
        }

        private void notifyWorkComplete(List<SortedRecordToWriterPump> writerList) {
            synchronized (completedWork) {
                log.debug("Notifying completed work. Tile: " + tileNum);
                completedWork.put(tileNum, writerList);
                completedWork.notifyAll();
            }
        }

        private synchronized void addRecord(final String barcode, final CLUSTER_OUTPUT_RECORD record) {
            SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection = this.barcodeToRecordCollection.get(barcode);

            if (recordCollection != null) {
                recordCollection.add(record);
            } else if (!ignoreUnexpectedBarcodes) {
                throw new PicardException(String.format("Read records with barcode %s, but this barcode was not expected.  (Is it referenced in the parameters file?)", barcode));
            }
        }

        private synchronized SortingCollection<CLUSTER_OUTPUT_RECORD> newSortingCollection() {
            final int maxRecordsInRam =
                    Math.max(1, maxReadsInRamPerTile /
                            barcodeRecordWriterMap.size());
            return SortingCollection.newInstanceFromPaths(
                    outputRecordClass,
                    codecPrototype.clone(),
                    outputRecordComparator,
                    maxRecordsInRam,
                    IOUtil.filesToPaths(tmpDirs));
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
            synchronized (completedWork) {
                while (currentTileIndex < tiles.size()) {
                    if (tileProcessingError) {
                        throw new InterruptedException("Tile processing error, shutting down completed work checker.");
                    }
                    // Wait only if tile processing is still occurring
                    if (!tileProcessingComplete) {
                        log.debug("Waiting for completed work.");
                        completedWork.wait();
                    }
                    final Integer currentTile = tiles.get(currentTileIndex);
                    if (completedWork.containsKey(currentTile)) {
                        log.debug("Writing out tile. Tile: " + currentTile);
                        completedWork.get(currentTile).forEach(SortedRecordToWriterPump::writeRecords);
                        currentTileIndex++;
                    }
                }
            }
        }
    }
}
