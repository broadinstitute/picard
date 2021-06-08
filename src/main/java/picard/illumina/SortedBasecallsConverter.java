package picard.illumina;

import htsjdk.io.AsyncWriterPool;
import htsjdk.io.Writer;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.ThreadPoolExecutorUtil;
import picard.util.ThreadPoolExecutorWithExceptions;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

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
    protected static final Log log = Log.getInstance(SortedBasecallsConverter.class);
    private final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator;
    private final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype;
    private final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass;
    private final int maxReadsInRamPerTile;
    private final List<File> tmpDirs;
    private final Map<Integer, List<? extends Runnable>> completedWork = new ConcurrentHashMap<>();
    private final ThreadPoolExecutorWithExceptions tileWriteExecutor;
    private final ThreadPoolExecutorWithExceptions tileReadExecutor;
    private final ProgressLogger readProgressLogger = new ProgressLogger(log, 1000000, "Read");
    private final ProgressLogger writeProgressLogger = new ProgressLogger(log, 1000000, "Write");
    private final AtomicInteger tileWriteJobs = new AtomicInteger(0);

    /**
     * Constructs a new SortedBaseCallsConverter.
     *
     * @param basecallsDir                 Where to read basecalls from.
     * @param barcodesDir                  Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lanes                        What lanes to process.
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
            final int[] lanes,
            final ReadStructure readStructure,
            final Map<String, ? extends Writer<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
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
            final boolean includeNonPfReads,
            final AsyncWriterPool writerPool
    ) {
        super(basecallsDir, barcodesDir, lanes, readStructure, barcodeRecordWriterMap, demultiplex,
                firstTile, tileLimit, bclQualityEvaluationStrategy,
                ignoreUnexpectedBarcodes, applyEamssFiltering, includeNonPfReads, writerPool);

        this.tmpDirs = tmpDirs;
        this.maxReadsInRamPerTile = maxReadsInRamPerTile;
        this.codecPrototype = codecPrototype;
        this.outputRecordComparator = outputRecordComparator;
        this.outputRecordClass = outputRecordClass;
        tileWriteExecutor = new ThreadPoolExecutorWithExceptions(numThreads);
        tileReadExecutor = new ThreadPoolExecutorWithExceptions(numThreads);
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
    public void processTilesAndWritePerSampleOutputs(final Set<String> barcodes) throws IOException {
        for (final Integer tile : tiles) {
            tileReadExecutor.submit(new TileProcessor(tile, barcodes));
        }
        awaitTileProcessingCompletion();
    }

    /**
     * SortedRecordToWriterPump takes a collection of output records and writes them using a
     * ConvertedClusterDataWriter.
     */
    private class SortedRecordToWriterPump implements Runnable {
        private final SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection;
        private final Writer<CLUSTER_OUTPUT_RECORD> writer;

        SortedRecordToWriterPump(final Writer<CLUSTER_OUTPUT_RECORD> writer,
                                 final SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection) {
            this.writer = writer;
            this.recordCollection = recordCollection;
        }

        @Override
        public void run() {
            tileWriteJobs.incrementAndGet();
            for (final CLUSTER_OUTPUT_RECORD record : recordCollection) {
                writer.write(record);
                writeProgressLogger.record(null, 0);
            }
            recordCollection.cleanup();
            int writeJobsRemaining = tileWriteJobs.decrementAndGet();

            if (writeJobsRemaining == 0) {
                synchronized (tileWriteJobs) {
                    tileWriteJobs.notifyAll();
                }
            }
        }
    }

    /**
     * TileProcessor is a Runnable that process all records for a given tile. It uses the underlying
     * IlluminaDataProvider to iterate over cluster data for a specific tile. Records are added to a
     * SortingCollection as they are read and decoded. This processor also optionally filters non-PF reads.
     * In addition, it will optionally demultiplex by barcode.
     * <p>
     * After the tile processing is complete it notifies the CompletedWorkChecker that data is ready
     * for writing.
     */
    private class TileProcessor implements Runnable {
        private final int tileNum;
        private final Map<String, SortingCollection<CLUSTER_OUTPUT_RECORD>> barcodeToRecordCollection;

        TileProcessor(final int tileNum, final Set<String> barcodes) {
            this.tileNum = tileNum;
            this.barcodeToRecordCollection = new HashMap<>(barcodes.size(), 1.0f);
            for (String barcode : barcodes) {
                SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection = createSortingCollection();
                this.barcodeToRecordCollection.put(barcode, recordCollection);
            }
        }

        @Override
        public void run() {
            final List<SortedRecordToWriterPump> writerList = new ArrayList<>();
            for (IlluminaDataProviderFactory laneFactory : laneFactories) {
                if (laneFactory.getAvailableTiles().contains(tileNum)) {
                    final BaseIlluminaDataProvider dataProvider = laneFactory.makeDataProvider(tileNum);

                    while (dataProvider.hasNext()) {
                        final ClusterData cluster = dataProvider.next();
                        readProgressLogger.record(null, 0);
                        if (includeNonPfReads || cluster.isPf()) {
                            addRecord(cluster.getMatchedBarcode(), converter.convertClusterToOutputRecord(cluster));
                        }
                    }
                    dataProvider.close();
                }
            }
            barcodeToRecordCollection.forEach((barcode, value) -> {
                value.doneAdding();
                final Writer<CLUSTER_OUTPUT_RECORD> writer = barcodeRecordWriterMap.get(barcode);
                log.debug("Writing out barcode " + barcode);
                writerList.add(new SortedRecordToWriterPump(writer, value));
            });

            completedWork.put(tileNum, writerList);

            log.debug("Finished processing tile " + tileNum);
        }

        private synchronized void addRecord(final String barcode, final CLUSTER_OUTPUT_RECORD record) {
            SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection = this.barcodeToRecordCollection.get(barcode);

            if (recordCollection != null) {
                recordCollection.add(record);
            } else if (!ignoreUnexpectedBarcodes) {
                throw new PicardException(String.format("Read records with barcode %s, but this barcode was not expected.  (Is it referenced in the parameters file?)", barcode));
            }
        }

        private synchronized SortingCollection<CLUSTER_OUTPUT_RECORD> createSortingCollection() {
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

    protected void awaitTileProcessingCompletion() throws IOException {
        tileReadExecutor.shutdown();
        // Wait for all the read threads to complete before checking for errors
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Reading executor", tileReadExecutor, Duration.ofMinutes(5));

        // Check for reading errors
        if (tileReadExecutor.hasError()) {
            interruptAndShutdownExecutors(tileReadExecutor, tileWriteExecutor);
        }

        int tileProcessingIndex = 0;

        while (tileProcessingIndex < tiles.size()) {
            if (tileWriteJobs.get() == 0) {
                completedWork.get(tiles.get(tileProcessingIndex)).forEach(tileWriteExecutor::submit);
                tileProcessingIndex++;
                try {
                    synchronized (tileWriteJobs) {
                        tileWriteJobs.wait();
                        // Short sleep to ensure data is flushed.
                        Thread.sleep(500);
                    }
                } catch (InterruptedException e) {
                    throw new PicardException("Error waiting for thread lock during tile processing.", e);
                }
            }
        }

        tileWriteExecutor.shutdown();
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Writing executor", tileWriteExecutor, Duration.ofMinutes(5));

        // Check for tile work synchronization errors
        if (tileWriteExecutor.hasError()) {
            interruptAndShutdownExecutors(tileWriteExecutor);
        }
        closeWriters();
    }
}
