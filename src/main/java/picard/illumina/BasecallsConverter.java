package picard.illumina;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.illumina.parser.*;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.ThreadPoolExecutorUtil;
import picard.util.ThreadPoolExecutorWithExceptions;

import java.io.File;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(BasecallsConverter.class);
    public static final IlluminaDataType[] DATA_TYPES_WITH_BARCODE = {
            IlluminaDataType.BaseCalls,
            IlluminaDataType.QualityScores,
            IlluminaDataType.Position,
            IlluminaDataType.PF,
            IlluminaDataType.Barcodes
    };
    public static final IlluminaDataType[] DATA_TYPES_WITHOUT_BARCODE =
            Arrays.copyOfRange(DATA_TYPES_WITH_BARCODE, 0, DATA_TYPES_WITH_BARCODE.length - 1);

    protected final BclQualityEvaluationStrategy bclQualityEvaluationStrategy;
    protected final IlluminaDataProviderFactory factory;
    final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator;
    final int maxReadsInRamPerTile;
    final boolean demultiplex;
    final List<File> tmpDirs;
    final boolean ignoreUnexpectedBarcodes;
    final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype;
    final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass;
    final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap;
    final ProgressLogger readProgressLogger = new ProgressLogger(log, 1000000, "Read");
    final ProgressLogger writeProgressLogger = new ProgressLogger(log, 1000000, "Write");
    final boolean includeNonPfReads;
    private final Map<String, ThreadPoolExecutorWithExceptions> barcodeWriterThreads = new HashMap<>();
    private final Map<Integer, List<RecordWriter>> completedWork = new HashMap<>();
    private final Map<Integer, File> barcodesFiles = new HashMap<>();
    protected List<Integer> tiles;
    int numThreads;
    ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter = null;
    private boolean tileProcessingComplete = false;


    /**
     * @param basecallsDir             Where to read basecalls from.
     * @param barcodesDir              Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lane                     What lane to process.
     * @param readStructure            How to interpret each cluster.
     * @param barcodeRecordWriterMap   Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                                 one writer stored with key=null.
     * @param demultiplex              If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param maxReadsInRamPerTile     Configures number of reads each tile will store in RAM before spilling to disk.
     * @param tmpDirs                  For SortingCollection spilling.
     * @param numProcessors            Controls number of threads.  If <= 0, the number of threads allocated is
     *                                 available cores - numProcessors.
     * @param firstTile                (For infoging) If non-null, start processing at this tile.
     * @param tileLimit                (For infoging) If non-null, process no more than this many tiles.
     * @param outputRecordComparator   For sorting output records within a single tile.
     * @param codecPrototype           For spilling output records to disk.
     * @param outputRecordClass        Inconveniently needed to create SortingCollections.
     * @param ignoreUnexpectedBarcodes If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap.
     * @param applyEamssFiltering      If true, apply EAMSS filtering if parsing BCLs for bases and quality scores.
     * @param includeNonPfReads        If true, will include ALL reads (including those which do not have PF set).
     *                                 This option does nothing for instruments that output cbcls (Novaseqs)
     */
    public BasecallsConverter(final File basecallsDir, final File barcodesDir, final int lane,
                              final ReadStructure readStructure,
                              final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
                              final boolean demultiplex,
                              final int maxReadsInRamPerTile,
                              final List<File> tmpDirs, final int numProcessors,
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
        this.barcodeRecordWriterMap = barcodeRecordWriterMap;
        this.maxReadsInRamPerTile = maxReadsInRamPerTile;
        this.tmpDirs = tmpDirs;
        this.codecPrototype = codecPrototype;
        this.ignoreUnexpectedBarcodes = ignoreUnexpectedBarcodes;
        this.demultiplex = demultiplex;
        this.outputRecordComparator = outputRecordComparator;
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;
        this.outputRecordClass = outputRecordClass;
        this.factory = new IlluminaDataProviderFactory(basecallsDir,
                barcodesDir, lane, readStructure, bclQualityEvaluationStrategy, getDataTypesFromReadStructure(readStructure, demultiplex));

        if (numProcessors == 0) {
            this.numThreads = Runtime.getRuntime().availableProcessors();
        } else if (numProcessors < 0) {
            this.numThreads = Runtime.getRuntime().availableProcessors() + numProcessors;
        } else {
            this.numThreads = numProcessors;
        }

        this.factory.setApplyEamssFiltering(applyEamssFiltering);
        this.includeNonPfReads = includeNonPfReads;
        this.tiles = new ArrayList<>();
        final File laneDir = new File(basecallsDir, IlluminaFileUtil.longLaneStr(lane));

        final Pattern filterRegex = Pattern.compile(ParameterizedFileUtil.escapePeriods(
                ParameterizedFileUtil.makeLaneTileRegex(".filter", lane)));
        File[] filterFiles = getTiledFiles(laneDir, filterRegex);
        for (final File filterFile : filterFiles) {
            final Matcher tileMatcher = filterRegex.matcher(filterFile.getName());
            if (tileMatcher.matches()) {
                tiles.add(Integer.valueOf(tileMatcher.group(1)));
            }
        }

        IOUtil.assertFilesAreReadable(Arrays.asList(filterFiles));
        tiles.sort(TILE_NUMBER_COMPARATOR);
        setTileLimits(firstTile, tileLimit);
        barcodeRecordWriterMap.keySet().forEach(barcode -> barcodeWriterThreads.put(barcode, new ThreadPoolExecutorWithExceptions(1)));

        if (demultiplex) {
            final Pattern barcodeRegex = Pattern.compile(ParameterizedFileUtil.escapePeriods(
                    ParameterizedFileUtil.makeBarcodeRegex(lane)));
            final File[] barcodeTileFiles = getTiledFiles(barcodesDir, barcodeRegex);
            if (barcodeTileFiles.length < tiles.size()) {
                throw new PicardException(String.format(
                        "Barcode files are required for each tile. Found %d expected %d.",
                        barcodeTileFiles.length, tiles.size()));
            }
            for (final File barcodeFile : barcodeTileFiles) {
                final Matcher tileMatcher = barcodeRegex.matcher(barcodeFile.getName());
                if (tileMatcher.matches()) {
                    barcodesFiles.put(Integer.valueOf(tileMatcher.group(1)), barcodeFile);
                }
            }
        }
    }

    public static File[] getTiledFiles(final File baseDirectory, final Pattern pattern) {
        return IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
    }

    /**
     * Given a read structure return the data types that need to be parsed for this run
     */
    protected static IlluminaDataType[] getDataTypesFromReadStructure(final ReadStructure readStructure,
                                                                      final boolean demultiplex) {
        if (!readStructure.hasSampleBarcode() || !demultiplex) {
            return DATA_TYPES_WITHOUT_BARCODE;
        } else {
            return DATA_TYPES_WITH_BARCODE;
        }
    }

    public IlluminaDataProviderFactory getFactory() {
        return factory;
    }

    /**
     * Must be called before doTileProcessing.  This is not passed in the ctor because often the
     * IlluminaDataProviderFactory is needed in order to construct the converter.
     *
     * @param converter Converts ClusterData to CLUSTER_OUTPUT_RECORD
     */
    void setConverter(final ClusterDataConverter<CLUSTER_OUTPUT_RECORD> converter) {
        this.converter = converter;
    }

    void setTileLimits(final Integer firstTile, final Integer tileLimit) {
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

    public void doTileProcessing() {

        final ThreadPoolExecutorWithExceptions completedWorkExecutor = new ThreadPoolExecutorWithExceptions(1);

        final CompletedWorkChecker workChecker = new CompletedWorkChecker();
        completedWorkExecutor.submit(workChecker);
        completedWorkExecutor.shutdown();

        //  Thread by surface tile
        final ThreadPoolExecutorWithExceptions tileProcessingExecutor = new ThreadPoolExecutorWithExceptions(numThreads);

        for (final Integer tile : tiles) {
            tileProcessingExecutor.submit(new TileProcessor(tile, barcodesFiles.get(tile)));
        }

        tileProcessingExecutor.shutdown();

        // Wait for all the threads to complete before checking for errors
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Reading executor", tileProcessingExecutor, Duration.ofSeconds(5));
        tileProcessingComplete = true;
        synchronized (completedWork) {
            log.debug("Final notification of work complete.");
            completedWork.notifyAll();
        }
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Tile completion executor", completedWorkExecutor, Duration.ofMinutes(5));

        // We are all done scheduling work. Now schedule the closers.
        barcodeRecordWriterMap.forEach((barcode, writer) -> barcodeWriterThreads.get(barcode).submit(new Closer(writer, barcode)));

        barcodeWriterThreads.values().forEach(ThreadPoolExecutor::shutdown);
        barcodeWriterThreads.forEach((barcode, executor) -> ThreadPoolExecutorUtil.awaitThreadPoolTermination(barcode + " writer", executor, Duration.ofMinutes(5)));


        if (tileProcessingExecutor.hasError() ||
                completedWorkExecutor.hasError() ||
                barcodeWriterThreads.values().stream().anyMatch(ThreadPoolExecutorWithExceptions::hasError)) {
            int tasksStillRunning = completedWorkExecutor.shutdownNow().size();
            tasksStillRunning += barcodeWriterThreads.values().stream().mapToLong(executor -> executor.shutdownNow().size()).sum();
            throw new PicardException("Exceptions in tile processing. There were " + tasksStillRunning
                    + " tasks were still running or queued and have been cancelled.");
        }
    }


    interface ClusterDataConverter<OUTPUT_RECORD> {
        /**
         * Creates the OUTPUT_RECORDs from the cluster
         */
        OUTPUT_RECORD convertClusterToOutputRecord(final ClusterData cluster);
    }

    interface ConvertedClusterDataWriter<OUTPUT_RECORD> {
        void write(final OUTPUT_RECORD rec);
        void close();
    }

    private class RecordWriter implements Runnable {
        private final SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection;
        private final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer;
        private final String barcode;

        RecordWriter(final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer,
                     final SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection, final String barcode) {
            this.writer = writer;
            this.recordCollection = recordCollection;
            this.barcode = barcode;
        }

        @Override
        public void run() {
            log.info("Writing out barcode " + barcode);
            for (final CLUSTER_OUTPUT_RECORD record : recordCollection) {
                writer.write(record);
                writeProgressLogger.record(null, 0);
            }
        }

        public String getBarcode() {
            return barcode;
        }
    }

    private class Closer implements Runnable {
        private final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer;
        private final String barcode;

        private Closer(final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer, final String barcode) {
            this.writer = writer;
            this.barcode = barcode;
        }

        @Override
        public void run() {
            log.info("Closing writer for barcode " + barcode);
            this.writer.close();
        }
    }

    private class TileProcessor implements Runnable {
        private final int tileNum;
        private final Map<String, SortingCollection<CLUSTER_OUTPUT_RECORD>> barcodeToRecordCollection = new HashMap<>();
        private final File barcodeFile;

        TileProcessor(final int tileNum, final File barcodeFile) {
            this.tileNum = tileNum;
            this.barcodeFile = barcodeFile;
        }

        @Override
        public void run() {
            final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(tileNum, barcodeFile);

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                final String barcode = (demultiplex ? cluster.getMatchedBarcode() : null);
                if (cluster.isPf() || includeNonPfReads) {
                    addRecord(barcode, converter.convertClusterToOutputRecord(cluster));
                }
            }

            dataProvider.close();

            final List<RecordWriter> writerList = new ArrayList<>();
            barcodeToRecordCollection.forEach((barcode, value) -> {
                value.doneAdding();
                final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer = barcodeRecordWriterMap.get(barcode);
                writerList.add(new RecordWriter(writer, value, barcode));
            });

            notifyWorkComplete(writerList);

            log.info("Finished processing tile " + tileNum);
        }

        private void notifyWorkComplete(List<RecordWriter> writerList) {
            synchronized (completedWork) {
                log.debug("Notifying completed work. Tile: " + tileNum);
                completedWork.put(tileNum, writerList);
                completedWork.notifyAll();
            }
        }

        private synchronized void addRecord(final String barcode, final CLUSTER_OUTPUT_RECORD record) {
            // Grab the existing collection, or initialize it if it doesn't yet exist
            SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection = this.barcodeToRecordCollection.get(barcode);
            if (recordCollection == null) {
                // TODO: The implementation here for supporting ignoreUnexpectedBarcodes is not efficient,
                // but the alternative is an extensive rewrite.  We are living with the inefficiency for
                // this special case for the time being.
                if (!barcodeRecordWriterMap.containsKey(barcode)) {
                    if (ignoreUnexpectedBarcodes) {
                        return;
                    }
                    throw new PicardException(String.format("Read records with barcode %s, but this barcode was not expected.  (Is it referenced in the parameters file?)", barcode));
                }
                recordCollection = newSortingCollection();
                this.barcodeToRecordCollection.put(barcode, recordCollection);
            }
            recordCollection.add(record);
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
                    // Wait only if tile processing is still occurring
                    if (!tileProcessingComplete) {
                        log.debug("Waiting for completed work.");
                        completedWork.wait();
                    }
                    final Integer currentTile = tiles.get(currentTileIndex);
                    if (completedWork.containsKey(currentTile)) {
                        log.debug("Writing out tile. Tile: " + currentTile);
                        completedWork.get(currentTile).forEach(writer -> barcodeWriterThreads.get(writer.getBarcode()).submit(writer));
                        currentTileIndex++;
                    }
                }
            }
        }
    }

    public static final Comparator<Integer> TILE_NUMBER_COMPARATOR = (integer1, integer2) -> {
        final String s1 = integer1.toString();
        final String s2 = integer2.toString();
        // Because a the tile number is followed by a colon, a tile number that
        // is a prefix of another tile number should sort after. (e.g. 10 sorts after 100).
        if (s1.length() < s2.length()) {
            if (s2.startsWith(s1)) {
                return 1;
            }
        } else if (s2.length() < s1.length() && s1.startsWith(s2)) {
            return -1;
        }
        return s1.compareTo(s2);
    };
}
