package picard.illumina;


import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaFileUtil;
import picard.illumina.parser.ParameterizedFileUtil;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.illumina.parser.readers.LocsFileReader;
import picard.util.ThreadPoolExecutorUtil;
import picard.util.ThreadPoolExecutorWithExceptions;

import java.io.File;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NewIlluminaBasecallsConverter<CLUSTER_OUTPUT_RECORD> extends BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(NewIlluminaBasecallsConverter.class);
    private final List<File> cbcls;
    private final List<AbstractIlluminaPositionFileReader.PositionInfo> locs = new ArrayList<>();
    private final File[] filterFiles;
    private final Map<String, ThreadPoolExecutorWithExceptions> barcodeWriterThreads = new HashMap<>();
    private final Map<Integer, List<RecordWriter>> completedWork = Collections.synchronizedMap(new HashMap<>());
    private final Map<Integer, File> barcodesFiles = new HashMap<>();

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
     * @param firstTile                (For debugging) If non-null, start processing at this tile.
     * @param tileLimit                (For debugging) If non-null, process no more than this many tiles.
     * @param outputRecordComparator   For sorting output records within a single tile.
     * @param codecPrototype           For spilling output records to disk.
     * @param outputRecordClass        Inconveniently needed to create SortingCollections.
     * @param ignoreUnexpectedBarcodes If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap,
     */
    public NewIlluminaBasecallsConverter(final File basecallsDir, final File barcodesDir, final int lane,
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
                                         final boolean ignoreUnexpectedBarcodes) {

        super(barcodeRecordWriterMap, maxReadsInRamPerTile, tmpDirs, codecPrototype, ignoreUnexpectedBarcodes,
                demultiplex, outputRecordComparator, bclQualityEvaluationStrategy,
                outputRecordClass, numProcessors, new IlluminaDataProviderFactory(basecallsDir,
                        barcodesDir, lane, readStructure, bclQualityEvaluationStrategy));
        this.tiles = new ArrayList<>();

        barcodeRecordWriterMap.keySet().forEach(barcode -> barcodeWriterThreads.put(barcode, new ThreadPoolExecutorWithExceptions(1)));

        final File laneDir = new File(basecallsDir, IlluminaFileUtil.longLaneStr(lane));

        final File[] cycleDirs = IOUtil.getFilesMatchingRegexp(laneDir, IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN);

        //CBCLs
        cbcls = new ArrayList<>();
        Arrays.asList(cycleDirs)
                .forEach(cycleDir -> cbcls.addAll(
                        Arrays.asList(IOUtil.getFilesMatchingRegexp(
                                cycleDir, "^" + IlluminaFileUtil.longLaneStr(lane) + "_(\\d{1,5}).cbcl$"))));

        if (cbcls.size() == 0) {
            throw new PicardException("No CBCL files found.");
        }

        IOUtil.assertFilesAreReadable(cbcls);

        //locs
        final File locsFile = new File(basecallsDir.getParentFile(), AbstractIlluminaPositionFileReader.S_LOCS_FILE);
        try (LocsFileReader locsFileReader = new LocsFileReader(locsFile)) {
            while (locsFileReader.hasNext()) {
                locs.add(locsFileReader.next());
            }
        }
        IOUtil.assertFileIsReadable(locsFile);
        //filter

        final Pattern filterRegex = Pattern.compile(ParameterizedFileUtil.escapePeriods(
                ParameterizedFileUtil.makeLaneTileRegex(".filter", lane)));
        filterFiles = getTiledFiles(laneDir, filterRegex);
        for (final File filterFile : filterFiles) {
            final Matcher tileMatcher = filterRegex.matcher(filterFile.getName());
            if (tileMatcher.matches()) {
                tiles.add(Integer.valueOf(tileMatcher.group(1)));
            }
        }
        IOUtil.assertFilesAreReadable(Arrays.asList(filterFiles));
        tiles.sort(TILE_NUMBER_COMPARATOR);

        if (demultiplex) {
            final Pattern barcodeRegex = Pattern.compile(ParameterizedFileUtil.escapePeriods(
                    ParameterizedFileUtil.makeBarcodeRegex(lane)));
            final File[] barcodeTileFiles = getTiledFiles(barcodesDir, barcodeRegex);
            if (barcodeTileFiles.length != tiles.size()) {
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

        setTileLimits(firstTile, tileLimit);
    }

    public static File[] getTiledFiles(final File baseDirectory, final Pattern pattern) {
        return IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
    }

    @Override
    public void doTileProcessing() {

        final ThreadPoolExecutor completedWorkExecutor = new ThreadPoolExecutorWithExceptions(1);

        final CompletedWorkChecker workChecker = new CompletedWorkChecker();
        completedWorkExecutor.submit(workChecker);
        completedWorkExecutor.shutdown();

        //thread by surface tile
        final ThreadPoolExecutorWithExceptions tileProcessingExecutor = new ThreadPoolExecutorWithExceptions(numThreads);

        for (final Integer tile : tiles) {
            tileProcessingExecutor.submit(new TileProcessor(tile, barcodesFiles.get(tile)));
        }

        tileProcessingExecutor.shutdown();

        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Reading executor", tileProcessingExecutor, Duration.ofMinutes(5));

        // if there was an exception reading then initiate an immediate shutdown.
        if (tileProcessingExecutor.exception != null) {
            int tasksStillRunning = completedWorkExecutor.shutdownNow().size();
            tasksStillRunning += barcodeWriterThreads.values().stream().mapToLong(executor -> executor.shutdownNow().size()).sum();
            throw new PicardException("Reading executor had exceptions. There were " + tasksStillRunning
                    + " tasks were still running or queued and have been cancelled.", tileProcessingExecutor.exception);
        } else {
            ThreadPoolExecutorUtil.awaitThreadPoolTermination("Tile completion executor", completedWorkExecutor, Duration.ofMinutes(5));
            barcodeWriterThreads.values().forEach(ThreadPoolExecutor::shutdown);
            barcodeWriterThreads.forEach((barcode, executor) -> ThreadPoolExecutorUtil.awaitThreadPoolTermination(barcode + " writer", executor, Duration.ofMinutes(5)));
        }
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
            log.debug("Closing writer for barcode " + barcode);
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
            final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(cbcls, locs, filterFiles, tileNum, barcodeFile);

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                final String barcode = (demultiplex ? cluster.getMatchedBarcode() : null);
                addRecord(barcode, converter.convertClusterToOutputRecord(cluster));
            }

            dataProvider.close();

            final List<RecordWriter> writerList = new ArrayList<>();
            barcodeToRecordCollection.forEach((barcode, value) -> {
                value.doneAdding();
                final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer = barcodeRecordWriterMap.get(barcode);
                writerList.add(new RecordWriter(writer, value, barcode));

            });
            completedWork.put(tileNum, writerList);

            log.info("Finished processing tile " + tileNum);
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
            while (currentTileIndex < tiles.size()) {
                final Integer currentTile = tiles.get(currentTileIndex);
                if (completedWork.containsKey(currentTile)) {
                    log.info("Writing out tile " + currentTile);
                    completedWork.get(currentTile).forEach(writer -> barcodeWriterThreads.get(writer.getBarcode()).submit(writer));
                    currentTileIndex++;
                } else {
                    try {
                        Thread.sleep(5000);
                    } catch (final InterruptedException e) {
                        throw new PicardException(e.getMessage(), e);
                    }
                }
            }

            //we are all done scheduling work.. now schedule the closes
            barcodeRecordWriterMap.forEach((barcode, writer) -> barcodeWriterThreads.get(barcode).submit(new Closer(writer, barcode)));
        }

    }
}
