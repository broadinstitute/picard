package picard.illumina;


import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaFileUtil;
import picard.illumina.parser.ParameterizedFileUtil;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.illumina.parser.readers.LocsFileReader;
import picard.util.IlluminaUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NewIlluminaBasecallsConverter<CLUSTER_OUTPUT_RECORD> extends BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(NewIlluminaBasecallsConverter.class);
    private final Map<String, BarcodeMetric> barcodesMetrics = new HashMap<>();
    private final BarcodeMetric noMatchMetric;
    private final List<File> cbcls;
    private final List<AbstractIlluminaPositionFileReader.PositionInfo> locs = new ArrayList<>();
    private final File[] filterFiles;
    private final int maxNoCalls;
    private final int maxMismatches;
    private final int minMismatchDelta;
    private final int minimumBaseQuality;
    private final MetricsFile<BarcodeMetric, Integer> metrics;
    private final File metricsFile;
    private final Map<String, ThreadPoolExecutorWithExceptions> barcodeWriterThreads = new HashMap<>();
    private final Map<Integer, List<RecordWriter>> completedWork = new HashMap<>();

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
     * @param includeNonPfReads        If true, will include ALL reads (including those which do not have PF set)
     * @param ignoreUnexpectedBarcodes If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap,
     * @param maxNoCalls               The maximum number of no calls to allow in a match.
     * @param maxMismatches            The maximum number of mismatched calls to allow in a match.
     * @param minMismatchDelta         The minimum mismatch difference between the best and second best matches.
     * @param minimumBaseQuality       The minimum base quality to allow for a matching base call.
     * @param metrics                  The metrics output for barcode metrics.
     * @param metricsFile              The metrics file to write out the barcode metrics to.
     */
    public NewIlluminaBasecallsConverter(final File basecallsDir, File barcodesDir, final int lane,
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
                                         final boolean applyEamssFiltering, final boolean includeNonPfReads,
                                         final boolean ignoreUnexpectedBarcodes, int maxNoCalls, int maxMismatches,
                                         int minMismatchDelta, int minimumBaseQuality, MetricsFile<BarcodeMetric, Integer> metrics, File metricsFile) {

        super(barcodeRecordWriterMap, maxReadsInRamPerTile, tmpDirs, codecPrototype, ignoreUnexpectedBarcodes,
                demultiplex, outputRecordComparator, includeNonPfReads, bclQualityEvaluationStrategy,
                outputRecordClass, numProcessors, new IlluminaDataProviderFactory(basecallsDir,
                        barcodesDir, lane, readStructure, bclQualityEvaluationStrategy));
        this.maxNoCalls = maxNoCalls;
        this.maxMismatches = maxMismatches;
        this.minMismatchDelta = minMismatchDelta;
        this.minimumBaseQuality = minimumBaseQuality;
        this.tiles = new ArrayList<>();
        this.metrics = metrics;
        this.metricsFile = metricsFile;
        int numBarcodes = readStructure.sampleBarcodes.length();

        barcodeRecordWriterMap.keySet().forEach(barcode -> {
            if (barcode != null) {
                int pos = 0;
                String[] bcStrings = new String[numBarcodes];
                for (int i = 0; i < numBarcodes; i++) {
                    int endIndex = readStructure.sampleBarcodes.getDescriptorLengths()[i];
                    bcStrings[i] = barcode.substring(pos, endIndex + pos);
                    pos += endIndex;
                }
                this.barcodesMetrics.put(barcode, new BarcodeMetric(null, null, barcode, bcStrings));
            }
            barcodeWriterThreads.put(barcode, new ThreadPoolExecutorWithExceptions(1));
        });

        File laneDir = new File(basecallsDir, IlluminaFileUtil.longLaneStr(lane));

        File[] cycleDirs = IOUtil.getFilesMatchingRegexp(laneDir, IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN);

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
        File locsFile = new File(basecallsDir.getParentFile(), "s.locs");
        LocsFileReader locsFileReader = new LocsFileReader(locsFile);
        while (locsFileReader.hasNext()) {
            locs.add(locsFileReader.next());
        }
        IOUtil.assertFileIsReadable(locsFile);
        //filter

        Pattern laneTileRegex = Pattern.compile(ParameterizedFileUtil.escapePeriods(
                ParameterizedFileUtil.makeLaneTileRegex(".filter", lane)));
        filterFiles = getTiledFiles(laneDir, laneTileRegex);
        for (File filterFile : filterFiles) {
            Matcher tileMatcher = laneTileRegex.matcher(filterFile.getName());
            if (tileMatcher.matches()) {
                tiles.add(Integer.valueOf(tileMatcher.group(1)));
            }
        }
        IOUtil.assertFilesAreReadable(Arrays.asList(filterFiles));
        tiles.sort(TILE_NUMBER_COMPARATOR);

        this.factory.setApplyEamssFiltering(applyEamssFiltering);
        setTileLimits(firstTile, tileLimit);

        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[readStructure.sampleBarcodes.length()];
        int index = 0;
        for (final ReadDescriptor d : readStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                noMatchBarcode[index++] = StringUtil.repeatCharNTimes('N', d.length);
            }
        }

        this.noMatchMetric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

    }

    private File[] getTiledFiles(File baseDirectory, Pattern pattern) {
        return IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
    }

    @Override
    public void doTileProcessing() {

        ThreadPoolExecutor completedWorkExecutor = new ThreadPoolExecutorWithExceptions(1);

        CompletedWorkChecker workChecker = new CompletedWorkChecker();
        completedWorkExecutor.submit(workChecker);
        completedWorkExecutor.shutdown();

        //thread by surface tile
        ThreadPoolExecutor tileProcessingExecutor = new ThreadPoolExecutorWithExceptions(numThreads);

        for (Integer tile : tiles) {
            tileProcessingExecutor.submit(new TileProcessor(tile));
        }

        tileProcessingExecutor.shutdown();

        awaitThreadPoolTermination("Reading executor", tileProcessingExecutor);
        awaitThreadPoolTermination("Tile completion executor", completedWorkExecutor);

        barcodeWriterThreads.values().forEach(ThreadPoolExecutor::shutdown);
        barcodeWriterThreads.forEach((barcode, executor) -> awaitThreadPoolTermination(barcode + " writer", executor));

        if (metricsFile != null) {
            ExtractIlluminaBarcodes.finalizeMetrics(barcodesMetrics, noMatchMetric);

            for (final BarcodeMetric barcodeMetric : barcodesMetrics.values()) {
                metrics.addMetric(barcodeMetric);
            }
            metrics.addMetric(noMatchMetric);
            metrics.write(metricsFile);
            CloserUtil.close(metricsFile);
        }
    }

    private void awaitThreadPoolTermination(String executorName, ThreadPoolExecutor executorService) {
        try {
            while (!executorService.awaitTermination(300, TimeUnit.SECONDS)) {
                log.info(String.format("%s waiting for job completion. Finished jobs - %d : Running jobs - %d : Queued jobs  - %d",
                        executorName, executorService.getCompletedTaskCount(), executorService.getActiveCount(),
                        executorService.getQueue().size()));
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private class RecordWriter implements Runnable {
        private final SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection;
        private final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer;
        private final String barcode;

        RecordWriter(ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer,
                     SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection, String barcode) {
            this.writer = writer;
            this.recordCollection = recordCollection;
            this.barcode = barcode;
        }

        @Override
        public void run() {
            for (CLUSTER_OUTPUT_RECORD record : recordCollection) {
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

        private Closer(ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer, String barcode) {
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

        TileProcessor(int tileNum) {
            this.tileNum = tileNum;
        }

        @Override
        public void run() {
            final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(cbcls, locs, filterFiles, tileNum, barcodesMetrics,
                    noMatchMetric, maxNoCalls, maxMismatches, minMismatchDelta, minimumBaseQuality);

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                final String barcode = (demultiplex ? cluster.getMatchedBarcode() : null);
                addRecord(barcode, converter.convertClusterToOutputRecord(cluster));
            }

            dataProvider.close();

            List<RecordWriter> writerList = new ArrayList<>();
            barcodeToRecordCollection.forEach((barcode, value) -> {
                value.doneAdding();
                ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer = barcodeRecordWriterMap.get(barcode);
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
            return SortingCollection.newInstance(
                    outputRecordClass,
                    codecPrototype.clone(),
                    outputRecordComparator,
                    maxRecordsInRam,
                    tmpDirs);
        }
    }

    private class ThreadPoolExecutorWithExceptions extends ThreadPoolExecutor {
        ThreadPoolExecutorWithExceptions(int threads) {
            super(threads, threads, 0, TimeUnit.SECONDS, new LinkedBlockingDeque<>());
        }

        @Override
        protected void afterExecute(Runnable r, Throwable t) {
            if (t == null && r instanceof Future<?>) {
                try {
                    Future<?> future = (Future<?>) r;
                    if (future.isDone()) {
                        future.get();
                    }
                } catch (CancellationException ce) {
                    t = ce;
                } catch (ExecutionException ee) {
                    t = ee.getCause();
                } catch (InterruptedException ie) {
                    Thread.currentThread().interrupt(); // ignore/reset
                }
            }
            if (t != null) {
                throw new PicardException(t.getMessage(), t);
            }
        }
    }


    private class CompletedWorkChecker implements Runnable {

        private int currentTileIndex = 0;

        @Override
        public void run() {
            while (currentTileIndex < tiles.size()) {
                Integer currentTile = tiles.get(currentTileIndex);
                if (completedWork.containsKey(currentTile)) {
                    log.info("Writing out tile " + currentTile);
                    completedWork.get(currentTile).forEach(writer -> barcodeWriterThreads.get(writer.getBarcode()).submit(writer));
                    currentTileIndex++;
                } else {
                    try {
                        Thread.sleep(5000);
                    } catch (InterruptedException e) {
                        throw new PicardException(e.getMessage(), e);
                    }
                }
            }

            //we are all done scheduling work.. now schedule the closes
            barcodeRecordWriterMap.forEach((barcode, writer) -> barcodeWriterThreads.get(barcode).submit(new Closer(writer, barcode)));
        }

    }
}
