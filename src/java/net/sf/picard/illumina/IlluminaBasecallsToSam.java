/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

package net.sf.picard.illumina;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.illumina.parser.*;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.*;
import net.sf.picard.util.IlluminaUtil.IlluminaAdapterPair;
import net.sf.samtools.*;
import net.sf.samtools.util.Iso8601Date;
import net.sf.samtools.util.PeekIterator;
import net.sf.samtools.util.SortingCollection;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicBoolean;

import static java.util.concurrent.TimeUnit.MILLISECONDS;

/**
 * IlluminaBasecallsToSam transforms a lane of Illumina data file formats (bcl, locs, clocs, qseqs, etc.) into
 * SAM or BAM file format.
 * <p/>
 * In this application, barcode data is read from Illumina data file groups, each of which is associated with a tile.
 * Each tile may contain data for any number of barcodes, and a single barcode's data may span multiple tiles.  Once the
 * barcode data is collected from files, each barcode's data is written to its own SAM/BAM.  The barcode data must be
 * written in order; this means that barcode data from each tile is sorted before it is written to file, and that if a
 * barcode's data does span multiple tiles, data collected from each tile must be written in the order of the tiles
 * themselves.
 * <p/>
 * This class employs a number of private subclasses to achieve this goal.  The TileReadAggregator controls the flow
 * of operation.  It is fed a number of Tiles which it uses to spawn TileReaders.  TileReaders are responsible for
 * reading Illumina data for their respective tiles from disk, and as they collect that data, it is fed back into the
 * TileReadAggregator.  When a TileReader completes a tile, it notifies the TileReadAggregator, which reviews what was
 * read and conditionally queues its writing to disk, baring in mind the requirements of write-order described in the
 * previous paragraph.  As writes complete, the TileReadAggregator re-evalutes the state of reads/writes and may queue
 * more writes.  When all barcodes for all tiles have been written, the TileReadAggregator shuts down.
 * <p/>
 * The TileReadAggregator controls task execution using a specialized ThreadPoolExecutor.  It accepts special Runnables
 * of type PriorityRunnable which allow a priority to be assigned to the runnable.  When the ThreadPoolExecutor is
 * assigning threads, it gives priority to those PriorityRunnables with higher priority values.  In this application,
 * TileReaders are assigned lowest priority, and write tasks are assigned high priority.  It is designed in this fashion
 * to minimize the amount of time data must remain in memory (write the data as soon as possible, then discard it from
 * memory) while maximizing CPU usage.
 *
 * @author jburke@broadinstitute.org
 * @author mccowan@broadinstitute.org
 */
public class IlluminaBasecallsToSam extends CommandLineProgram {
    // The following attributes define the command-line arguments
    @Usage
    public String USAGE =
            getStandardUsagePreamble() + "Generate a SAM or BAM file from data in an Illumina basecalls output directory.\n";

    @Option(doc = "The basecalls output directory. ", shortName = "B")
    public File BASECALLS_DIR;
    @Option(doc = "Lane number. ", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;
    @Option(doc = "Deprecated (use LIBRARY_PARAMS).  The output SAM or BAM file. Format is determined by extension.",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            mutex = {"BARCODE_PARAMS", "LIBRARY_PARAMS"})
    public File OUTPUT;

    @Option(doc = "Prefixed to read names.")
    public String RUN_BARCODE;
    @Option(doc = "Deprecated (use LIBRARY_PARAMS).  The name of the sequenced sample",
            shortName = StandardOptionDefinitions.SAMPLE_ALIAS_SHORT_NAME,
            mutex = {"BARCODE_PARAMS", "LIBRARY_PARAMS"})
    public String SAMPLE_ALIAS;
    @Option(doc = "ID used to link RG header record with RG tag in SAM record.  " +
            "If these are unique in SAM files that get merged, merge performance is better.  " +
            "If not specified, READ_GROUP_ID = <first 5 chars of RUN_BARCODE>.<LANE> .",
            shortName = StandardOptionDefinitions.READ_GROUP_ID_SHORT_NAME, optional = true)
    public String READ_GROUP_ID;
    @Option(doc = "Deprecated (use LIBRARY_PARAMS).  The name of the sequenced library",
            shortName = StandardOptionDefinitions.LIBRARY_NAME_SHORT_NAME,
            optional = true,
            mutex = {"BARCODE_PARAMS", "LIBRARY_PARAMS"})
    public String LIBRARY_NAME;
    @Option(doc = "The name of the sequencing center that produced the reads to fill in the RG.CN tag.", optional = true)
    public String SEQUENCING_CENTER = "BI";
    @Option(doc = "The start date of the run.", optional = true)
    public Date RUN_START_DATE;
    @Option(doc = "The name of the sequencing technology that produced the read.", optional = true)
    public String PLATFORM = "illumina";

    @Option(doc = ReadStructure.PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Option(doc = "Deprecated (use LIBRARY_PARAMS).  Tab-separated file for creating all output BAMs for barcoded run " +
            "with single IlluminaBasecallsToSam invocation.  Columns are BARCODE, OUTPUT, SAMPLE_ALIAS, and " +
            "LIBRARY_NAME.  Row with BARCODE=N is used to specify a file for no barcode match",
            mutex = {"OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "LIBRARY_PARAMS"})
    public File BARCODE_PARAMS;

    @Option(doc = "Tab-separated file for creating all output BAMs for a run with single IlluminaBasecallsToSam " +
            "invocation.  TheColumns are OUTPUT, SAMPLE_ALIAS, and LIBRARY_NAME, BARCODE_1, BARCODE_2 ... BARCODE_X " +
            "where X = number of barcodes per cluster (optional).  Row with BARCODE_1=N is used to specify a file " +
            "for no barcode match.  You may also provide any 2 letter RG header attributes (excluding PU, CN, PL, and" +
            " DT)  as columns in this file and the values for those columns will be inserted into the RG tag for the" +
            " BAM file created for a given row.",
            mutex = {"OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "BARCODE_PARAMS"})
    public File LIBRARY_PARAMS;

    @Option(doc = "Which adapters to look for in the read.")
    public List<IlluminaAdapterPair> ADAPTERS_TO_CHECK = new ArrayList<IlluminaAdapterPair>(
            Arrays.asList(IlluminaAdapterPair.INDEXED,
                    IlluminaAdapterPair.DUAL_INDEXED, IlluminaAdapterPair.NEXTERA_V2));

    @Option(doc = "Run this many threads in parallel. If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0, then the number of cores used will" +
            " be the number available on the machine less NUM_PROCESSORS.")
    public Integer NUM_PROCESSORS = 0;
    @Option(doc = "If set, this is the first tile to be processed (for debugging).  Note that tiles are not processed" +
            " in numerical order.",
            optional = true)
    public Integer FIRST_TILE;
    @Option(doc = "If set, process no more than this many tiles (for debugging).", optional = true)
    public Integer TILE_LIMIT;

    @Option(doc = "If true, call System.gc() periodically.  This is useful in cases in which the -Xmx value passed " +
            "is larger than the available memory.  Default: true.", optional = true)
    public Boolean FORCE_GC = true;
    @Option(doc = "Configure SortingCollections to store this many records before spilling to disk. For an indexed" +
            " run, each SortingCollection gets this value/number of indices.")
    public int MAX_READS_IN_RAM_PER_TILE = 1200000;

    private final Map<String, SAMFileWriter> barcodeSamWriterMap = new HashMap<String, SAMFileWriter>();
    private IlluminaBasecallsToSamConverter converter;
    private IlluminaDataProviderFactory factory;
    private ReadStructure readStructure;
    private int numThreads;
    private final Comparator<SAMRecord> samRecordQueryNameComparator = new SAMRecordQueryNameComparator();
    private List<Integer> tiles;
    public static final IlluminaDataType[] DATA_TYPES_NO_BARCODE =
            {IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.Position, IlluminaDataType.PF};
    private static final IlluminaDataType[] DATA_TYPES_WITH_BARCODE = Arrays.copyOf(DATA_TYPES_NO_BARCODE, DATA_TYPES_NO_BARCODE.length + 1);
    private static final Log log = Log.getInstance(IlluminaBasecallsToSam.class);
    private final ProgressLogger readProgressLogger = new ProgressLogger(log, 1000000, "Read");
    private final ProgressLogger writeProgressLogger = new ProgressLogger(log, 1000000, "Write");

    /**
     * Describes the state of a barcode's data's processing in the context of a tile.  It is either not available in
     * that tile, has been read, has been queued to be written to file, or has been written to file.  A barcode only
     * takes on a state once the tile (which is serving as the context of this state) has been read.
     */
    private enum TileBarcodeProcessingState {
        NA, READ, QUEUED_FOR_WRITE, WRITTEN
    }

    /**
     * Describes the state of a tile being processed.  It is either not yet completely read, or read.
     */
    private enum TileProcessingState {
        NOT_DONE_READING, DONE_READING
    }

    /**
     * A comparator for tile numbers, which are not necessarily ordered by the number's value.
     */
    public static final Comparator<Integer> TILE_NUMBER_COMPARATOR = new Comparator<Integer>() {
        @Override
        public int compare(final Integer integer1, final Integer integer2) {
            final String s1 = integer1.toString();
            final String s2 = integer2.toString();
            // Because a the tile number is followed by a colon, a tile number that
            // is a prefix of another tile number should sort after. (e.g. 10 sorts after 100).
            if (s1.length() < s2.length()) {
                if (s2.startsWith(s1)) {
                    return 1;
                }
            } else if (s2.length() < s1.length()) {
                if (s1.startsWith(s2)) {
                    return -1;
                }
            }
            return s1.compareTo(s2);
        }
    };

    static {
        DATA_TYPES_WITH_BARCODE[DATA_TYPES_WITH_BARCODE.length - 1] = IlluminaDataType.Barcodes;
    }

    /**
     * Simple representation of a tile
     */
    private class Tile implements Comparable<Tile> {
        private final int tileNumber;

        public Tile(final int i) {
            tileNumber = i;
        }

        public int getNumber() {
            return tileNumber;
        }

        @Override
        public boolean equals(final Object o) {
            return o instanceof Tile && this.getNumber() == ((Tile) o).getNumber();
        }

        @Override
        public int compareTo(final Tile o) {
            return TILE_NUMBER_COMPARATOR.compare(this.getNumber(), o.getNumber());
        }
    }

    /**
     * Reads the information from a tile via an IlluminaDataProvider and feeds red information into a processingRecord
     * managed by the TileReadAggregator.
     */
    private class TileReader {
        private final Tile tile;
        private final TileReadAggregator handler;
        private final TileProcessingRecord processingRecord;

        public TileReader(final Tile tile, final TileReadAggregator handler, final TileProcessingRecord processingRecord) {
            this.tile = tile;
            this.handler = handler;
            this.processingRecord = processingRecord;
        }

        /**
         * Reads the data from the appropriate IlluminaDataProvider and feeds it into the TileProcessingRecord for
         * this tile.
         */
        public void process() {
            final IlluminaDataProvider dataProvider = factory.makeDataProvider(Arrays.asList(this.tile.getNumber()));
            final SAMRecord[] recordContainer = new SAMRecord[converter.getNumRecordsPerCluster()];

            log.debug(String.format("Reading data from tile %s ...", tile.getNumber()));

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                final String barcode = cluster.getMatchedBarcode();

                converter.createSamRecords(cluster, null, recordContainer);

                for (final SAMRecord record : recordContainer) {
                    readProgressLogger.record(record);
                    this.processingRecord.addRecord(barcode, record);
                }
            }

            this.handler.completeTile(this.tile);
        }
    }

    /**
     * Represents the state of a tile's processing and encapsulates the data collected from that tile
     */
    private class TileProcessingRecord {
        final private Map<String, SortingCollection<SAMRecord>> barcodeToRecordCollection = new HashMap<String, SortingCollection<SAMRecord>>();
        final private Map<String, TileBarcodeProcessingState> barcodeToProcessingState = new HashMap<String, TileBarcodeProcessingState>();
        private TileProcessingState state = TileProcessingState.NOT_DONE_READING;
        private long recordCount = 0;

        /**
         * Returns the state of this tile's processing.
         */
        public TileProcessingState getState() {
            return this.state;
        }

        /**
         * Sets the state of this tile's processing.
         */
        public void setState(final TileProcessingState state) {
            this.state = state;
        }

        /**
         * Adds the provided recoded to this tile.
         */
        public void addRecord(final String barcode, final SAMRecord... records) {
            this.recordCount += records.length;

            // Grab the existing collection, or initialize it if it doesn't yet exist
            SortingCollection<SAMRecord> recordCollection = this.barcodeToRecordCollection.get(barcode);
            if (recordCollection == null) {
                recordCollection = this.newSortingCollection();
                this.barcodeToRecordCollection.put(barcode, recordCollection);
                this.barcodeToProcessingState.put(barcode, null);
            }

            for (final SAMRecord record : records) {
                recordCollection.add(record);
            }

        }

        private SortingCollection<SAMRecord> newSortingCollection() {
            final int maxRecordsInRam =
                    IlluminaBasecallsToSam.this.MAX_READS_IN_RAM_PER_TILE /
                            IlluminaBasecallsToSam.this.barcodeSamWriterMap.size();
            return SortingCollection.newInstance(
                    SAMRecord.class,
                    new BAMRecordCodec(null),
                    new SAMRecordQueryNameComparator(),
                    maxRecordsInRam,
                    IlluminaBasecallsToSam.this.TMP_DIR);
        }

        /**
         * Returns the number of unique barcodes read.
         */
        public long getBarcodeCount() {
            return this.barcodeToRecordCollection.size();
        }

        /**
         * Returns the number of records read.
         */
        public long getRecordCount() {
            return recordCount;
        }

        /**
         * Returns the mapping of barcodes to records associated with them.
         */
        public Map<String, SortingCollection<SAMRecord>> getBarcodeRecords() {
            return barcodeToRecordCollection;
        }

        /**
         * Gets the state of the provided barcode's data's processing progress.  Only invoke this query if this tile
         * is in a DONE_READING state.
         *
         * @throws IllegalStateException When a barcode is queried before the tile is in the DONE_READING state
         */
        public TileBarcodeProcessingState getBarcodeState(final String barcode) {
            if (this.getState() == TileProcessingState.NOT_DONE_READING) {
                throw new IllegalStateException(
                        "A tile's barcode data's state cannot be queried until the tile has been completely read.");
            }

            if (this.barcodeToProcessingState.containsKey(barcode)) {
                return this.barcodeToProcessingState.get(barcode);
            } else {
                return TileBarcodeProcessingState.NA;
            }
        }

        /**
         * Sets the processing state of the provided barcode in this record.
         *
         * @throws NoSuchElementException When the provided barcode is not one associated with this record.
         */
        public void setBarcodeState(final String barcode, final TileBarcodeProcessingState state) {
            if (this.barcodeToProcessingState.containsKey(barcode)) {
                this.barcodeToProcessingState.put(barcode, state);
            } else {
                throw new NoSuchElementException(String.format("No record of the provided barcode, %s.", barcode));
            }
        }

        /**
         * Returns the distinct set of barcodes for which data has been collected in this record.
         *
         * @return
         */
        public Set<String> getBarcodes() {
            return this.getBarcodeRecords().keySet();
        }
    }

    /**
     * Aggregates data collected from tiles and writes them to file. Accepts records from TileReaders and maps
     * them to the appropriate BAM writers.
     */
    private class TileReadAggregator {
        /**
         * The collection of records associated with a particular tile.
         * <p/>
         * Implemented as a TreeMap to guarantee tiles are iterated over in natural order.
         */
        private final Map<Tile, TileProcessingRecord> tileRecords = new TreeMap<Tile, TileProcessingRecord>();

        /**
         * The executor responsible for doing work.
         * <p/>
         * Implemented as a ThreadPoolExecutor with a PriorityBlockingQueue which orders submitted Runnables by their
         * priority.
         */
        private final ExecutorService prioritizingThreadPool = new ThreadPoolExecutor(
                IlluminaBasecallsToSam.this.numThreads,
                IlluminaBasecallsToSam.this.numThreads,
                0L,
                MILLISECONDS,
                new PriorityBlockingQueue<Runnable>(5, new Comparator<Runnable>() {
                    @Override
                    /**
                     * Compare the two Runnables, and assume they are PriorityRunnable; if not something strange is
                     * going on, so allow a ClassCastException be thrown.
                     */
                    public int compare(final Runnable o1, final Runnable o2) {
                        // Higher priority items go earlier in the queue, so reverse the "natural" comparison.
                        return ((PriorityRunnable) o2).getPriority() - ((PriorityRunnable) o1).getPriority();
                    }
                }));

        /**
         * The object acting as a latch to notify when the aggregator completes its work.
         */
        private final Object completionLatch = new Object();

        /**
         * Stores the thread that is executing this work so that it can be interrupted upon failure.
         */
        private Thread parentThread;
        private final Object workEnqueueMonitor = new Object();
        private final AtomicBoolean submitted = new AtomicBoolean(false);


        /**
         * Creates a TileReadAggregator that reads from the provided tiles.
         *
         * @param tiles
         */
        public TileReadAggregator(final Collection<Tile> tiles) {
            for (final Tile t : tiles) {
                tileRecords.put(t, new TileProcessingRecord());
            }
        }

        /**
         * Execute the tile aggregator's work.  Creates a thread pool to read data from tiles and write them to file.
         * Invoke this method only once.
         *
         * @throws IllegalStateException If submit was called more than once.
         */
        public void submit() {
            // Ensure the aggregator as not yet been submitted
            if (!this.submitted.compareAndSet(false, true)) {
                throw new IllegalStateException("The submit() method may not be called more than once.");
            }

            // Set the thread that is executing this work
            this.parentThread = Thread.currentThread();

            /**
             * For each tile, create and submit a tile processor.  Give it a negative execution priority (so that
             * prioritized tasks with a positive execution priority execute first), and give later tiles a lesser
             * (more negative) priority.
             */
            int priority = 0;
            for (final Tile tile : this.tileRecords.keySet()) {
                final TileReader reader = new TileReader(tile, this, this.tileRecords.get(tile));
                this.prioritizingThreadPool.execute(new PriorityRunnable(--priority) {
                    @Override
                    public void run() {
                        reader.process();
                    }
                });
            }
        }

        /**
         * Signals that a tile's processing is complete.  This must be invoked exactly once per tile, and only after
         * all of that tile has been processed.
         *
         * @throws IllegalStateException When the tile is already in the completed state.
         */
        private void completeTile(final Tile tile) {
            final TileProcessingRecord tileRecord = this.tileRecords.get(tile);

            if (tileRecord.getState() == TileProcessingState.DONE_READING) {
                throw new IllegalStateException("This tile is already in the completed state.");
            }

            // Update all of the barcodes and the tile to be marked as read
            for (final String barcode : tileRecord.getBarcodes()) {
                tileRecord.setBarcodeState(barcode, TileBarcodeProcessingState.READ);
            }
            tileRecord.setState(TileProcessingState.DONE_READING);

            log.debug(String.format("Completed reading tile %s; collected %s reads spanning %s barcodes.",
                    tile.getNumber(), tileRecord.getRecordCount(), tileRecord.getBarcodeCount()));

            //noinspection SynchronizationOnLocalVariableOrMethodParameter
            this.findAndEnqueueWorkOrSignalCompletion();
        }

        /**
         * Blocks until this aggregator completes its work.
         *
         * @throws InterruptedException
         */
        public void awaitWorkComplete() throws InterruptedException {
            synchronized (this.completionLatch) {
                this.completionLatch.wait();
            }
        }

        /**
         * Signals to any thread awaiting via awaitWorkComplete() that no work remains. Called
         * when this aggregator has reached its completed state.
         */
        private void signalWorkComplete() {
            synchronized (this.completionLatch) {
                this.completionLatch.notifyAll();
            }
        }

        /**
         * Poll the aggregator to find more tasks for it to enqueue.  Specifically, searches for un-written data
         * read from tiles for each barcode and enqueues it for writing.
         */
        private void findAndEnqueueWorkOrSignalCompletion() {
            synchronized (this.workEnqueueMonitor) {
                /**
                 * If there is work remaining to be done in this aggregator, walk through all of the barcodes and find
                 * tiles which have not yet written their barcode data but are in a state where they are able to.
                 */
                if (this.isWorkCompleted()) {
                    this.signalWorkComplete();
                } else {
                    final Queue<Runnable> tasks = new LinkedList<Runnable>();
                    for (final String barcode : barcodeSamWriterMap.keySet()) {
                        NEXT_BARCODE:
                        for (final Map.Entry<Tile, TileProcessingRecord> entry : this.tileRecords.entrySet()) {
                            final Tile tile = entry.getKey();
                            final TileProcessingRecord tileRecord = entry.getValue();

                            /**
                             * If this tile has not been read, we cannot write this or later tiles' barcode data;
                             * move to the next barcode.
                             */
                            if (tileRecord.getState() != TileProcessingState.DONE_READING) {
                                break NEXT_BARCODE;
                            }
                            switch (tileRecord.getBarcodeState(barcode)) {
                                case NA:
                                case WRITTEN:
                                    /**
                                     * There is no data for this barcode for this tile, or it is already written; in
                                     * either scenario, this barcode will not be processed further for this tile, so
                                     * move onto the next tile as a possible candidate.
                                     */
                                    continue;
                                case QUEUED_FOR_WRITE:
                                    /**
                                     * The write for this barcode is in progress for this tile, so skip to the next
                                     * barcode.
                                     */
                                    break NEXT_BARCODE;
                                case READ:
                                    /**
                                     * This barcode has beenr read, and all of the earlier tiles have been written
                                     * for this barcode, so queue its writing.
                                     */
                                    tileRecord.setBarcodeState(barcode, TileBarcodeProcessingState.QUEUED_FOR_WRITE);
                                    log.debug(String.format("Enqueuing work for tile %s and barcode %s.", tile.getNumber(), barcode));
                                    tasks.add(this.newBarcodeWorkInstance(tile, tileRecord, barcode));
                                    break NEXT_BARCODE;
                            }
                        }
                    }

                    for (final Runnable task : tasks) {
                        this.prioritizingThreadPool.execute(task);
                    }
                }
            }
        }

        /**
         * Returns a PriorityRunnable that encapsulates the work involved with writing the provided tileRecord's data
         * for the given barcode to disk.
         *
         * @param tile       The tile from which the record was read
         * @param tileRecord The processing record associated with the tile
         * @param barcode    The barcode whose data within the tileRecord is to be written
         * @return The runnable that upon invocation writes the barcode's data from the tileRecord to disk
         */
        private PriorityRunnable newBarcodeWorkInstance(final Tile tile, final TileProcessingRecord tileRecord, final String barcode) {
            return new PriorityRunnable() {
                @Override
                public void run() {
                    try {
                        final SortingCollection<SAMRecord> records = tileRecord.getBarcodeRecords().get(barcode);
                        final SAMFileWriter writer = barcodeSamWriterMap.get(barcode);

                        log.debug(String.format("Writing records from tile %s with barcode %s ...", tile.getNumber(), barcode));

                        final PeekIterator<SAMRecord> it = new PeekIterator<SAMRecord>(records.iterator());
                        while (it.hasNext()) {
                            final SAMRecord rec = it.next();

                            /**
                             * PIC-330 Sometimes there are two reads with the same cluster coordinates, and thus
                             * the same read name.  Discard both of them.  This code assumes that the two first of pairs
                             * will come before the two second of pairs, so it isn't necessary to look ahead a different
                             * distance for paired end.  It also assumes that for paired ends there will be duplicates
                             * for both ends, so there is no need to be PE-aware.
                             */
                            if (it.hasNext()) {
                                final SAMRecord lookAhead = it.peek();
                                if (!rec.getReadUnmappedFlag() || !lookAhead.getReadUnmappedFlag()) {
                                    throw new IllegalStateException("Should not have mapped reads.");
                                }

                                if (samRecordQueryNameComparator.compare(rec, lookAhead) == 0) {
                                    it.next();
                                    log.info("Skipping reads with identical read names: " + rec.getReadName());
                                    continue;
                                }
                            }

                            writer.addAlignment(rec);
                            writeProgressLogger.record(rec);
                        }

                        tileRecord.setBarcodeState(barcode, TileBarcodeProcessingState.WRITTEN);
                        findAndEnqueueWorkOrSignalCompletion();

                    } catch (RuntimeException e) {
                        /**
                         * In the event of an internal failure, signal to the parent thread that something has gone
                         * wrong.  This is necessary because if an item of work fails to complete, the aggregator will
                         * will never reach its completed state, and it will never terminate.
                         */
                        parentThread.interrupt();
                        throw e;
                    }
                }

            };
        }

        /**
         * Returns true if this aggregator has not completed its work.  Specifically, returns false iff
         * any tile's barcode data yas not yet been written.
         *
         * @return True if more work remains to be done, false otherwise
         */
        public boolean isWorkCompleted() {
            for (final Map.Entry<Tile, TileProcessingRecord> entry : this.tileRecords.entrySet()) {
                final TileProcessingRecord tileProcessingRecord = entry.getValue();

                if (tileProcessingRecord.getState() != TileProcessingState.DONE_READING) {
                    return false;
                } else {
                    for (final TileBarcodeProcessingState barcodeProcessingState : tileProcessingRecord.barcodeToProcessingState.values()) {
                        if (barcodeProcessingState != TileBarcodeProcessingState.WRITTEN) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        /**
         * Terminates the threads currently exiting in the thread pool abruptly via ThreadPoolExecutor.shutdownNow().
         */
        public void shutdown() {
            this.prioritizingThreadPool.shutdownNow();
        }
    }

    /**
     * A Runnable that carries a priority which is used to compare and order other PriorityRunnables in a task queue.
     */
    private abstract class PriorityRunnable implements Runnable {
        private final int priority;

        /**
         * Create a new priority runnable with a default priority of 1.
         */
        public PriorityRunnable() {
            this(1);
        }

        public PriorityRunnable(final int priority) {
            this.priority = priority;
        }

        /**
         * Returns the priority level.  Higher priorities are run earlier.
         *
         * @return
         */
        int getPriority() {
            return this.priority;
        }
    }

    @Override
    protected int doWork() {
        initialize();

        doTileProcessing();

        finalise();

        return 0;
    }

    /**
     * Closes the SAMFileWriters in barcodeSamWriterMap.
     * <p/>
     * This method is intentionally misspelled to avoid conflict with Object.finalize();
     */
    private void finalise() {
        // Close the writers
        for (final Map.Entry<String, SAMFileWriter> entry : barcodeSamWriterMap.entrySet()) {
            final SAMFileWriter writer = entry.getValue();
            log.debug(String.format("Closing file for barcode %s.", entry.getKey()));
            writer.close();
        }
    }

    /**
     * Prepares loggers, initiates garbage collection thread, parses arguments and initialized variables appropriately/
     */
    private void initialize() {
        if (OUTPUT != null) {
            IoUtil.assertFileIsWritable(OUTPUT);
        }

        if (LIBRARY_PARAMS != null) {
            IoUtil.assertFileIsReadable(LIBRARY_PARAMS);
        }

        // If we're forcing garbage collection, collect every 5 minutes in a daemon thread.
        if (this.FORCE_GC) {
            final Timer gcTimer = new Timer(true);
            final long delay = 5 * 1000 * 60;
            gcTimer.scheduleAtFixedRate(new TimerTask() {
                @Override
                public void run() {
                    System.out.println("Before explicit GC, Runtime.totalMemory()=" + Runtime.getRuntime().totalMemory());
                    System.gc();
                    System.runFinalization();
                    System.out.println("After explicit GC, Runtime.totalMemory()=" + Runtime.getRuntime().totalMemory());
                }
            }, delay, delay);
        }

        readStructure = new ReadStructure(READ_STRUCTURE);
        factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, getDataTypesFromReadStructure(readStructure));

        log.info("DONE_READING STRUCTURE IS " + readStructure.toString());

        this.tiles = new ArrayList<Integer>(factory.getAvailableTiles());
        // Since the first non-fixed part of the read name is the tile number, without preceding zeroes,
        // and the output is sorted by read name, process the tiles in this order.
        Collections.sort(tiles, TILE_NUMBER_COMPARATOR);
        if (FIRST_TILE != null) {
            int i;
            for (i = 0; i < tiles.size(); ++i) {
                if (tiles.get(i).intValue() == FIRST_TILE.intValue()) {
                    tiles = tiles.subList(i, tiles.size());
                    break;
                }
            }
            if (tiles.get(0).intValue() != FIRST_TILE.intValue()) {
                throw new PicardException("FIRST_TILE=" + FIRST_TILE + ", but that tile was not found.");
            }
        }
        if (TILE_LIMIT != null && tiles.size() > TILE_LIMIT) {
            tiles = tiles.subList(0, TILE_LIMIT);
        }

        if (OUTPUT != null) {
            barcodeSamWriterMap.put(null, buildSamFileWriter(OUTPUT, SAMPLE_ALIAS, LIBRARY_NAME, buildSamHeaderParameters(null)));
        } else {
            populateWritersFromLibraryParams();
        }

        /**
         * Be sure to pass the outputReadStructure to IlluminaBasecallsToSamConverter, which reflects the structure of the output cluster
         * data which may be different from the input read structure (specifically if there are skips).
         */
        converter = new IlluminaBasecallsToSamConverter(RUN_BARCODE, READ_GROUP_ID, factory.getOutputReadStructure(), ADAPTERS_TO_CHECK);

        if (NUM_PROCESSORS == 0) {
            this.numThreads = Runtime.getRuntime().availableProcessors();
        } else if (NUM_PROCESSORS < 0) {
            this.numThreads = Runtime.getRuntime().availableProcessors() + NUM_PROCESSORS;
        } else {
            this.numThreads = NUM_PROCESSORS;
        }
        this.numThreads = Math.max(1, Math.min(this.numThreads, tiles.size()));
    }

    private void doTileProcessing() {
        // TODO: Eliminate this when switch to JDK 7
        FileChannelJDKBugWorkAround.doBugWorkAround();

        // Generate the list of tiles that will be processed
        final List<Tile> tiles = new ArrayList<Tile>();
        for (final Integer tileNumber : this.tiles) {
            tiles.add(new Tile(tileNumber));
        }

        final TileReadAggregator tileReadAggregator = new TileReadAggregator(tiles);
        tileReadAggregator.submit();
        try {
            tileReadAggregator.awaitWorkComplete();
        } catch (InterruptedException e) {
            log.error(e, "Attempting to shut down worker threads ...");
            tileReadAggregator.shutdown();
            throw new PicardException(String.format("Main thread interrupted: %s.", e.getMessage()));
        }
    }

    /**
     * Assert that expectedCols are present and return actualCols - expectedCols
     *
     * @param actualCols   The columns present in the LIBRARY_PARAMS file
     * @param expectedCols The columns that are REQUIRED
     * @return actualCols - expectedCols
     */
    private Set<String> findAndFilterExpectedColumns(final Set<String> actualCols, final Set<String> expectedCols) {
        final Set<String> missingColumns = new HashSet<String>(expectedCols);
        missingColumns.removeAll(actualCols);

        if (missingColumns.size() > 0) {
            throw new PicardException(String.format(
                    "LIBRARY_PARAMS file %s is missing the following columns: %s.",
                    LIBRARY_PARAMS.getAbsolutePath(), StringUtil.join(", ", missingColumns
            )));
        }

        final Set<String> remainingColumns = new HashSet<String>(actualCols);
        remainingColumns.removeAll(expectedCols);
        return remainingColumns;
    }

    /**
     * Given a set of columns assert that all columns conform to the format of an RG header attribute (i.e. 2 letters)
     * the attribute is NOT a member of the rgHeaderTags that are built by default in buildSamHeaderParameters
     *
     * @param rgTagColumns A set of columns that should conform to the rg header attribute format
     */
    private void checkRgTagColumns(final Set<String> rgTagColumns) {
        final Set<String> forbiddenHeaders = buildSamHeaderParameters(null).keySet();
        forbiddenHeaders.retainAll(rgTagColumns);

        if (forbiddenHeaders.size() > 0) {
            throw new PicardException("Illegal ReadGroup tags in library params(barcode params) file(" + LIBRARY_PARAMS.getAbsolutePath() + ") Offending headers = " + StringUtil.join(", ", forbiddenHeaders));
        }

        for (final String column : rgTagColumns) {
            if (column.length() > 2) {
                throw new PicardException("Column label (" + column + ") unrecognized.  Library params(barcode params) can only contain the columns " +
                        "(OUTPUT, LIBRARY_NAME, SAMPLE_ALIAS, BARCODE, BARCODE_<X> where X is a positive integer) OR two letter RG tags!");
            }
        }
    }

    /**
     * For each line in the LIBRARY_PARAMS file create a SamFileWriter and put it in the barcodeSamWriterMap map, where
     * the key to the map is the concatenation of all barcodes in order for the given line
     */
    private void populateWritersFromLibraryParams() {
        final TabbedTextFileWithHeaderParser libraryParamsParser = new TabbedTextFileWithHeaderParser(LIBRARY_PARAMS);

        final Set<String> expectedColumnLabels = CollectionUtil.makeSet("OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME");
        final List<String> barcodeColumnLabels = new ArrayList<String>();
        if (readStructure.barcodes.length() == 1) {
            //For the single barcode read case, the barcode label name can either by BARCODE or BARCODE_1
            if (libraryParamsParser.hasColumn("BARCODE")) {
                barcodeColumnLabels.add("BARCODE");
            } else if (libraryParamsParser.hasColumn("BARCODE_1")) {
                barcodeColumnLabels.add("BARCODE_1");
            } else {
                throw new PicardException("LIBRARY_PARAMS(BARCODE_PARAMS) file " + LIBRARY_PARAMS + " does not have column BARCODE or BARCODE_1.");
            }
        } else {
            for (int i = 1; i <= readStructure.barcodes.length(); i++) {
                barcodeColumnLabels.add("BARCODE_" + i);
            }
        }

        expectedColumnLabels.addAll(barcodeColumnLabels);
        final Set<String> rgTagColumns = findAndFilterExpectedColumns(libraryParamsParser.columnLabels(), expectedColumnLabels);
        checkRgTagColumns(rgTagColumns);

        for (final TabbedTextFileWithHeaderParser.Row row : libraryParamsParser) {
            List<String> barcodeValues = null;

            if (barcodeColumnLabels.size() > 0) {
                barcodeValues = new ArrayList<String>();
                for (final String barcodeLabel : barcodeColumnLabels) {
                    barcodeValues.add(row.getField(barcodeLabel));
                }
            }

            final String key = (barcodeValues == null || barcodeValues.contains("N")) ? null : StringUtil.join("", barcodeValues);
            if (barcodeSamWriterMap.containsKey(key)) {    //This will catch the case of having more than 1 line in a non-barcoded LIBRARY_PARAMS file
                throw new PicardException("Row for barcode " + key + " appears more than once in LIBRARY_PARAMS or BARCODE_PARAMS file " +
                        LIBRARY_PARAMS);
            }

            final Map<String, String> samHeaderParams = buildSamHeaderParameters(barcodeValues);

            for (final String tagName : rgTagColumns) {
                samHeaderParams.put(tagName, row.getField(tagName));
            }

            final SAMFileWriter writer = buildSamFileWriter(new File(row.getField("OUTPUT")),
                    row.getField("SAMPLE_ALIAS"), row.getField("LIBRARY_NAME"), samHeaderParams);
            barcodeSamWriterMap.put(key, writer);
        }
        if (barcodeSamWriterMap.isEmpty()) {
            throw new PicardException("LIBRARY_PARAMS(BARCODE_PARAMS) file " + LIBRARY_PARAMS + " does have any data rows.");
        }
    }

    /**
     * Create the list of headers that will be added to the SAMFileHeader for a library with the given barcodes (or
     * the entire run if barcodes == NULL).  Note that any value that is null will NOT be added via buildSamFileWriter
     * but is placed in the map in order to be able to query the tags that we automatically add.
     *
     * @param barcodes The list of barcodes that uniquely identify the read group we are building parameters for
     * @return A Map of ReadGroupHeaderTags -> Values
     */
    private Map<String, String> buildSamHeaderParameters(final List<String> barcodes) {
        final Map<String, String> params = new HashMap<String, String>();

        String platformUnit = RUN_BARCODE + "." + LANE;
        if (barcodes != null) platformUnit += ("." + IlluminaUtil.barcodeSeqsToString(barcodes));
        params.put("PU", platformUnit);

        params.put("CN", SEQUENCING_CENTER);
        params.put("PL", PLATFORM);
        if (RUN_START_DATE != null) {
            final Iso8601Date date = new Iso8601Date(RUN_START_DATE);
            params.put("DT", date.toString());
        } else {
            params.put("DT", null);
        }

        return params;
    }

    /**
     * Build a SamFileWriter that will output it's contents to output.
     *
     * @param output           The file to which to write
     * @param sampleAlias      The sample alias set in the read group header
     * @param libraryName      The name of the library to which this read group belongs
     * @param headerParameters Header parameters that will be added to the RG header for this SamFile
     * @return A SAMFileWriter
     */
    private SAMFileWriter buildSamFileWriter(final File output, final String sampleAlias, final String libraryName, final Map<String, String> headerParameters) {
        IoUtil.assertFileIsWritable(output);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(READ_GROUP_ID);
        rg.setSample(sampleAlias);

        if (libraryName != null) rg.setLibrary(libraryName);
        for (final Map.Entry<String, String> tagNameToValue : headerParameters.entrySet()) {
            if (tagNameToValue.getValue() != null) {
                rg.setAttribute(tagNameToValue.getKey(), tagNameToValue.getValue());
            }
        }

        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        header.addReadGroup(rg);
        return new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, output);
    }

    public static void main(final String[] args) {
        System.exit(new IlluminaBasecallsToSam().instanceMain(args));
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access args.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     *         to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        if (BARCODE_PARAMS != null) {
            LIBRARY_PARAMS = BARCODE_PARAMS;
        }

        final ArrayList<String> messages = new ArrayList<String>();

        readStructure = new ReadStructure(READ_STRUCTURE);
        if (!readStructure.barcodes.isEmpty()) {
            if (LIBRARY_PARAMS == null) {
                messages.add("BARCODE_PARAMS or LIBRARY_PARAMS is missing.  If READ_STRUCTURE contains a B (barcode)" +
                        " then either LIBRARY_PARAMS or BARCODE_PARAMS(deprecated) must be provided!");
            }
        }

        if (READ_GROUP_ID == null) {
            READ_GROUP_ID = RUN_BARCODE.substring(0, 5) + "." + LANE;
        }
        if (messages.size() == 0) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }

    /**
     * Given a read structure return the data types that need to be parsed for this run
     */
    private static IlluminaDataType[] getDataTypesFromReadStructure(final ReadStructure readStructure) {
        if (readStructure.barcodes.isEmpty()) {
            return DATA_TYPES_NO_BARCODE;
        } else {
            return DATA_TYPES_WITH_BARCODE;
        }
    }


}
