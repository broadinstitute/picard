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

import net.sf.picard.illumina.parser.*;
import net.sf.picard.util.FileChannelJDKBugWorkAround;
import net.sf.picard.util.TabbedTextFileWithHeaderParser;
import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;
import net.sf.samtools.util.*;

import java.io.File;
import java.util.*;

/**
 * @author jburke@broadinstitute.org
 */
public class IlluminaBasecallsToSam extends CommandLineProgram {

    private static final Log log = Log.getInstance(IlluminaBasecallsToSam.class);
    private static final boolean PRINT_TIMING = false;

    public static final String READ_STRUCTURE_DOC =
            "A description of the logical configuration of an Illumina Run, i.e. a description of what structure IlluminaBasecallsToSam "      +
            "assumes the  data to be in. It should consist of integer/character pairs describing the number of cycles and the type of those "  +
            "cycles (B for Barcode, T for Template, and S for skip).  E.g. If the input data consists of 80 base clusters and we provide a "   +
            "run configuration of \"36T8B8S30T\" then, before being converted to SAM records those bases will be split into 4 reads where "    +
            "read one consists of 36 cycles of template, read two consists of 8 cycles of barcode, read three will be an 8 base read of "      +
            "skipped cycles and read four is another 30 cycle template read.  The read consisting of skipped cycles would NOT be included "    +
            "in output SAM files.";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE =
        getStandardUsagePreamble() +  "Generate a SAM or BAM file from data in an Illumina basecalls output directory.\n";
    
    @Option(doc="The basecalls output directory. ", shortName="B")
    public File BASECALLS_DIR;
    @Option(doc="Lane number. ", shortName= StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;
    @Option(doc="The output SAM or BAM file. Format is determined by extension.",
            shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            mutex = {"BARCODE_PARAMS"})
    public File OUTPUT;
    
    @Option(doc = "Prefixed to read names.")
    public String RUN_BARCODE;
    @Option(doc="The name of the sequenced sample",
            shortName=StandardOptionDefinitions.SAMPLE_ALIAS_SHORT_NAME,
            mutex = {"BARCODE_PARAMS"})
    public String SAMPLE_ALIAS;
    @Option(doc="ID used to link RG header record with RG tag in SAM record.  " +
            "If these are unique in SAM files that get merged, merge performance is better.  " +
            "If not specified, READ_GROUP_ID = <first 5 chars of RUN_BARCODE>.<LANE> .",
        shortName = StandardOptionDefinitions.READ_GROUP_ID_SHORT_NAME, optional = true)
    public String READ_GROUP_ID;
    @Option(doc="The name of the sequenced library",
            shortName=StandardOptionDefinitions.LIBRARY_NAME_SHORT_NAME,
            optional=true,
            mutex = {"BARCODE_PARAMS"})
    public String LIBRARY_NAME;
    @Option(doc="The name of the sequencing center that produced the reads to fill in the RG.CN tag.", optional=true)
    public String SEQUENCING_CENTER = "BI";
    @Option(doc="The start date of the run.", optional=true)
    public Date RUN_START_DATE;
    @Option(doc="The name of the sequencing technology that produced the read.", optional=true)
    public String PLATFORM = "illumina";

    @Option(doc=READ_STRUCTURE_DOC, //see above
            optional = true, shortName="RS", mutex = {"BARCODE_CYCLE", "BARCODE_LENGTH"})
    public String READ_STRUCTURE;

    @Option(doc="Deprecated, use READ_STRUCTURE.  1-based cycle number of the start of the barcode.  If BARCODE or BARCODE_LENGTH are specified, this must be specified.",
            optional = true, shortName = "BARCODE_POSITION", mutex={"READ_STRUCTURE"})
    public Integer BARCODE_CYCLE;
    @Option(doc="Deprecated, use READ_STRUCTURE.  Length of the barcode.  If BARCODE_CYCLE or BARCODE is specified, this must be specified.",
            optional = true, mutex={"READ_STRUCTURE"})
    public Integer BARCODE_LENGTH;
    @Option(doc="Tab-separated file for creating all output BAMs for barcoded run with single IlluminaBasecallsToSam invocation.  " +
            "Columns are BARCODE, OUTPUT, SAMPLE_ALIAS, and LIBRARY_NAME.  Row with BARCODE=N is used to specify a file for no barcode match",
            mutex = {"OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME"})
    public File BARCODE_PARAMS;

    @Option(doc="Whether to mark the position of the adapter in the read")
    public boolean MARK_ADAPTER = true;
    @Option(doc = "Run this many TileProcessors in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.")
    public Integer NUM_PROCESSORS = 0;
    @Option(doc="If set, this is the first tile to be processed (for debugging).  Note that tiles are not processed in numerical order.",
    optional = true)
    public Integer FIRST_TILE;
    @Option(doc="If set, process no more than this many tiles (for debugging).", optional=true)
    public Integer TILE_LIMIT;

    @Option(doc="If true, call System.gc() periodically.  This is useful in cases in which the -Xmx value passed " +
            "is larger than the available memory.  Default: enable FORCE_GC if multi-threading.", optional = true)
    public Boolean FORCE_GC;
    @Option(doc="Configure SortingCollections in TileProcessors to store this many records before spilling to disk.  " + "" +
            "For an indexed run, each SortingCollection gets this value/number of indices.")
    public int MAX_READS_IN_RAM_PER_TILE = 1200000;

    private int recordsWritten = 0;
    private IlluminaBasecallsToSamConverter converter;
    private IlluminaDataProviderFactory factory;
    private IlluminaRunConfiguration runConfig;

    // key == barcode, value == corresponding SAMFileWriter
    // If not barcoded run, key == null.
    private final Map<String, SAMFileWriter> writersByBarcode = new HashMap<String, SAMFileWriter>();

    @Override
	protected int doWork() {
        
        log.info("POSITION  IS " + BARCODE_CYCLE);
        
        if (OUTPUT != null) {
            IoUtil.assertFileIsWritable(OUTPUT);
        }
        if (BARCODE_PARAMS != null) {
            IoUtil.assertFileIsReadable(BARCODE_PARAMS);
        }

        if (READ_STRUCTURE != null) { //Possibly barcoded, possibly not depending on configuration
            final IlluminaRunConfiguration rc = new IlluminaRunConfiguration(READ_STRUCTURE); //TODO: Next round this all gets prettier
            if(rc.numBarcodes > 0) {
                factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, rc,
                    IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF, IlluminaDataType.Barcodes);
            } else {
                factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, rc,
                    IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF);
            }
        } else if(BARCODE_CYCLE != null) { //Barcoded run whose configuration is determined (using the Barcode params and the qseqs)
            factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, BARCODE_CYCLE, BARCODE_LENGTH,
                IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF, IlluminaDataType.Barcodes);
        } else { //non-barcoded run whose configuration is determined (using available qseqs)
            factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, null,
                IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF);
        }

        //In the next round of changes only the first constructor above will be used and this will be replaced by the construction
        //of an IlluminaRunConfiguration object in the same manner as the first branch of the if else block above
        runConfig = factory.getRunConfig();
        log.info("RUN CONFIGURATION IS " + runConfig.toString());
        
        List<Integer> tiles = factory.getTiles();
        // Since the first non-fixed part of the read name is the tile number, without preceding zeroes,
        // and the output is sorted by read name, process the tiles in this order.
        Collections.sort(tiles, new Comparator<Integer>() {
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
        });
        if (FIRST_TILE != null) {
            int i;
            for (i = 0; i < tiles.size(); ++i) {
                if (tiles.get(i).intValue() == FIRST_TILE.intValue()) {
                    tiles = tiles.subList(i, tiles.size());
                    break;
                }
            }
            if (tiles.get(0).intValue() != FIRST_TILE.intValue()) {
                throw new PicardException("FIRST_TILE=" + FIRST_TILE +", but that tile was not found.");
            }
        }
        if (TILE_LIMIT != null && tiles.size() > TILE_LIMIT) {
            tiles = tiles.subList(0, TILE_LIMIT);
        }

        if (OUTPUT != null) {
            writersByBarcode.put(null, buildSamFileWriter(OUTPUT, SAMPLE_ALIAS, LIBRARY_NAME, null));
        } else {
            populateWritersByBarcode();
        }

        converter = new IlluminaBasecallsToSamConverter(RUN_BARCODE, READ_GROUP_ID, runConfig);

        // Process each tile separately, so the order of tile processing is determined by the
        // order of tiles list.  The SAMRecords produced will be cached and sorted a tile
        // at a time, rather than using SortingCollection to sort all of them.
        int numProcessors;
        if (NUM_PROCESSORS == 0) {
            numProcessors = Runtime.getRuntime().availableProcessors();
        }
        else if (NUM_PROCESSORS < 0) {
            numProcessors = Runtime.getRuntime().availableProcessors() + NUM_PROCESSORS;
        }
        else {
            numProcessors = NUM_PROCESSORS;
        }

        numProcessors = Math.min(numProcessors, tiles.size());

        if (FORCE_GC == null) {
            FORCE_GC = (numProcessors > 1);
        }

        if (numProcessors > 1) {
            // TODO: Eliminate this when switch to JDK 7
            FileChannelJDKBugWorkAround.doBugWorkAround();
            
            // Multiple processors -- multithread
            log.info("Creating " + numProcessors + " TileProcessors.");
            final LinkedList<TileProcessor> activeProcessors = new LinkedList<TileProcessor>();
            final Iterator<Integer> tileIterator = tiles.iterator();
            for (int i = 0; i < numProcessors; ++i) {
                final TileProcessor tileProcessor = new TileProcessor();
                tileProcessor.initialize(tileIterator.next());
                tileProcessor.startThread();
                activeProcessors.addLast(tileProcessor);
            }
            while (!activeProcessors.isEmpty()) {
                final TileProcessor earliestProcessor = activeProcessors.removeFirst();
                maybeGC();
                earliestProcessor.waitUntilDone();

                // Writing is done on main thread in order to produce consistent order.
                earliestProcessor.sortAndFlush();

                // Recycle the TileProcessor rather than creating a new one, because the sorter containing
                // in it holds a big array, and it probably helps memory situation not to release and reallocate that.
                if (tileIterator.hasNext()) {
                    earliestProcessor.initialize(tileIterator.next());
                    earliestProcessor.startThread();
                    activeProcessors.addLast(earliestProcessor);
                }
            }
        } else {
            // Only one processor, so single-thread
            log.info("Single TileProcessor mode.");
            final TileProcessor tileProcessor = new TileProcessor();
            for (final int tile : tiles) {
                tileProcessor.initialize(tile);
                tileProcessor.run();
                maybeGC();
                tileProcessor.sortAndFlush();
            }
        }
        for (final SAMFileWriter writer : writersByBarcode.values()) {
            writer.close();
        }

        log.info("Wrote " + recordsWritten + " read records total.");
        return 0;
    }

    private void populateWritersByBarcode() {
        final TabbedTextFileWithHeaderParser barcodeParamsParser = new TabbedTextFileWithHeaderParser(BARCODE_PARAMS);
        for (final String columnLabel : new String[]{"BARCODE", "OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME"}) {
            if (!barcodeParamsParser.hasColumn(columnLabel)) {
                throw new PicardException("BARCODE_PARAMS file " + BARCODE_PARAMS + " does not have column " +
                        columnLabel + ".");
            }
        }
        for (final TabbedTextFileWithHeaderParser.Row row : barcodeParamsParser) {
            final String barcode = row.getField("BARCODE");
            final String key = barcode.equals("N")? null: barcode;
            if (writersByBarcode.containsKey(key)) {
                throw new PicardException("Row for barcode " + barcode + " appears more than once in BARCODE_PARAMS file " +
                        BARCODE_PARAMS);
            }
            final SAMFileWriter writer = buildSamFileWriter(new File(row.getField("OUTPUT")),
                    row.getField("SAMPLE_ALIAS"), row.getField("LIBRARY_NAME"), barcode);
            writersByBarcode.put(key, writer);
        }
        if (writersByBarcode.isEmpty()) {
            throw new PicardException("BARCODE_PARAMS file " + BARCODE_PARAMS + " does have any data rows.");
        }
    }

    private void maybeGC() {
        if (FORCE_GC) {
            System.out.println("Before explicit GC, Runtime.totalMemory()=" + Runtime.getRuntime().totalMemory());
            System.gc();
            System.runFinalization();
            System.out.println("After explicit GC, Runtime.totalMemory()=" + Runtime.getRuntime().totalMemory());
        }
    }

    /**
     * Synchronized because called by multiple threads.  Might this be a source of contention?  This could
     * be called during sortAndFlush(), which is all on the main thread.
     */
    private synchronized void incrementRecordsWritten() {
        recordsWritten++;
        if (recordsWritten % 1000000 == 0) {
            log.info(recordsWritten + " read records written...");
        }

    }
    
    private SAMFileWriter buildSamFileWriter(final File output, final String sampleAlias, final String libraryName, final String barcode) {
        IoUtil.assertFileIsWritable(output);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(READ_GROUP_ID);
        rg.setSample(sampleAlias);
        String platformUnit = RUN_BARCODE + "." + LANE;
        if (barcode != null) platformUnit += ("." + barcode);
        rg.setAttribute("PU", platformUnit);

        if (libraryName != null) rg.setLibrary(libraryName);
        if (SEQUENCING_CENTER != null) rg.setAttribute("CN", SEQUENCING_CENTER);
        if (PLATFORM != null) rg.setAttribute("PL", PLATFORM);
        if (RUN_START_DATE != null) {
            final Iso8601Date date = new Iso8601Date(RUN_START_DATE);
            rg.setAttribute("DT", date.toString());
        }
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        header.addReadGroup(rg);
        return new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, output);
    }

    public static void main(final String[] argv) {
        System.exit(new IlluminaBasecallsToSam().instanceMain(argv));
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     *         to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> messages = new ArrayList<String>();

        //TODO: Check this!
        //Can have only a READ_STRUCTURE, a READ_STRUCTURE and a Barcode Params, or a BARCODE CYCLE, BARCODE LENGTH and BARCODE PARAMS
        if(READ_STRUCTURE == null) {
            if(BARCODE_CYCLE != null || BARCODE_LENGTH != null || BARCODE_PARAMS != null) {
                if(BARCODE_CYCLE == null) {
                    messages.add("BARCODE_CYCLE is missing.  Either each parameter BARCODE_CYCLE, BARCODE_LENGTH, and BARCODE_PARAMS or none must be specified.");
                } else {
                    if(BARCODE_CYCLE < 1)
                        messages.add("BARCODE_CYCLE must be >= 1");
                }
                if(BARCODE_LENGTH == null) {
                    messages.add("BARCODE_LENGTH is missing.  Either each parameter BARCODE_CYCLE, BARCODE_LENGTH, and BARCODE_PARAMS or none must be specified.");
                }
                if(BARCODE_PARAMS == null) {
                    messages.add("BARCODE_PARAMS is missing.  Either each parameter BARCODE_CYCLE, BARCODE_LENGTH, and BARCODE_PARAMS or none must be specified.");
                }
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

    private boolean isBarcoded() {
        return BARCODE_CYCLE != null;
    }

    /**
     * Accumulate the SAMRecord for a tile in RAM, so they they can be query_name sorted,
     * and written to a SAMWriter with presorted==true.
     */
    private static class SamRecordSorter {
        private final Comparator<SAMRecord> comparator = new SAMRecordQueryNameComparator();
        private final SortingCollection.Codec<SAMRecord> codec;
        private final SAMFileWriter writer;
        private final int maxRecordsInRam;
        private final List<File> tmpDirs;
        private SortingCollection<SAMRecord> records;
        private int numRecords = 0;

        /**
         * Object to capture some SAMRecords (a tile's worth), sort them in RAM, and write them
         * to a SAMFileWriter.
         * @param writer Where to write the records.
         * @param maxRecordsInRam Passed through to sorting collection
         */
        private SamRecordSorter(final int maxRecordsInRam, final SAMFileWriter writer, final List<File> tmpDirs) {
            this.maxRecordsInRam = maxRecordsInRam;
            this.tmpDirs = tmpDirs;
            codec = new BAMRecordCodec(writer.getFileHeader());
            this.writer = writer;
            createSortingCollection();
        }

        private void createSortingCollection() {
			if (records != null) records.cleanup();
            records = SortingCollection.newInstance(SAMRecord.class, codec, comparator, maxRecordsInRam, tmpDirs);
        }

        void addAlignment(final SAMRecord rec) {
            records.add(rec);
            ++numRecords;
        }

        void sortAndFlush() {
            final StopWatch writeWatch;
            if (PRINT_TIMING) {
                writeWatch = new StopWatch();
                writeWatch.start();
            }

            records.doneAdding();
            final PeekIterator<SAMRecord> it = new PeekIterator<SAMRecord>(records.iterator());
            while (it.hasNext()) {
                final SAMRecord rec = it.next();

                // PIC-330 Sometimes there are two reads with the same cluster coords, and thus
                // the same read name.  Discard both of them.  This code assumes that the two first of pairs
                // will come before the two second of pairs, so it isn't necessary to look ahead a different
                // distance for paired end.  It also assumes that for paired ends there will be duplicates
                // for both ends, so there is no need to be PE-aware.
                if (it.hasNext()) {
                    final SAMRecord lookAhead = it.peek();
                    if (!rec.getReadUnmappedFlag() || !lookAhead.getReadUnmappedFlag()) {
                        throw new IllegalStateException("Should not have mapped reads.");
                    }
                    if (comparator.compare(rec, lookAhead) == 0) {
                        it.next();
                        log.info("Skipping reads with identical read names: " + rec.getReadName());
                        continue;
                    }
                }

                writer.addAlignment(rec);
            }
            if (PRINT_TIMING) {
                writeWatch.stop();
                System.err.println("msec to sort and write: " + writeWatch.getElapsedTime());
            }

            // Reset for reuse
            createSortingCollection();
            numRecords = 0;
        }

        int size() {
            return numRecords;
        }

        SAMFileWriter getWriter() {
            return writer;
        }
    }

    /**
     * Converts a single tile's reads into SAMRecords.
     *
     * Note that instances of this class are designed to be recycled after a tile is processed,
     * because the SamRecordSorter holds a big array and it probably makes the GC situation better not to keep
     * harvesting and allocating that anew.
     *
     * Because TileProcessors are designed to run in parallel, there are several synchronization points to be
     * addressed:
     *
     * * The main thread waits until the TileProcessor has finished converting SAMRecords, after which
     * it writes the SAMRecords and possibly recycles the instance.
     *
     * * If a later TileProcessor does not have a good secondary base caller (typically because there are not enough
     * reads in that tile), then it asks its predecessor for its caller.  The predecessor may not be done training
     * its caller yet, so the later TileProcessor may have to wait until the predecessor is finished training.
     *
     * * When a TileProcessor T1 is recycled, the most recent good secondary base caller is pushed into the
     * next TileProcessor T2.  This may happen before or after T2 has trained its own caller.  T2 may already
     * have asked T1 for its caller.  Care must be taken to avoid replacing T2's own caller with T1's caller,
     * if T2's caller is good. 
     */
    private class TileProcessor implements Runnable {
        // Map from barcode string to SamRecordSorter.  barcode may be null.
        private final Map<String, SamRecordSorter> sorters = new HashMap<String, SamRecordSorter>();

        // Main thread waits on this monitor, which is notified when the tile has been processed.
        private final Object tileProcessedMonitor = new Object();
        // If this is true, main thread should not wait on tileProcessedMonitor.
        private boolean tileProcessed;

        private int tile;

        // When in multi-threaded mode, this is the thread on which the TileProcessor is running.
        private Thread thread;

        // Set to true if in multi-threaded mode and the thread threw an exception.  waitUntilDone() checks
        // this and throws an exception so that the program aborts itself.
        private boolean threadFailed;

        private SAMRecord [] records;

        private TileProcessor() {
            final int maxRecordsInRam = MAX_READS_IN_RAM_PER_TILE / writersByBarcode.size();
            for (final Map.Entry<String, SAMFileWriter> entry : writersByBarcode.entrySet()) {
                sorters.put(entry.getKey(), new SamRecordSorter(maxRecordsInRam, entry.getValue(), TMP_DIR));
            }
            records = new SAMRecord[converter.getNumRecordsPerCluster()];
        }

        // For constructing or recycling the instance.
        void initialize(final int tile) {
            this.tile = tile;
            this.thread = null;
            this.threadFailed = false;
            synchronized (this.tileProcessedMonitor) {
                this.tileProcessed = false;
            }
        }

        /**
         * Called by the main program in multithreaded mode on the earliest active TileProcessor.  Blocks until
         * the TileProcessor has finished reading all the reads and converting to SAMRecords.  If the thread
         * threw an exception, that would have been written to the log.  This method detects that an exception was
         * thrown and also throws.
         */
        void waitUntilDone() {
            synchronized (this.tileProcessedMonitor) {
                if (this.threadFailed) {
                    throw new PicardException("TileProcessor thread terminated with an exception");
                }
                if (this.tileProcessed) {
                    return;
                }
                try {
                    this.tileProcessedMonitor.wait();
                    if (this.threadFailed) {
                        throw new PicardException("TileProcessor thread terminated with an exception");
                    }
                } catch (InterruptedException e) {
                    throw new PicardException("Waiting for tile to finish.", e);
                }
            }
        }

        public void run() {
            try {
                final List<Integer> oneTile = Arrays.asList(tile);
                final IlluminaDataProvider dataProvider = factory.makeDataProvider(oneTile);
                processTile(dataProvider);
            } catch (Throwable e) {
                if (this.thread != null) {
                    // In multi-thread mode, log the error and set a flag so that the main thread also throws.
                    this.threadFailed = true;
                    log.error(e, "Exception in TileProcessor");
                } else {
                    throw new PicardException("Exception in TileProcessor", e);
                }
            } finally {
                // Notify the main thread even if this thread died with an exception.
                synchronized (this.tileProcessedMonitor) {
                    this.tileProcessed = true;
                    this.tileProcessedMonitor.notifyAll();
                }
            }
        }

        /**
         * Assign this instance to a new Thread and start it.
         */
        public void startThread() {
            this.thread = new Thread(this);
            this.thread.start();
        }

        /**
         * Read all the reads and convert to SAMRecords.
         * @param dataProvider
         */
        private void processTile(final IlluminaDataProvider dataProvider) {
            final StopWatch readWatch;
            final StopWatch sqCallWatch;
            if (PRINT_TIMING) {
                sqCallWatch = new StopWatch();
                readWatch = new StopWatch();
            } else {
                sqCallWatch = null;
                readWatch = null;
            }
            while (dataProvider.hasNext()) {
                if (PRINT_TIMING) readWatch.start();
                final ClusterData cluster = dataProvider.next();
                if (PRINT_TIMING) readWatch.stop();

                final SamRecordSorter sorter = sorters.get(cluster.getMatchedBarcode());
                if (sorter == null) {
                    throw new PicardException("Barcode encountered in that was not specified in BARCODE_PARAMS: " +
                    cluster.getMatchedBarcode());
                }

                final IlluminaRunConfiguration runConfig = dataProvider.getRunConfig();
                if(runConfig.numTemplates != 1 && runConfig.numTemplates != 2) {
                    throw new PicardException("Number of templates(" + runConfig.numTemplates + ") specified by run configuration was greater than 2");
                }

                final SAMFileHeader header = sorter.getWriter().getFileHeader();
                converter.createSamRecords(cluster, header, MARK_ADAPTER, records);

                for(final SAMRecord sam : records) {
                    sorter.addAlignment(sam);
                    incrementRecordsWritten();
                }
            }
            if (PRINT_TIMING) {
                System.err.println("msec to read and train: " + readWatch.getElapsedTime());
                System.err.println("msec to call SQ tag: " + sqCallWatch.getElapsedTime());
            }
        }

        /**
         * Push the SAMRecords into the appropriate output SAM file in query order.
         */
        void sortAndFlush() {
            for (final Map.Entry<String, SamRecordSorter> entry : sorters.entrySet()) {
                log.info("Writing " + entry.getValue().size() + " records for tile " + tile + " and barcode " + entry.getKey());
                entry.getValue().sortAndFlush();
            }
        }
    }
}
