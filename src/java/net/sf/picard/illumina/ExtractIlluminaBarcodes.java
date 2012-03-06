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
import net.sf.picard.util.IlluminaUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.TabbedTextFileWithHeaderParser;
import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricBase;
import net.sf.picard.metrics.MetricsFile;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.text.NumberFormat;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Determine the barcode for each read in an Illumina lane.
 * For each tile, a file is written to the basecalls directory of the form s_<lane>_<tile>_barcode.txt.
 * An output file contains a line for each read in the tile, aligned with the regular basecall output
 * The output file contains the following tab-separated columns:
 * - read subsequence at barcode position
 * - Y or N indicating if there was a barcode match
 * - matched barcode sequence (empty if read did not match one of the barcodes).  If there is no match
 *   but we're close to the threshold of calling it a match we output the barcode that would have been
 *   matched but in lower case
 *
 * @author jburke@broadinstitute.org
 */
public class ExtractIlluminaBarcodes extends CommandLineProgram {

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE =
        getStandardUsagePreamble() +  "Determine the barcode for each read in an Illumina lane.\n" +
                "For each tile, a file is written to the basecalls directory of the form s_<lane>_<tile>_barcode.txt." +
                "An output file contains a line for each read in the tile, aligned with the regular basecall output\n" +
                "The output file contains the following tab-separated columns: \n" +
                "    * read subsequence at barcode position\n" +
                "    * Y or N indicating if there was a barcode match\n" +
                "    * matched barcode sequence\n" +
                "Note that the order of specification of barcodes can cause arbitrary differences in output for poorly matching barcodes.\n\n";

    @Option(doc="The Illumina basecalls output directory. ", shortName="B")
    public File BASECALLS_DIR;

    @Option(doc="Where to write _barcode.txt files.  By default, these are written to BASECALLS_DIR.", optional = true)
    public File OUTPUT_DIR;

    @Option(doc="Lane number. ", shortName= StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Option(doc=IlluminaBasecallsToSam.READ_STRUCTURE_DOC, shortName="RS")
    public String READ_STRUCTURE;

    @Option(doc="Barcode sequence.  These must be unique, and all the same length.  This cannot be used with reads that " +
            "have more than one barcode; use BARCODE_FILE in that case. ", mutex = {"BARCODE_FILE"})
    public List<String> BARCODE = new ArrayList<String>();

    @Option(doc="Tab-delimited file of barcode sequences, barcode name and and optionally library name.  " +
            "Barcodes must be unique, and all the same length.  Column headers must be 'barcode_sequence1', " +
            "'barcode_sequence2' (optional), 'barcode_name', and 'library_name'.", mutex = {"BARCODE"})
    public File BARCODE_FILE;

    @Option(doc="Per-barcode and per-lane metrics written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME)
    public File METRICS_FILE;

    @Option(doc="Maximum mismatches for a barcode to be considered a match.")
    public int MAX_MISMATCHES = 1;

    @Option(doc="Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.")
    public int MIN_MISMATCH_DELTA = 1;

    @Option(doc="Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.")
    public int MAX_NO_CALLS = 2;

    @Option(shortName="GZIP", doc="Compress output s_l_t_barcode.txt files using gzip and append a .gz extension to the filenames.")
    public boolean COMPRESS_OUTPUTS = false;

    @Option(doc = "Run this many PerTileBarcodeExtractors in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.")
    public int NUM_PROCESSORS = 1;

    private final Log log = Log.getInstance(ExtractIlluminaBarcodes.class);

    private ReadStructure readStructure;
    private IlluminaDataProviderFactory factory;

    private final Map<String,BarcodeMetric> barcodeToMetrics = new LinkedHashMap<String,BarcodeMetric>();
    private BarcodeMetric noMatchMetric = null;

    private final NumberFormat tileNumberFormatter = NumberFormat.getNumberInstance();

    public ExtractIlluminaBarcodes() {
        tileNumberFormatter.setMinimumIntegerDigits(4);
        tileNumberFormatter.setGroupingUsed(false);
    }

    @Override
	protected int doWork() {

        IoUtil.assertDirectoryIsWritable(BASECALLS_DIR);
        IoUtil.assertFileIsWritable(METRICS_FILE);
        if (OUTPUT_DIR == null) {
            OUTPUT_DIR = BASECALLS_DIR;
        }
        IoUtil.assertDirectoryIsWritable(OUTPUT_DIR);

        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[readStructure.barcodes.length()];
        int index = 0;
        for (final ReadDescriptor d : readStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                noMatchBarcode[index++] = StringUtil.repeatCharNTimes('N', d.length);
            }
        }
        noMatchMetric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

        final int numProcessors;
        if (NUM_PROCESSORS == 0) {
            numProcessors = Runtime.getRuntime().availableProcessors();
        }
        else if (NUM_PROCESSORS < 0) {
            numProcessors = Runtime.getRuntime().availableProcessors() + NUM_PROCESSORS;
        }
        else {
            numProcessors = NUM_PROCESSORS;
        }

        log.info("Processing with " + numProcessors + " PerTileBarcodeExtractor(s).");
        final ExecutorService pool = Executors.newFixedThreadPool(numProcessors);

        final List<PerTileBarcodeExtractor> extractors = new ArrayList<PerTileBarcodeExtractor>(factory.getAvailableTiles().size());
        for (final int tile : factory.getAvailableTiles()) {

            final PerTileBarcodeExtractor extractor = new PerTileBarcodeExtractor(
                tile, factory.makeDataProvider(Arrays.asList(tile)), getBarcodeFile(tile));
            pool.submit(extractor);
            extractors.add(extractor);
        }
        pool.shutdown();
        try {
            // Wait a while for existing tasks to terminate
            if (!pool.awaitTermination(6, TimeUnit.HOURS)) {
                pool.shutdownNow(); // Cancel any still-executing tasks
                // Wait a while for tasks to respond to being cancelled
                if (!pool.awaitTermination(60, TimeUnit.SECONDS))
                    log.error("Pool did not terminate");
                return 1;
            }
        }
        catch (InterruptedException ie) {
            // (Re-)Cancel if current thread also interrupted
            pool.shutdownNow();
            return 2;
        }

        log.info("Processed " + extractors.size() + " tiles.");
        for (final PerTileBarcodeExtractor extractor : extractors) {
            for (final String key : barcodeToMetrics.keySet()) {
                barcodeToMetrics.get(key).merge(extractor.getMetrics().get(key));
            }
            noMatchMetric.merge(extractor.getNoMatchMetric());
            if (extractor.getException() != null) {
                log.error("Abandoning metrics calculation because one or more PerTileBarcodeExtractors failed.");
                return 4;
            }
        }


        // Finish metrics tallying.
        int totalReads = noMatchMetric.READS;
        int totalPfReads = noMatchMetric.PF_READS;
        int totalPfReadsAssigned = 0;
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            totalReads += barcodeMetric.READS;
            totalPfReads += barcodeMetric.PF_READS;
            totalPfReadsAssigned += barcodeMetric.PF_READS;
        }

        if (totalReads > 0) {
            noMatchMetric.PCT_MATCHES = noMatchMetric.READS/(double)totalReads;
            double bestPctOfAllBarcodeMatches = 0;
            for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                barcodeMetric.PCT_MATCHES = barcodeMetric.READS/(double)totalReads;
                if (barcodeMetric.PCT_MATCHES > bestPctOfAllBarcodeMatches) {
                    bestPctOfAllBarcodeMatches = barcodeMetric.PCT_MATCHES;
                }
            }
            if (bestPctOfAllBarcodeMatches > 0) {
                noMatchMetric.RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                        noMatchMetric.PCT_MATCHES/bestPctOfAllBarcodeMatches;
                for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                    barcodeMetric.RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                            barcodeMetric.PCT_MATCHES/bestPctOfAllBarcodeMatches;
                }
            }
        }

        if (totalPfReads > 0) {
            noMatchMetric.PF_PCT_MATCHES = noMatchMetric.PF_READS/(double)totalPfReads;
            double bestPctOfAllBarcodeMatches = 0;
            for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                barcodeMetric.PF_PCT_MATCHES = barcodeMetric.PF_READS/(double)totalPfReads;
                if (barcodeMetric.PF_PCT_MATCHES > bestPctOfAllBarcodeMatches) {
                    bestPctOfAllBarcodeMatches = barcodeMetric.PF_PCT_MATCHES;
                }
            }
            if (bestPctOfAllBarcodeMatches > 0) {
                noMatchMetric.PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                        noMatchMetric.PF_PCT_MATCHES/bestPctOfAllBarcodeMatches;
                for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                    barcodeMetric.PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                            barcodeMetric.PF_PCT_MATCHES/bestPctOfAllBarcodeMatches;
                }
            }
        }

        // Calculate the normalized matches
        if (totalPfReadsAssigned > 0) {
            final double mean = (double) totalPfReadsAssigned / (double) barcodeToMetrics.values().size();
            for (final BarcodeMetric m : barcodeToMetrics.values()) {
                m.PF_NORMALIZED_MATCHES = m.PF_READS  / mean;
            }
        }

        final MetricsFile<BarcodeMetric, Integer> metrics = getMetricsFile();
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            metrics.addMetric(barcodeMetric);
        }
        metrics.addMetric(noMatchMetric);
        metrics.write(METRICS_FILE);
        return 0;
    }

    /**
     * Create a barcode filename corresponding to the given tile qseq file.
     */
    private File getBarcodeFile(final int tile) {
        return new File(OUTPUT_DIR,
                        "s_" + LANE + "_" + tileNumberFormatter.format(tile) + "_barcode.txt" + (COMPRESS_OUTPUTS ? ".gz" : ""));
    }

    /**
     * Validate that POSITION >= 1, and that all BARCODEs are the same length and unique
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     *         to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> messages = new ArrayList<String>();

        readStructure = new ReadStructure(READ_STRUCTURE);
        factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, IlluminaDataType.BaseCalls, IlluminaDataType.PF);

        if (BARCODE_FILE != null) {
            parseBarcodeFile(messages);
        } else {
            final Set<String> barcodes = new HashSet<String>();
            for (final String barcode : BARCODE) {
                if (barcodes.contains(barcode)) {
                    messages.add("Barcode " + barcode + " specified more than once.");
                }
                barcodes.add(barcode);
                final BarcodeMetric metric = new BarcodeMetric(null, null, barcode, new String[]{barcode});
                barcodeToMetrics.put(barcode, metric);
            }
        }
        if (barcodeToMetrics.keySet().size() == 0) {
            messages.add("No barcodes have been specified.");
        }
        if (messages.size() == 0) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }

    public static void main(final String[] argv) {
        System.exit(new ExtractIlluminaBarcodes().instanceMain(argv));
    }

    private static final String BARCODE_SEQUENCE_COLUMN = "barcode_sequence";
    private static final String BARCODE_SEQUENCE_1_COLUMN = "barcode_sequence_1";
    private static final String BARCODE_NAME_COLUMN = "barcode_name";
    private static final String LIBRARY_NAME_COLUMN = "library_name";

    private void parseBarcodeFile(final ArrayList<String> messages) {
        final TabbedTextFileWithHeaderParser barcodesParser = new TabbedTextFileWithHeaderParser(BARCODE_FILE);
        final String sequenceColumn = barcodesParser.hasColumn(BARCODE_SEQUENCE_COLUMN)
                ? BARCODE_SEQUENCE_COLUMN : barcodesParser.hasColumn(BARCODE_SEQUENCE_1_COLUMN)
                ? BARCODE_SEQUENCE_1_COLUMN : null;
        if (sequenceColumn == null) {
            messages.add(BARCODE_FILE + " does not have " + BARCODE_SEQUENCE_COLUMN + " or " +
                    BARCODE_SEQUENCE_1_COLUMN + " column header");
            return;
        }
        final boolean hasBarcodeName = barcodesParser.hasColumn(BARCODE_NAME_COLUMN);
        final boolean hasLibraryName = barcodesParser.hasColumn(LIBRARY_NAME_COLUMN);
        final int numBarcodes = readStructure.barcodes.length();
        final Set<String> barcodes = new HashSet<String>();
        for (final TabbedTextFileWithHeaderParser.Row row : barcodesParser) {
            final String bcStrings[] = new String[numBarcodes];
            int barcodeNum = 1;
            for (final ReadDescriptor rd : readStructure.descriptors) {
                if (rd.type != ReadType.Barcode) continue;
                final String header = barcodeNum == 1 ? sequenceColumn : "barcode_sequence_" + String.valueOf(barcodeNum);
                bcStrings[barcodeNum-1] = row.getField(header);
                barcodeNum++;
            }
            final String bcStr = IlluminaUtil.barcodeSeqsToString(bcStrings);
            if (barcodes.contains(bcStr)) {
                messages.add("Barcode " + bcStr + " specified more than once in " + BARCODE_FILE);
            }
            barcodes.add(bcStr);
            final String barcodeName = (hasBarcodeName? row.getField(BARCODE_NAME_COLUMN): "");
            final String libraryName = (hasLibraryName? row.getField(LIBRARY_NAME_COLUMN): "");
            final BarcodeMetric metric = new BarcodeMetric(barcodeName, libraryName, bcStr, bcStrings);
            barcodeToMetrics.put(StringUtil.join("", bcStrings), metric);
        }
        barcodesParser.close();
    }


    /**
     * Metrics produced by the ExtractIlluminaBarcodes program that is used to parse data in
     * the basecalls directory and determine to which barcode each read should be assigned.
     */
    public static class BarcodeMetric extends MetricBase {
        /**
         * The barcode (from the set of expected barcodes) for which the following metrics apply.
         * Note that the "symbolic" barcode of NNNNNN is used to report metrics for all reads that
         * do not match a barcode.
         */
        public String BARCODE;
        public String BARCODE_NAME = "";
        public String LIBRARY_NAME = "";
        /** The total number of reads matching the barcode. */
        public int READS = 0;
        /** The number of PF reads matching this barcode (always less than or equal to READS). */
        public int PF_READS = 0;
        /** The number of all reads matching this barcode that matched with 0 errors or no-calls. */
        public int PERFECT_MATCHES = 0;
        /** The number of PF reads matching this barcode that matched with 0 errors or no-calls. */
        public int PF_PERFECT_MATCHES = 0;
        /** The number of all reads matching this barcode that matched with 1 error or no-call. */
        public int ONE_MISMATCH_MATCHES = 0;
        /** The number of PF reads matching this barcode that matched with 1 error or no-call. */
        public int PF_ONE_MISMATCH_MATCHES = 0;
        /** The percentage of all reads in the lane that matched to this barcode. */
        public double PCT_MATCHES = 0d;
        /**
         * The rate of all reads matching this barcode to all reads matching the most prevelant barcode. For the
         * most prevelant barcode this will be 1, for all others it will be less than 1.  One over the lowest
         * number in this column gives you the fold-difference in representation between barcodes.
         */
        public double RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT = 0d;
        /** The percentage of PF reads in the lane that matched to this barcode. */
        public double PF_PCT_MATCHES = 0d;

        /**
         * The rate of PF reads matching this barcode to PF reads matching the most prevelant barcode. For the
         * most prevelant barcode this will be 1, for all others it will be less than 1.  One over the lowest
         * number in this column gives you the fold-difference in representation of PF reads between barcodes.
         */
        public double PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT = 0d;

        /**
         * The "normalized" matches to each barcode. This is calculated as the number of pf reads matching
         * this barcode over the sum of all pf reads matching any barcode (excluding orphans). If all barcodes
         * are represented equally this will be 1.
         */
        public double PF_NORMALIZED_MATCHES;

        protected byte[][] barcodeBytes;

        public BarcodeMetric(final String barcodeName, final String libraryName,
                             final String barcodeDisplay, final String[] barcodeSeqs) {

            this.BARCODE = barcodeDisplay;
            this.BARCODE_NAME = barcodeName;
            this.LIBRARY_NAME = libraryName;
            this.barcodeBytes = new byte[barcodeSeqs.length][];
            for (int i = 0; i < barcodeSeqs.length; i++) {
                barcodeBytes[i] = net.sf.samtools.util.StringUtil.stringToBytes(barcodeSeqs[i]);
            }
        }

        /**
         * This ctor is necessary for when reading metrics from file
         */
        public BarcodeMetric() {
            barcodeBytes = null;
        }

        /**
         * Creates a copy of metric initialized with only non-accumulated and non-calculated values set
         */
        public static BarcodeMetric copy(final BarcodeMetric metric) {
            final BarcodeMetric result = new BarcodeMetric();
            result.BARCODE = metric.BARCODE;
            result.BARCODE_NAME = metric.BARCODE_NAME;
            result.LIBRARY_NAME = metric.LIBRARY_NAME;
            result.barcodeBytes = metric.barcodeBytes;
            return result;
        }


        /**
         * Adds the non-calculated 
         * @param metric
         */
        public void merge(final BarcodeMetric metric) {
            this.READS += metric.READS;
            this.PF_READS += metric.PF_READS;
            this.PERFECT_MATCHES += metric.PERFECT_MATCHES;
            this.PF_PERFECT_MATCHES += metric.PF_PERFECT_MATCHES;
            this.ONE_MISMATCH_MATCHES += metric.ONE_MISMATCH_MATCHES;
            this.PF_ONE_MISMATCH_MATCHES += metric.PF_ONE_MISMATCH_MATCHES;
        }

    }

    /**
     * Extracts barcodes and accumulates metrics for an entire tile.
     */
    private class PerTileBarcodeExtractor
            implements Runnable {

        private final int tile;
        private final IlluminaDataProvider provider;
        private final File barcodeFile;
        private final Map<String,BarcodeMetric> metrics;
        private final BarcodeMetric noMatch;
        private Exception exception = null;

        /**
         * Utility class to hang onto data about the best match for a given barcode
         */
        class BarcodeMatch {
            boolean matched;
            String barcode;
            int mismatches;
            int mismatchesToSecondBest;
        }


        /**
         * Constructor
         * @param tile              The number of the tile being processed; used for logging only.
         * @param provider          The IlluminaDataProvider for the tile.
         * @param barcodeFile       The file to write the barcodes to
         */
        public PerTileBarcodeExtractor(final int tile, final IlluminaDataProvider provider, final File barcodeFile) {
            this.tile = tile;
            this.provider = provider;
            this.barcodeFile = barcodeFile;
            this.metrics = new LinkedHashMap<String,BarcodeMetric>(barcodeToMetrics.size());
            for (final String key : barcodeToMetrics.keySet()) {
                this.metrics.put(key, BarcodeMetric.copy(barcodeToMetrics.get(key)));
            }
            this.noMatch = BarcodeMetric.copy(noMatchMetric);
        }


        // These methods return the results of the extraction
        public synchronized Map<String,BarcodeMetric> getMetrics() { return this.metrics; }
        public synchronized BarcodeMetric getNoMatchMetric() { return this.noMatch; }
        public synchronized Exception getException() { return this.exception; }

        /**
         * run method which extracts barcodes and accumulates metrics for an entire tile
         */
        synchronized public void run() {
            log.info("Extracting barcodes for tile " + tile);
            final int [] barcodeIndices = readStructure.barcodes.getIndices();
            final BufferedWriter writer = IoUtil.openFileForBufferedWriting(barcodeFile);
            try {
                final byte barcodeSubsequences[][] = new byte[barcodeIndices.length][];
                while (provider.hasNext()) {
                    // Extract the barcode from the cluster and write it to the file for the tile
                    final ClusterData cluster = provider.next();
                    for (int i = 0; i < barcodeIndices.length; i++) {
                        barcodeSubsequences[i] = cluster.getRead(barcodeIndices[i]).getBases();
                    }
                    final boolean passingFilter = cluster.isPf();
                    final BarcodeMatch match = findBestBarcodeAndUpdateMetrics(barcodeSubsequences, passingFilter, metrics, noMatchMetric);

                    final String yOrN = (match.matched ? "Y" : "N");

                    for (final byte[] bc : barcodeSubsequences) {
                        writer.write(StringUtil.bytesToString(bc) + "\t");
                    }
                    writer.write(yOrN + "\t" + match.barcode + "\t" + String.valueOf(match.mismatches) +
                                 "\t" + String.valueOf(match.mismatchesToSecondBest));
                    writer.newLine();
                }
                writer.close();
            }
            catch (Exception e) {
                log.error("Error processing tile " + this.tile, e);
                this.exception = e;
            }
        }

        /**
         * Find the best barcode match for the given read sequence, and accumulate metrics
         * @param readSubsequences portion of read containing barcode
         * @param passingFilter PF flag for the current read
         * @return perfect barcode string, if there was a match within tolerance, or null if not.
         */
        private BarcodeMatch findBestBarcodeAndUpdateMetrics(final byte readSubsequences[][],
                                                             final boolean passingFilter,
                                                             final Map<String, BarcodeMetric> metrics,
                                                             final BarcodeMetric noMatchBarcodeMetric) {
            BarcodeMetric bestBarcodeMetric = null;
            int totalBarcodeReadBases = 0;
            int numNoCalls = 0; // NoCalls are calculated for all the barcodes combined

            for (final byte[] bc : readSubsequences) {
                totalBarcodeReadBases += bc.length;
                for (final byte b : bc) if (SequenceUtil.isNoCall(b)) ++numNoCalls;
            }

            // PIC-506 When forcing all reads to match a single barcode, allow a read to match even if every
            // base is a mismatch.
            int numMismatchesInBestBarcode = totalBarcodeReadBases + 1;
            int numMismatchesInSecondBestBarcode = totalBarcodeReadBases + 1;


            for (final BarcodeMetric barcodeMetric : metrics.values()) {
                final int numMismatches = countMismatches(barcodeMetric.barcodeBytes, readSubsequences);
                if (numMismatches < numMismatchesInBestBarcode) {
                    if (bestBarcodeMetric != null) {
                        numMismatchesInSecondBestBarcode = numMismatchesInBestBarcode;
                    }
                    numMismatchesInBestBarcode = numMismatches;
                    bestBarcodeMetric = barcodeMetric;
                } else if (numMismatches < numMismatchesInSecondBestBarcode) {
                    numMismatchesInSecondBestBarcode = numMismatches;
                }
            }

            final boolean matched = bestBarcodeMetric != null &&
                    numNoCalls <= MAX_NO_CALLS &&
                    numMismatchesInBestBarcode <= MAX_MISMATCHES &&
                    numMismatchesInSecondBestBarcode - numMismatchesInBestBarcode >= MIN_MISMATCH_DELTA;

            final BarcodeMatch match = new BarcodeMatch();

            // If we have something that's not a "match" but matches one barcode
            // slightly, we output that matching barcode in lower case
            if (numNoCalls + numMismatchesInBestBarcode < totalBarcodeReadBases) {
                match.mismatches = numMismatchesInBestBarcode;
                match.mismatchesToSecondBest = numMismatchesInSecondBestBarcode;
                match.barcode = bestBarcodeMetric.BARCODE.toLowerCase().replaceAll(IlluminaUtil.BARCODE_DELIMITER, "");
            }
            else {
                match.mismatches = totalBarcodeReadBases;
                match.barcode = "";
            }

            if (matched) {
                ++bestBarcodeMetric.READS;
                if (passingFilter) {
                    ++bestBarcodeMetric.PF_READS;
                }
                if (numMismatchesInBestBarcode == 0) {
                    ++bestBarcodeMetric.PERFECT_MATCHES;
                    if (passingFilter) {
                        ++bestBarcodeMetric.PF_PERFECT_MATCHES;
                    }
                } else if (numMismatchesInBestBarcode == 1) {
                    ++bestBarcodeMetric.ONE_MISMATCH_MATCHES;
                    if (passingFilter) {
                        ++bestBarcodeMetric.PF_ONE_MISMATCH_MATCHES;
                    }
                }

                match.matched = true;
                match.barcode = bestBarcodeMetric.BARCODE.replaceAll(IlluminaUtil.BARCODE_DELIMITER, "");
            }
            else {
                ++noMatchBarcodeMetric.READS;
                if (passingFilter) {
                    ++noMatchBarcodeMetric.PF_READS;
                }
            }

            return match;
        }

        /**
         * Compare barcode sequence to bases from read
         * @return how many bases did not match
         */
        private int countMismatches(final byte[][] barcodeBytes, final byte[][] readSubsequence) {
            int numMismatches = 0;
            // Read sequence and barcode length may not be equal, so we just use the shorter of the two
            for (int j = 0; j < barcodeBytes.length; j++) {
                final int basesToCheck = Math.min(barcodeBytes[j].length, readSubsequence[j].length);
                for (int i = 0; i < basesToCheck; ++i) {
                    if (!SequenceUtil.isNoCall(readSubsequence[j][i]) && !SequenceUtil.basesEqual(barcodeBytes[j][i], readSubsequence[j][i])) {
                        ++numMismatches;
                    }
                }
            }
            return numMismatches;
        }


    }

}
