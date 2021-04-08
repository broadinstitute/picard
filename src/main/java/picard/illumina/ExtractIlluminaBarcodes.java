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
package picard.illumina;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.BaseCallingProgramGroup;
import picard.illumina.parser.*;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.*;

import java.io.BufferedWriter;
import java.io.File;
import java.text.NumberFormat;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Determine the barcode for each read in an Illumina lane.
 * For each tile, a file is written to the basecalls directory of the form s_<lane>_<tile>_barcode.txt.
 * An output file contains a line for each read in the tile, aligned with the regular basecall output
 * The output file contains the following tab-separated columns:
 * - read subsequence at barcode position
 * - Y or N indicating if there was a barcode match
 * - matched barcode sequence (empty if read did not match one of the barcodes).  If there is no match
 * but we're close to the threshold of calling it a match we output the barcode that would have been
 * matched but in lower case
 * - distance to best matching barcode, "mismatches" (*)
 * - distance to second-best matching barcode, "mismatchesToSecondBest" (*)
 *
 * NOTE (*): Due to an optimization the reported mismatches & mismatchesToSecondBest values may be inaccurate as long as
 * the conclusion (match vs. no-match) isn't affected. For example, reported mismatches and
 * mismatchesToSecondBest may be smaller than their true value if mismatches is truly larger than MAX_MISMATCHES.
 * Also, mismatchesToSecondBest might be smaller than its true value if its true value is greater than
 * mismatches + MIN_MISMATCH_DELTA.
 */
@CommandLineProgramProperties(

        summary = ExtractIlluminaBarcodes.USAGE_SUMMARY + ExtractIlluminaBarcodes.USAGE_DETAILS,
        oneLineSummary = ExtractIlluminaBarcodes.USAGE_SUMMARY,
        programGroup = BaseCallingProgramGroup.class
)
@DocumentedFeature
public class ExtractIlluminaBarcodes extends CommandLineProgram {

    /**
     * Column header for the first barcode sequence (preferred).
     */
    public static final String BARCODE_SEQUENCE_COLUMN = "barcode_sequence";
    /**
     * Column header for the first barcode sequence.
     */
    public static final String BARCODE_SEQUENCE_1_COLUMN = "barcode_sequence_1";
    /**
     * Column header for the barcode name.
     */
    public static final String BARCODE_NAME_COLUMN = "barcode_name";
    /**
     * Column header for the library name.
     */
    public static final String LIBRARY_NAME_COLUMN = "library_name";

    static final String USAGE_SUMMARY = "Tool determines the barcode for each read in an Illumina lane.  ";
    static final String USAGE_DETAILS = "<p>This tool determines the numbers of reads containing barcode-matching sequences and provides " +
            "statistics on the quality of these barcode matches.</p> " +
            "<p>Illumina sequences can contain at least two types of barcodes, sample and molecular (index).  Sample barcodes " +
            "(B in the read structure) are used to demultiplex pooled samples while index barcodes (M in the read structure) are used " +
            "to differentiate multiple reads of a template when carrying out paired-end sequencing.  Note that this tool only extracts " +
            "sample (B) and not molecular barcodes (M).</p>" +
            "" +
            "<p>Barcodes can be provided in the form of a list (BARCODE_FILE) or a string representing the barcode (BARCODE).  " +
            "The BARCODE_FILE contains multiple fields including '" + BARCODE_SEQUENCE_COLUMN + "' (or '" + BARCODE_SEQUENCE_1_COLUMN + "'), " +
            "'barcode_sequence_2' (optional), '" + BARCODE_NAME_COLUMN + "', and '" + LIBRARY_NAME_COLUMN + "'. " +
            "In contrast, the BARCODE argument is used for runs with reads containing a single " +
            "barcode (nonmultiplexed) and can be added directly as a string of text e.g. BARCODE=CAATAGCG.</p>" +
            "" +
            "<p>Data is output per lane/tile within the BaseCalls directory with the file name format of 's_{lane}_{tile}_barcode.txt'.  " +
            "These files contain the following tab-separated columns:" +
            "<ul> " +
            "<li>Read subsequence at barcode position</li>" +
            "<li>Y or N indicating if there was a barcode match</li>" +
            "<li>Matched barcode sequence (empty if read did not match one of the barcodes)</li>  " +
            "<li>The number of mismatches if there was a barcode match</li>  " +
            "<li>The number of mismatches to the second best barcode if there was a barcode match</li>  " +
            "</ul>" +
            "If there is no match but we're close to the threshold of calling it a match, we output the barcode that would have been " +
            "matched but in lower case.  Threshold values can be adjusted to accommodate barcode sequence mismatches from the reads." +
            "  The metrics file produced by the ExtractIlluminaBarcodes program indicates the number of matches (and mismatches)" +
            " between the barcode reads and the actual barcodes.  These metrics are provided both per-barcode and per lane and can be " +
            "found in the BaseCalls directory.</p>" +
            "<p>For poorly matching barcodes, the order of specification of barcodes can cause arbitrary output differences.</p>" +
            "" +
            "<h4>Usage example:</h4> " +
            "<pre>" +
            "java -jar picard.jar ExtractIlluminaBarcodes \\<br />" +
            "              BASECALLS_DIR=/BaseCalls/ \\<br />" +
            "              LANE=1 \\<br />" +
            "          READ_STRUCTURE=25T8B25T \\<br />" +
            "              BARCODE_FILE=barcodes.txt \\<br />" +
            "              METRICS_FILE=metrics_output.txt " +
            "</pre>" +
            "" +
            "Please see the ExtractIlluminaBarcodes.BarcodeMetric " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#ExtractIlluminaBarcodes.BarcodeMetric'>definitions</a> " +
            "for a complete description of the metrics produced by this tool.</p>" +
            "" +
            "<hr />";

    // The following attributes define the command-line arguments

    @Argument(doc = "The Illumina basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(doc = "Where to write _barcode.txt files.  By default, these are written to BASECALLS_DIR.", optional = true)
    public File OUTPUT_DIR;

    @Argument(doc = "Lane number. ", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Argument(doc = ReadStructure.PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Argument(doc = "Barcode sequence.  These must be unique, and all the same length.  This cannot be used with reads that " +
            "have more than one barcode; use BARCODE_FILE in that case. ", mutex = {"BARCODE_FILE"})
    public List<String> BARCODE = new ArrayList<>();

    @Argument(doc = "Tab-delimited file of barcode sequences, barcode name and, optionally, library name.  " +
            "Barcodes must be unique and all the same length.  Column headers must be '" + BARCODE_SEQUENCE_COLUMN + "' (or '" + BARCODE_SEQUENCE_1_COLUMN + "'), " +
            "'barcode_sequence_2' (optional), '" + BARCODE_NAME_COLUMN + "', and '" + LIBRARY_NAME_COLUMN + "'.", mutex = {"BARCODE"})
    public File BARCODE_FILE;

    @Argument(doc = "Per-barcode and per-lane metrics written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME)
    public File METRICS_FILE;

    @Argument(doc = "Maximum mismatches for a barcode to be considered a match.")
    public int MAX_MISMATCHES = 1;

    @Argument(doc = "Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.")
    public int MIN_MISMATCH_DELTA = 1;

    @Argument(doc = "Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.")
    public int MAX_NO_CALLS = 2;

    @Argument(shortName = "Q", doc = "Minimum base quality. Any barcode bases falling below this quality will be considered a mismatch even if the bases match.")
    public int MINIMUM_BASE_QUALITY = 0;

    @Argument(doc = "The minimum quality (after transforming 0s to 1s) expected from reads.  If qualities are lower than this value, an error is thrown." +
            "The default of 2 is what the Illumina's spec describes as the minimum, but in practice the value has been observed lower.")
    public int MINIMUM_QUALITY = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

    @Argument(shortName = "GZIP", doc = "Compress output s_l_t_barcode.txt files using gzip and append a .gz extension to the file names.")
    public boolean COMPRESS_OUTPUTS = false;

    @Argument(doc = "Run this many PerTileBarcodeExtractors in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.")
    public int NUM_PROCESSORS = 1;

    @Argument(doc = "The distance metric that should be used to compare the barcode-reads and the provided barcodes for finding the best and second-best assignments.")
    public DistanceMetric DISTANCE_MODE = DistanceMetric.HAMMING;

    private static final Log LOG = Log.getInstance(ExtractIlluminaBarcodes.class);

    /**
     * The read structure of the actual Illumina Run, i.e. the readStructure of the input data
     */
    private ReadStructure readStructure;

    private IlluminaDataProviderFactory factory;

    private final Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<>();
    private final ConcurrentHashMap<String, PerTileBarcodeExtractor.BarcodeMatch> barcodeLookupMap = new ConcurrentHashMap<>();

    private final NumberFormat tileNumberFormatter = NumberFormat.getNumberInstance();
    private BclQualityEvaluationStrategy bclQualityEvaluationStrategy;

    public ExtractIlluminaBarcodes() {
        tileNumberFormatter.setMinimumIntegerDigits(4);
        tileNumberFormatter.setGroupingUsed(false);
    }

    @Override
    protected int doWork() {

        IOUtil.assertFileIsWritable(METRICS_FILE);
        if (OUTPUT_DIR == null) {
            OUTPUT_DIR = BASECALLS_DIR;
        }
        IOUtil.assertDirectoryIsWritable(OUTPUT_DIR);

        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[readStructure.sampleBarcodes.length()];
        int index = 0;
        for (final ReadDescriptor d : readStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                noMatchBarcode[index++] = StringUtil.repeatCharNTimes('N', d.length);
            }
        }

        final BarcodeMetric noMatchMetric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

        final int numProcessors;
        if (NUM_PROCESSORS == 0) {
            numProcessors = Runtime.getRuntime().availableProcessors();
        } else if (NUM_PROCESSORS < 0) {
            numProcessors = Runtime.getRuntime().availableProcessors() + NUM_PROCESSORS;
        } else {
            numProcessors = NUM_PROCESSORS;
        }

        LOG.info("Processing with " + numProcessors + " PerTileBarcodeExtractor(s).");
        final ThreadPoolExecutorWithExceptions pool = new ThreadPoolExecutorWithExceptions(numProcessors);

        final List<PerTileBarcodeExtractor> extractors = new ArrayList<>(factory.getAvailableTiles().size());
        // TODO: This is terribly inefficient; we're opening a huge number of files via the extractor constructor and we never close them.
        for (final int tile : factory.getAvailableTiles()) {
            final PerTileBarcodeExtractor extractor = new PerTileBarcodeExtractor(
                    tile,
                    getBarcodeFile(tile),
                    barcodeToMetrics,
                    barcodeLookupMap,
                    noMatchMetric,
                    factory,
                    MINIMUM_BASE_QUALITY,
                    MAX_NO_CALLS,
                    MAX_MISMATCHES,
                    MIN_MISMATCH_DELTA,
                    DISTANCE_MODE
            );
            extractors.add(extractor);
        }

        for (final PerTileBarcodeExtractor extractor : extractors) {
            pool.submit(extractor);
        }
        pool.shutdown();
        ThreadPoolExecutorUtil.awaitThreadPoolTermination("Per tile extractor executor", pool, Duration.ofMinutes(5));

        if (pool.hasError()) {
            throw new PicardException("Exceptions in tile processing. There were " + pool.shutdownNow().size()
                    + " tasks that were still running or queued and have been cancelled. Errors: " + pool.exception.toString());
        }

        LOG.info("Processed " + extractors.size() + " tiles.");
        for (final PerTileBarcodeExtractor extractor : extractors) {
            for (final String key : barcodeToMetrics.keySet()) {
                barcodeToMetrics.get(key).merge(extractor.getMetrics().get(key));
            }
            noMatchMetric.merge(extractor.getNoMatchMetric());
            if (extractor.getException() != null) {
                LOG.error("Abandoning metrics calculation because one or more PerTileBarcodeExtractors failed.");
                return 4;
            }
        }

        // Finish metrics tallying.
        finalizeMetrics(barcodeToMetrics, noMatchMetric);

        // Warn about minimum qualities and assert that we've achieved the minimum.
        for (final Map.Entry<Byte, Integer> entry : bclQualityEvaluationStrategy.getPoorQualityFrequencies().entrySet()) {
            LOG.warn(String.format("Observed low quality of %s %s times.", entry.getKey(), entry.getValue()));
        }
        bclQualityEvaluationStrategy.assertMinimumQualities();

        final MetricsFile<BarcodeMetric, Integer> metrics = getMetricsFile();
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            metrics.addMetric(barcodeMetric);
        }
        metrics.addMetric(noMatchMetric);
        metrics.write(METRICS_FILE);
        return 0;
    }

    public static void finalizeMetrics(final Map<String, BarcodeMetric> barcodeToMetrics,
                                       final BarcodeMetric noMatchMetric) {
        // Finish metrics tallying.
        long totalReads = noMatchMetric.READS;
        long totalPfReads = noMatchMetric.PF_READS;
        long totalPfReadsAssigned = 0;
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            totalReads += barcodeMetric.READS;
            totalPfReads += barcodeMetric.PF_READS;
            totalPfReadsAssigned += barcodeMetric.PF_READS;
        }

        if (totalReads > 0) {
            noMatchMetric.PCT_MATCHES = noMatchMetric.READS / (double) totalReads;
            double bestPctOfAllBarcodeMatches = 0;
            for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                barcodeMetric.PCT_MATCHES = barcodeMetric.READS / (double) totalReads;
                if (barcodeMetric.PCT_MATCHES > bestPctOfAllBarcodeMatches) {
                    bestPctOfAllBarcodeMatches = barcodeMetric.PCT_MATCHES;
                }
            }
            if (bestPctOfAllBarcodeMatches > 0) {
                noMatchMetric.RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                        noMatchMetric.PCT_MATCHES / bestPctOfAllBarcodeMatches;
                for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                    barcodeMetric.RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                            barcodeMetric.PCT_MATCHES / bestPctOfAllBarcodeMatches;
                }
            }
        }

        if (totalPfReads > 0) {
            noMatchMetric.PF_PCT_MATCHES = noMatchMetric.PF_READS / (double) totalPfReads;
            double bestPctOfAllBarcodeMatches = 0;
            for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                barcodeMetric.PF_PCT_MATCHES = barcodeMetric.PF_READS / (double) totalPfReads;
                if (barcodeMetric.PF_PCT_MATCHES > bestPctOfAllBarcodeMatches) {
                    bestPctOfAllBarcodeMatches = barcodeMetric.PF_PCT_MATCHES;
                }
            }
            if (bestPctOfAllBarcodeMatches > 0) {
                noMatchMetric.PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                        noMatchMetric.PF_PCT_MATCHES / bestPctOfAllBarcodeMatches;
                for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                    barcodeMetric.PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                            barcodeMetric.PF_PCT_MATCHES / bestPctOfAllBarcodeMatches;
                }
            }
        }

        // Calculate the normalized matches
        if (totalPfReadsAssigned > 0) {
            final double mean = (double) totalPfReadsAssigned / (double) barcodeToMetrics.values().size();
            for (final BarcodeMetric m : barcodeToMetrics.values()) {
                m.PF_NORMALIZED_MATCHES = m.PF_READS / mean;
            }
        }
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
     * to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> messages = new ArrayList<>();

        this.bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);

        /**
         * In extract illumina barcodes we NEVER want to look at the template reads nor the molecular barcodes, therefore replace them with
         * skips because IlluminaDataProvider and its factory will neither open these nor produce ClusterData with the template reads in them,
         * thus reducing the file IO and value copying done by the data provider
         */
        readStructure = new ReadStructure(READ_STRUCTURE.replaceAll("T|M", "S"));
        final Set<IlluminaDataType> datatypes = (MINIMUM_BASE_QUALITY > 0) ?
                new HashSet<>(Arrays.asList(IlluminaDataType.BaseCalls, IlluminaDataType.PF, IlluminaDataType.QualityScores)) :
                new HashSet<>(Arrays.asList(IlluminaDataType.BaseCalls, IlluminaDataType.PF));
        factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, bclQualityEvaluationStrategy, datatypes);

        if (BARCODE_FILE != null) {
            parseBarcodeFile(messages);
        } else {
            final Set<String> barcodes = new HashSet<>();
            for (final String barcode : BARCODE) {
                if (barcodes.contains(barcode)) {
                    messages.add("Barcode " + barcode + " specified more than once.");
                }
                barcodes.add(barcode);
                final BarcodeMetric metric = new BarcodeMetric(null, null, barcode, new String[]{barcode});
                barcodeToMetrics.put(barcode, metric);
            }
        }
        if (barcodeToMetrics.keySet().isEmpty()) {
            messages.add("No barcodes have been specified.");
        }
        if (messages.isEmpty()) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }

    private void parseBarcodeFile(final ArrayList<String> messages) {
        try (final TabbedTextFileWithHeaderParser barcodesParser = new TabbedTextFileWithHeaderParser(BARCODE_FILE)) {
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
            final int numBarcodes = readStructure.sampleBarcodes.length();
            final Set<String> barcodes = new HashSet<>();
            for (final TabbedTextFileWithHeaderParser.Row row : barcodesParser) {
                final String[] bcStrings = new String[numBarcodes];
                int barcodeNum = 0;
                for (final ReadDescriptor rd : readStructure.descriptors) {
                    if (rd.type != ReadType.Barcode) {
                        continue;
                    }
                    final String header = barcodeNum == 0 ? sequenceColumn : "barcode_sequence_" + (1 + barcodeNum);
                    final String field = row.getField(header);
                    if (field == null) {
                        messages.add(String.format("Null barcode in column %s of row: %s", header, row.getCurrentLine()));
                        bcStrings[barcodeNum] = "";
                    } else {
                        bcStrings[barcodeNum] = field;
                    }
                    ++barcodeNum;
                }
                final String bcStr = IlluminaUtil.barcodeSeqsToString(bcStrings);
                if (barcodes.contains(bcStr)) {
                    messages.add("Barcode " + bcStr + " specified more than once in " + BARCODE_FILE);
                }
                barcodes.add(bcStr);
                final String barcodeName = (hasBarcodeName ? row.getField(BARCODE_NAME_COLUMN) : "");
                final String libraryName = (hasLibraryName ? row.getField(LIBRARY_NAME_COLUMN) : "");
                final BarcodeMetric metric = new BarcodeMetric(barcodeName, libraryName, bcStr, bcStrings);
                barcodeToMetrics.put(StringUtil.join("", bcStrings), metric);
            }
        }
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

        public String BARCODE_WITHOUT_DELIMITER;
        /**
         * The barcode name.
         */
        public String BARCODE_NAME = "";
        /**
         * The name of the library
         */
        public String LIBRARY_NAME = "";
        /**
         * The total number of reads matching the barcode.
         */
        public long READS = 0;
        /**
         * The number of PF reads matching this barcode (always less than or equal to READS).
         */
        public long PF_READS = 0;
        /**
         * The number of all reads matching this barcode that matched with 0 errors or no-calls.
         */
        public long PERFECT_MATCHES = 0;
        /**
         * The number of PF reads matching this barcode that matched with 0 errors or no-calls.
         */
        public long PF_PERFECT_MATCHES = 0;
        /**
         * The number of all reads matching this barcode that matched with 1 error or no-call.
         */
        public long ONE_MISMATCH_MATCHES = 0;
        /**
         * The number of PF reads matching this barcode that matched with 1 error or no-call.
         */
        public long PF_ONE_MISMATCH_MATCHES = 0;
        /**
         * The fraction of all reads in the lane that matched to this barcode.
         */
        public double PCT_MATCHES = 0d;
        /**
         * The rate of all reads matching this barcode to all reads matching the most prevelant barcode. For the
         * most prevelant barcode this will be 1, for all others it will be less than 1 (except for the possible
         * exception of when there are more orphan reads than for any other barcode, in which case the value
         * may be arbitrarily large).  One over the lowest number in this column gives you the fold-difference
         * in representation between barcodes.
         */
        public double RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT = 0d;
        /**
         * The fraction of PF reads in the lane that matched to this barcode.
         */
        public double PF_PCT_MATCHES = 0d;

        /**
         * The rate of PF reads matching this barcode to PF reads matching the most prevelant barcode. For the
         * most prevelant barcode this will be 1, for all others it will be less than 1 (except for the possible
         * exception of when there are more orphan reads than for any other barcode, in which case the value
         * may be arbitrarily large).  One over the lowest number in this column gives you the fold-difference
         * in representation of PF reads between barcodes.
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
            this.BARCODE_WITHOUT_DELIMITER = barcodeDisplay.replaceAll(IlluminaUtil.BARCODE_DELIMITER, "");
            this.BARCODE_NAME = barcodeName;
            this.LIBRARY_NAME = libraryName;
            this.barcodeBytes = new byte[barcodeSeqs.length][];
            for (int i = 0; i < barcodeSeqs.length; i++) {
                barcodeBytes[i] = htsjdk.samtools.util.StringUtil.stringToBytes(barcodeSeqs[i]);
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
            result.BARCODE_WITHOUT_DELIMITER = metric.BARCODE_WITHOUT_DELIMITER;
            result.BARCODE_NAME = metric.BARCODE_NAME;
            result.LIBRARY_NAME = metric.LIBRARY_NAME;
            result.barcodeBytes = metric.barcodeBytes;
            return result;
        }

        /**
         * Adds the non-calculated
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
    public static class PerTileBarcodeExtractor implements Runnable {
        private final int tile;
        private final File barcodeFile;
        private final Map<String, BarcodeMetric> metrics;
        private final BarcodeMetric noMatch;
        private Exception exception = null;
        private final boolean usingQualityScores;
        private BaseIlluminaDataProvider provider = null;
        private final ReadStructure outputReadStructure;
        private final int maxNoCalls, maxMismatches, minMismatchDelta, minimumBaseQuality;
        private final IlluminaDataProviderFactory factory;
        private final DistanceMetric distanceMode;
        private final ConcurrentHashMap<String, BarcodeMatch> barcodeLookupMap;
        private final static int maxLookupSize = 100000;

        /**
         * Utility class to hang onto data about the best match for a given barcode
         */
        public static class BarcodeMatch {
            boolean matched;
            String barcode;
            int mismatches;
            int mismatchesToSecondBest;

            public boolean isMatched() {
                return matched;
            }

            public String getBarcode() {
                return barcode;
            }
        }

        /**
         * Constructor
         *
         * @param tile             The number of the tile being processed; used for logging only.
         * @param barcodeFile      The file to write the barcodes to
         * @param noMatchMetric    A "template" metric that is cloned and the clone is stored internally for accumulating data
         * @param barcodeToMetrics A "template" metric map whose metrics are cloned, and the clones are stored internally for accumulating data
         */
        public PerTileBarcodeExtractor(
                final int tile,
                final File barcodeFile,
                final Map<String, BarcodeMetric> barcodeToMetrics,
                final ConcurrentHashMap<String, BarcodeMatch> barcodeLookupMap,
                final BarcodeMetric noMatchMetric,
                final IlluminaDataProviderFactory factory,
                final int minimumBaseQuality,
                final int maxNoCalls,
                final int maxMismatches,
                final int minMismatchDelta,
                final DistanceMetric distanceMode
        ) {
            this.tile = tile;
            this.barcodeFile = barcodeFile;
            this.usingQualityScores = minimumBaseQuality > 0;
            this.maxNoCalls = maxNoCalls;
            this.maxMismatches = maxMismatches;
            this.minMismatchDelta = minMismatchDelta;
            this.minimumBaseQuality = minimumBaseQuality;
            this.metrics = new LinkedHashMap<>(barcodeToMetrics.size());
            for (final String key : barcodeToMetrics.keySet()) {
                this.metrics.put(key, BarcodeMetric.copy(barcodeToMetrics.get(key)));
            }
            this.barcodeLookupMap = barcodeLookupMap;
            this.noMatch = BarcodeMetric.copy(noMatchMetric);
            this.outputReadStructure = factory.getOutputReadStructure();
            this.distanceMode = distanceMode;
            this.factory = factory;
        }

        // These methods return the results of the extraction
        public synchronized Map<String, BarcodeMetric> getMetrics() {
            return this.metrics;
        }

        public synchronized BarcodeMetric getNoMatchMetric() {
            return this.noMatch;
        }

        public synchronized Exception getException() {
            return this.exception;
        }

        /**
         * run method which extracts barcodes and accumulates metrics for an entire tile
         */
        public synchronized void run() {
            try {
                //delayed instantiation for new provider
                if (this.provider == null) {
                    this.provider = factory.makeDataProvider(tile);
                }
                LOG.info("Extracting barcodes for tile " + tile);

                // Sometimes makeDataProvider takes a while waiting for slow file IO, for each tile the needed set of files
                // is non-overlapping sets of files so make the  data providers in the individual threads for PerTileBarcodeExtractors
                // so they are not all waiting for each others file operations

                // Most likely we have SKIPS in our read structure since we replace all template reads with skips in the input data structure
                // (see customCommnandLineValidation), therefore we must use the outputReadStructure to index into the output cluster data
                final int[] barcodeIndices = outputReadStructure.sampleBarcodes.getIndices();
                final BufferedWriter writer = IOUtil.openFileForBufferedWriting(barcodeFile);
                final byte[][] barcodeSubsequences = new byte[barcodeIndices.length][];
                final byte[][] qualityScores = usingQualityScores ? new byte[barcodeIndices.length][] : null;
                while (provider.hasNext()) {
                    // Extract the barcode from the cluster and write it to the file for the tile
                    final ClusterData cluster = provider.next();

                    for (int i = 0; i < barcodeIndices.length; i++) {
                        barcodeSubsequences[i] = cluster.getRead(barcodeIndices[i]).getBases();
                        if (usingQualityScores) {
                            qualityScores[i] = cluster.getRead(barcodeIndices[i]).getQualities();
                        }
                    }
                    final boolean passingFilter = cluster.isPf();
                    final BarcodeMatch match = findBestBarcode(barcodeSubsequences, qualityScores,
                            metrics, maxNoCalls, maxMismatches,
                            minMismatchDelta, minimumBaseQuality);
                    updateMetrics(match, passingFilter, metrics, noMatch);

                    final String yOrN = (match.matched ? "Y" : "N");

                    for (final byte[] bc : barcodeSubsequences) {
                        writer.write(StringUtil.bytesToString(bc));
                    }
                    writer.write("\t" + yOrN + "\t" + match.barcode + "\t" + String.valueOf(match.mismatches) +
                            "\t" + String.valueOf(match.mismatchesToSecondBest));
                    writer.newLine();
                }
                writer.close();
            } catch (final Exception e) {
                LOG.error(e, "Error processing tile ", this.tile);
                this.exception = e;
            } finally {
                CloserUtil.close(provider);
                provider = null;
            }
        }

        private static boolean ensureLookupMinimumValue(final byte[][] qualityScores, final int minimumBaseQuality) {
            if (qualityScores != null) {
                for (final byte[] qs : qualityScores) {
                    for (final byte q : qs) {
                        if (q < minimumBaseQuality) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        /**
         * Find the best barcode match for the given read sequence, and accumulate metrics
         *
         * NOTE: the returned BarcodeMatch object will contain mismatches mismatchesToSecondBest values that may be
         * inaccurate as long as the conclusion match/no-match isn't affected. for example, mismatches and mismatchesToSecondBest
         * may be smaller than their true value if mismatches is truly larger than maxMismatches.
         * Also, mismatchesToSecondBest might be smaller than its true value if its true value is greater than
         * mismatches + minMismatchDelta. This is due to an optimization which allows the distance calculation to stop once
         * the conclusion (Match or no-Match) can be reached.
         *
         * @param readSubsequences portion of read containing barcode
         * @return perfect barcode string, if there was a match within tolerance, or null if not.
         */
        private BarcodeMatch findBestBarcode(final byte[][] readSubsequences,
                                            final byte[][] qualityScores,
                                            final Map<String, BarcodeMetric> metrics,
                                            final int maxNoCalls,
                                            final int maxMismatches,
                                            final int minMismatchDelta,
                                            final int minimumBaseQuality) {
            final boolean canUseLookupTable = ensureLookupMinimumValue(qualityScores, minimumBaseQuality);
            final BarcodeMatch match;
            final String barcodesAsString = IlluminaUtil.barcodeSeqsToString(readSubsequences);

            // this implementation is optimized for barcodeLookupMap being a ConcurrentHashMap for which this
            // pattern is faster than using computeIfAbsent (or rather, it locks the map
            // for a shorter time, allowing other threads to access it).

            // Also, a ConcurrentHashMap was used rather than a Cache since a high-performance, thread-safe
            // Cache was not found, hence the "poor man's cache" of using the first maxLookupSize distinct
            // barcode reads.

            if (canUseLookupTable && barcodeLookupMap.containsKey(barcodesAsString)) {
                match = barcodeLookupMap.get(barcodesAsString);
            } else {
                match = calculateBarcodeMatch(readSubsequences, qualityScores, metrics, maxNoCalls,
                        maxMismatches, minMismatchDelta,
                        minimumBaseQuality, distanceMode);

                if (canUseLookupTable && barcodeLookupMap.size() < maxLookupSize) {
                    barcodeLookupMap.put(barcodesAsString, match);
                }
            }

            return match;
        }

        static BarcodeMatch calculateBarcodeMatch(final byte[][] readSubsequences,
                                                  final byte[][] qualityScores,
                                                  final Map<String, BarcodeMetric> metrics,
                                                  final int maxNoCalls, final int maxMismatches,
                                                  final int minMismatchDelta, final int minimumBaseQuality,
                                                  final DistanceMetric distanceMode) {
            final BarcodeMatch match;
            BarcodeMetric bestBarcodeMetric = null;
            match = new BarcodeMatch();

            int totalBarcodeReadBases = 0;
            int numNoCalls = 0; // NoCalls are calculated for all the barcodes combined

            for (final byte[] bc : readSubsequences) {
                totalBarcodeReadBases += bc.length;
                for (final byte b : bc) {
                    if (SequenceUtil.isNoCall(b)) {
                        ++numNoCalls;
                    }
                }
            }

            // PIC-506 When forcing all reads to match a single barcode, allow a read to match even if every
            // base is a mismatch.
            int numMismatchesInBestBarcode = totalBarcodeReadBases + 1;
            int numMismatchesInSecondBestBarcode = totalBarcodeReadBases + 1;

            for (final BarcodeMetric barcodeMetric : metrics.values()) {
                // need to add maxMismatches + minMismatchDelta together since the result might get used as numMismatchesInSecondBestBarcode
                final BarcodeEditDistanceQuery barcodeEditDistanceQuery = new BarcodeEditDistanceQuery(barcodeMetric.barcodeBytes, readSubsequences, qualityScores,
                        minimumBaseQuality, Math.min(maxMismatches, numMismatchesInBestBarcode) + minMismatchDelta);
                final int numMismatches = distanceMode.distance(barcodeEditDistanceQuery);

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

            match.matched = bestBarcodeMetric != null &&
                    numNoCalls <= maxNoCalls &&
                    numMismatchesInBestBarcode <= maxMismatches &&
                    numMismatchesInSecondBestBarcode - numMismatchesInBestBarcode >= minMismatchDelta;

            // If we have something that's not a "match" but matches one barcode
            // slightly, we output that matching barcode in lower case
            if (numNoCalls + numMismatchesInBestBarcode < totalBarcodeReadBases && bestBarcodeMetric != null) {
                match.mismatches = numMismatchesInBestBarcode;
                match.mismatchesToSecondBest = numMismatchesInSecondBestBarcode;
                match.barcode = bestBarcodeMetric.BARCODE_WITHOUT_DELIMITER.toLowerCase();
            } else {
                match.mismatches = totalBarcodeReadBases;
                match.barcode = "";
            }

            if (match.matched) {
                match.barcode = bestBarcodeMetric.BARCODE_WITHOUT_DELIMITER;
            }
            return match;
        }

        private static void updateMetrics(final BarcodeMatch match, final boolean passingFilter,
                                          final Map<String, BarcodeMetric> metrics, final BarcodeMetric noMatchBarcodeMetric) {
            if (match.matched) {
                final BarcodeMetric matchMetric = metrics.get(match.barcode);
                ++matchMetric.READS;
                if (passingFilter) {
                    ++matchMetric.PF_READS;
                }
                if (match.mismatches == 0) {
                    ++matchMetric.PERFECT_MATCHES;
                    if (passingFilter) {
                        ++matchMetric.PF_PERFECT_MATCHES;
                    }
                } else if (match.mismatches == 1) {
                    ++matchMetric.ONE_MISMATCH_MATCHES;
                    if (passingFilter) {
                        ++matchMetric.PF_ONE_MISMATCH_MATCHES;
                    }
                }
            } else {
                ++noMatchBarcodeMetric.READS;
                if (passingFilter) {
                    ++noMatchBarcodeMetric.PF_READS;
                }
            }
        }
    }
}
