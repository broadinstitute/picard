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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Illumina;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.IlluminaUtil;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.BufferedWriter;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
 * but we're close to the threshold of calling it a match we output the barcode that would have been
 * matched but in lower case
 *
 * @author jburke@broadinstitute.org
 */
@CommandLineProgramProperties(

        usage = ExtractIlluminaBarcodes.USAGE_SUMMARY + ExtractIlluminaBarcodes.USAGE_DETAILS,
        usageShort = ExtractIlluminaBarcodes.USAGE_SUMMARY,
        programGroup = Illumina.class
)
public class ExtractIlluminaBarcodes extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Tool determines the barcode for each read in an Illumina lane.  ";
    static final String USAGE_DETAILS = "<p>This tool determines the numbers of reads containing barcode-matching sequences and provides " +
            "statistics on the quality of these barcode matches.</p> " +
            "<p>Illumina sequences can contain at least two types of barcodes, sample and molecular (index).  Sample barcodes " +
            "(B in the read structure) are used to demultiplex pooled samples while index barcodes (M in the read structure) are used " +
            "to differentiate multiple reads of a template when carrying out paired-end sequencing.  Note that this tool only extracts " +
            "sample (B) and not molecular barcodes (M).</p>" +
            "" +
            "<p>Barcodes can be provided in the form of a list (BARCODE_FILE) or a string representing the barcode (BARCODE).  " +
            "The BARCODE_FILE contains multiple fields including 'barcode_sequence_1', 'barcode_sequence_2' (optional), " +
            "'barcode_name', and 'library_name'.  In contrast, the BARCODE argument is used for runs with reads containing a single " +
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

    @Option(doc = "The Illumina basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Option(doc = "Where to write _barcode.txt files.  By default, these are written to BASECALLS_DIR.", optional = true)
    public File OUTPUT_DIR;

    @Option(doc = "Lane number. ", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Option(doc = ReadStructure.PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Option(doc = "Barcode sequence.  These must be unique, and all the same length.  This cannot be used with reads that " +
            "have more than one barcode; use BARCODE_FILE in that case. ", mutex = {"BARCODE_FILE"})
    public List<String> BARCODE = new ArrayList<String>();

    @Option(doc = "Tab-delimited file of barcode sequences, barcode name and, optionally, library name.  " +
            "Barcodes must be unique and all the same length.  Column headers must be 'barcode_sequence_1', " +
            "'barcode_sequence_2' (optional), 'barcode_name', and 'library_name'.", mutex = {"BARCODE"})
    public File BARCODE_FILE;

    @Option(doc = "Per-barcode and per-lane metrics written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME)
    public File METRICS_FILE;

    @Option(doc = "Maximum mismatches for a barcode to be considered a match.")
    public int MAX_MISMATCHES = 1;

    @Option(doc = "Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.")
    public int MIN_MISMATCH_DELTA = 1;

    @Option(doc = "Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.")
    public int MAX_NO_CALLS = 2;

    @Option(shortName = "Q", doc = "Minimum base quality. Any barcode bases falling below this quality will be considered a mismatch even in the bases match.")
    public int MINIMUM_BASE_QUALITY = 0;

    @Option(doc = "The minimum quality (after transforming 0s to 1s) expected from reads.  If qualities are lower than this value, an error is thrown." +
            "The default of 2 is what the Illumina's spec describes as the minimum, but in practice the value has been observed lower.")
    public int MINIMUM_QUALITY = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

    @Option(shortName = "GZIP", doc = "Compress output s_l_t_barcode.txt files using gzip and append a .gz extension to the file names.")
    public boolean COMPRESS_OUTPUTS = false;

    @Option(doc = "Run this many PerTileBarcodeExtractors in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.")
    public int NUM_PROCESSORS = 1;

    private static final Log LOG = Log.getInstance(ExtractIlluminaBarcodes.class);

    /**
     * The read structure of the actual Illumina Run, i.e. the readStructure of the input data
     */
    private ReadStructure readStructure;

    private IlluminaDataProviderFactory factory;

    private final Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<String, BarcodeMetric>();

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
        final ExecutorService pool = Executors.newFixedThreadPool(numProcessors);

        // TODO: This is terribly inefficient; we're opening a huge number of files via the extractor constructor and we never close them.
        final List<PerTileBarcodeExtractor> extractors = new ArrayList<PerTileBarcodeExtractor>(factory.getAvailableTiles().size());
        for (final int tile : factory.getAvailableTiles()) {
            final PerTileBarcodeExtractor extractor = new PerTileBarcodeExtractor(
                    tile,
                    getBarcodeFile(tile),
                    barcodeToMetrics,
                    noMatchMetric,
                    factory,
                    MINIMUM_BASE_QUALITY,
                    MAX_NO_CALLS,
                    MAX_MISMATCHES,
                    MIN_MISMATCH_DELTA
            );
            extractors.add(extractor);
        }
        try {
            for (final PerTileBarcodeExtractor extractor : extractors) {
                pool.submit(extractor);
            }
            pool.shutdown();
            // Wait a while for existing tasks to terminate
            if (!pool.awaitTermination(6, TimeUnit.HOURS)) {
                pool.shutdownNow(); // Cancel any still-executing tasks
                // Wait a while for tasks to respond to being cancelled
                if (!pool.awaitTermination(60, TimeUnit.SECONDS))
                    LOG.error("Pool did not terminate");
                return 1;
            }
        } catch (final Throwable e) {
            // (Re-)Cancel if current thread also interrupted
            LOG.error(e, "Parent thread encountered problem submitting extractors to thread pool or awaiting shutdown of threadpool.  Attempting to kill threadpool.");
            pool.shutdownNow();
            return 2;
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
        for (Map.Entry<Byte, Integer> entry : bclQualityEvaluationStrategy.getPoorQualityFrequencies().entrySet()) {
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
        int totalReads = noMatchMetric.READS;
        int totalPfReads = noMatchMetric.PF_READS;
        int totalPfReadsAssigned = 0;
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
        final ArrayList<String> messages = new ArrayList<String>();

        this.bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);

        /**
         * In extract illumina barcodes we NEVER want to look at the template reads nor the molecular barcodes, therefore replace them with
         * skips because IlluminaDataProvider and its factory will neither open these nor produce ClusterData with the template reads in them,
         * thus reducing the file IO and value copying done by the data provider
         */
        readStructure = new ReadStructure(READ_STRUCTURE.replaceAll("T|M", "S"));
        final IlluminaDataType[] datatypes = (MINIMUM_BASE_QUALITY > 0) ?
                new IlluminaDataType[]{IlluminaDataType.BaseCalls, IlluminaDataType.PF, IlluminaDataType.QualityScores} :
                new IlluminaDataType[]{IlluminaDataType.BaseCalls, IlluminaDataType.PF};
        factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, bclQualityEvaluationStrategy, datatypes);

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
        if (barcodeToMetrics.keySet().isEmpty()) {
            messages.add("No barcodes have been specified.");
        }
        if (messages.isEmpty()) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }

    public static void main(final String[] argv) {
        new ExtractIlluminaBarcodes().instanceMainWithExit(argv);
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
        final int numBarcodes = readStructure.sampleBarcodes.length();
        final Set<String> barcodes = new HashSet<String>();
        for (final TabbedTextFileWithHeaderParser.Row row : barcodesParser) {
            final String[] bcStrings = new String[numBarcodes];
            int barcodeNum = 1;
            for (final ReadDescriptor rd : readStructure.descriptors) {
                if (rd.type != ReadType.Barcode) continue;
                final String header = barcodeNum == 1 ? sequenceColumn : "barcode_sequence_" + String.valueOf(barcodeNum);
                bcStrings[barcodeNum - 1] = row.getField(header);
                barcodeNum++;
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
        barcodesParser.close();
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
        private final BaseIlluminaDataProvider provider;
        private final ReadStructure outputReadStructure;
        private final int maxNoCalls, maxMismatches, minMismatchDelta, minimumBaseQuality;

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
                final BarcodeMetric noMatchMetric,
                final IlluminaDataProviderFactory factory,
                final int minimumBaseQuality,
                final int maxNoCalls,
                final int maxMismatches,
                final int minMismatchDelta
        ) {
            this.tile = tile;
            this.barcodeFile = barcodeFile;
            this.usingQualityScores = minimumBaseQuality > 0;
            this.maxNoCalls = maxNoCalls;
            this.maxMismatches = maxMismatches;
            this.minMismatchDelta = minMismatchDelta;
            this.minimumBaseQuality = minimumBaseQuality;
            this.metrics = new LinkedHashMap<String, BarcodeMetric>(barcodeToMetrics.size());
            for (final String key : barcodeToMetrics.keySet()) {
                this.metrics.put(key, BarcodeMetric.copy(barcodeToMetrics.get(key)));
            }
            this.noMatch = BarcodeMetric.copy(noMatchMetric);
            this.provider = factory.makeDataProvider(Arrays.asList(tile));
            this.outputReadStructure = factory.getOutputReadStructure();

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
        synchronized public void run() {
            try {
                LOG.info("Extracting barcodes for tile " + tile);

                //Sometimes makeDataProvider takes a while waiting for slow file IO, for each tile the needed set of files
                //is non-overlapping sets of files so make the  data providers in the individual threads for PerTileBarcodeExtractors
                //so they are not all waiting for each others file operations

                //Most likely we have SKIPS in our read structure since we replace all template reads with skips in the input data structure
                //(see customCommnandLineValidation), therefore we must use the outputReadStructure to index into the output cluster data
                final int[] barcodeIndices = outputReadStructure.sampleBarcodes.getIndices();
                final BufferedWriter writer = IOUtil.openFileForBufferedWriting(barcodeFile);
                final byte[][] barcodeSubsequences = new byte[barcodeIndices.length][];
                final byte[][] qualityScores = usingQualityScores ? new byte[barcodeIndices.length][] : null;
                while (provider.hasNext()) {
                    // Extract the barcode from the cluster and write it to the file for the tile
                    final ClusterData cluster = provider.next();
                    for (int i = 0; i < barcodeIndices.length; i++) {
                        barcodeSubsequences[i] = cluster.getRead(barcodeIndices[i]).getBases();
                        if (usingQualityScores) qualityScores[i] = cluster.getRead(barcodeIndices[i]).getQualities();
                    }
                    final boolean passingFilter = cluster.isPf();
                    final BarcodeMatch match = BarcodeMatch.findBestBarcodeAndUpdateMetrics(barcodeSubsequences,
                            qualityScores, passingFilter, metrics, noMatch, maxNoCalls, maxMismatches,
                            minMismatchDelta, minimumBaseQuality);

                    final String yOrN = (match.isMatched() ? "Y" : "N");

                    for (final byte[] bc : barcodeSubsequences) {
                        writer.write(StringUtil.bytesToString(bc));
                    }
                    writer.write("\t" + yOrN + "\t" + match.getBarcode() + "\t" + String.valueOf(match.getMismatches()) +
                            "\t" + String.valueOf(match.getMismatchesToSecondBest()));
                    writer.newLine();
                }
                writer.close();
            } catch (final Exception e) {
                LOG.error(e, "Error processing tile ", this.tile);
                this.exception = e;
            } finally {
                provider.close();
            }
        }
    }
}
