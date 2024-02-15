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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.programgroups.BaseCallingProgramGroup;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaDataType;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.util.IlluminaUtil;
import picard.util.ThreadPoolExecutorUtil;
import picard.util.ThreadPoolExecutorWithExceptions;

import java.io.BufferedWriter;
import java.io.File;
import java.text.NumberFormat;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
public class ExtractIlluminaBarcodes extends ExtractBarcodesProgram {
    static final String USAGE_SUMMARY = "Tool determines the barcode for each read in an Illumina lane.  ";
    static final String USAGE_DETAILS = "<p>This tool determines the numbers of reads containing barcode-matching sequences and provides " +
            "statistics on the quality of these barcode matches.</p> " +
            "<p>Illumina sequences can contain at least two types of barcodes, sample and molecular (index).  Sample barcodes " +
            "(B in the read structure) are used to demultiplex pooled samples while index barcodes (M in the read structure) are used " +
            "to differentiate multiple reads of a template when carrying out paired-end sequencing.  Note that this tool only extracts " +
            "sample (B) and not molecular barcodes (M).</p>" +
            "" +
            "<p>Barcodes can be provided in the form of a list (BARCODE_FILE) or a string representing the barcode (BARCODE).  " +
            "The BARCODE_FILE contains multiple fields including '" + BARCODE_SEQUENCE_COLUMN + "' (or '" + BARCODE_SEQUENCE_COLUMN + "_1'), '"
            + BARCODE_SEQUENCE_COLUMN + "_2' (optional), '" + BARCODE_NAME_COLUMN + "', and '" + LIBRARY_NAME_COLUMN + "'. " +
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

    @Argument(doc = "Tab-delimited file of barcode sequences, barcode name and, optionally, library name.  " +
            "Barcodes must be unique and all the same length.  Column headers must be '" +  BARCODE_SEQUENCE_COLUMN + "' (or '" + BARCODE_SEQUENCE_COLUMN + "_1'), '"
            + BARCODE_SEQUENCE_COLUMN + "_2' (optional), '"+ BARCODE_NAME_COLUMN + "', and '" + LIBRARY_NAME_COLUMN + "'.", mutex = {"BARCODE"})
    public File BARCODE_FILE;

    @Argument(doc = "Barcode sequence.  These must be unique, and all the same length.  This cannot be used with reads that " +
            "have more than one barcode; use BARCODE_FILE in that case. ", mutex = {"BARCODE_FILE"})
    public List<String> BARCODE = new ArrayList<>();

    @Argument(doc = "Run this many PerTileBarcodeExtractors in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.")
    public int NUM_PROCESSORS = 1;

    @Argument(doc = "Where to write _barcode.txt files.  By default, these are written to BASECALLS_DIR.", optional = true)
    public File OUTPUT_DIR;

    private static final Log LOG = Log.getInstance(ExtractIlluminaBarcodes.class);

    private final NumberFormat tileNumberFormatter = NumberFormat.getNumberInstance();

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

        final BarcodeExtractor barcodeExtractor = createBarcodeExtractor();
        final Set<IlluminaDataType> datatypes = (MINIMUM_BASE_QUALITY > 0) ?
                new HashSet<>(Arrays.asList(IlluminaDataType.BaseCalls, IlluminaDataType.PF, IlluminaDataType.QualityScores)) :
                new HashSet<>(Arrays.asList(IlluminaDataType.BaseCalls, IlluminaDataType.PF));
        final List<PerTileBarcodeExtractor> extractors = new ArrayList<>();
        for(Integer lane: LANE) {
            IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(BASECALLS_DIR, lane, inputReadStructure, bclQualityEvaluationStrategy, datatypes);

            for (final int tile : factory.getAvailableTiles()) {
                final PerTileBarcodeExtractor extractor = new PerTileBarcodeExtractor(
                        tile,
                        getBarcodeFile(lane, tile),
                        factory,
                        barcodeExtractor);

                extractors.add(extractor);
            }

            for (final PerTileBarcodeExtractor extractor : extractors) {
                pool.submit(extractor);
            }
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

        outputMetrics();
        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> messages = new ArrayList<>();

        INPUT_PARAMS_FILE = BARCODE_FILE;
        /*
          In extract illumina barcodes we NEVER want to look at the template reads nor the molecular barcodes, therefore replace them with
          skips because IlluminaDataProvider and its factory will neither open these nor produce ClusterData with the template reads in them,
          thus reducing the file IO and value copying done by the data provider
        */
        this.inputReadStructure = new ReadStructure(READ_STRUCTURE.replaceAll("[TM]", "S"));

        if (INPUT_PARAMS_FILE == null) {
            final int numBarcodes = inputReadStructure.sampleBarcodes.length();
            final Set<String> barcodes = new HashSet<>();

            for (final String barcode : BARCODE) {
                if (barcodes.contains(barcode)) {
                    messages.add("Barcode " + barcode + " specified more than once.");
                }
                barcodes.add(barcode);
                int barcodeNum = 0;
                int pos = 0;
                final String[] bcStrings = new String[numBarcodes];
                for (final ReadDescriptor rd : inputReadStructure.descriptors) {
                    if (rd.type != ReadType.Barcode) {
                        continue;
                    }
                    bcStrings[barcodeNum] = barcode.substring(pos, pos + rd.length);
                    pos += rd.length;
                    ++barcodeNum;
                }

                final BarcodeMetric metric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(bcStrings), bcStrings);
                barcodeToMetrics.put(barcode, metric);
            }
        }

        String[] superErrors = super.customCommandLineValidation();

        if ((INPUT_PARAMS_FILE != null || !BARCODE.isEmpty()) && barcodeToMetrics.keySet().isEmpty()) {
            messages.add("No barcodes have been specified.");
        }

        return collectErrorMessages(messages, superErrors);
    }



    /**
     * Create a barcode filename corresponding to the given tile qseq file.
     */
    private File getBarcodeFile(final int lane, final int tile) {
        return new File(OUTPUT_DIR,
                "s_" + lane + "_" + tileNumberFormatter.format(tile) + "_barcode.txt" + (COMPRESS_OUTPUTS ? ".gz" : ""));
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
        private BaseIlluminaDataProvider provider;
        private final IlluminaDataProviderFactory factory;
        private final ReadStructure outputReadStructure;
        private final BarcodeExtractor barcodeExtractor;

        /**
         * Constructor
         *  @param tile               The number of the tile being processed; used for logging only.
         *  @param barcodeFile        The file to write the barcodes to
         */
        public PerTileBarcodeExtractor(
                final int tile,
                final File barcodeFile,
                final IlluminaDataProviderFactory factory,
                final BarcodeExtractor extractor
        ) {
            this.barcodeExtractor = extractor;
            this.tile = tile;
            this.barcodeFile = barcodeFile;
            this.usingQualityScores = barcodeExtractor.getMinimumBaseQuality() > 0;
            this.metrics = new LinkedHashMap<>(barcodeExtractor.getMetrics().size());
            for (final String key : barcodeExtractor.getMetrics().keySet()) {
                this.metrics.put(key, barcodeExtractor.getMetrics().get(key).copy());
            }

            this.noMatch = barcodeExtractor.getNoMatchMetric().copy();
            this.factory = factory;
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
        public synchronized void run() {
            try {
                LOG.info("Extracting barcodes for tile " + tile);

                // Sometimes makeDataProvider takes a while waiting for slow file IO, for each tile the needed set of files
                // is non-overlapping sets of files so make the  data providers in the individual threads for PerTileBarcodeExtractors
                // so they are not all waiting for each others file operations. This also avoids opening numerous files
                // to create the provider prior to its use.
                this.provider = factory.makeDataProvider(tile);

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
                    final BarcodeExtractor.BarcodeMatch match = barcodeExtractor.findBestBarcode(barcodeSubsequences, qualityScores, false);

                    BarcodeExtractor.updateMetrics(match, passingFilter, metrics, noMatch);

                    final String yOrN = (match.isMatched() ? "Y" : "N");

                    for (final byte[] bc : barcodeSubsequences) {
                        writer.write(StringUtil.bytesToString(bc));
                    }
                    writer.write("\t" + yOrN + "\t" + match.getBarcode() + "\t" + match.getMismatches() + "\t" + match.getMismatchesToSecondBest());
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
    }
}
