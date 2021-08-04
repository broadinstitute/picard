package picard.illumina;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.util.Tuple;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.IlluminaUtil;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public abstract class ExtractBarcodesProgram extends CommandLineProgram {
    @Argument(doc = "The distance metric that should be used to compare the barcode-reads and the provided barcodes for finding the best and second-best assignments.")
    public DistanceMetric DISTANCE_MODE = DistanceMetric.HAMMING;

    @Argument(doc = "Maximum mismatches for a barcode to be considered a match.")
    public int MAX_MISMATCHES = 1;

    @Argument(doc = "Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.")
    public int MIN_MISMATCH_DELTA = 1;

    @Argument(doc = "Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.")
    public int MAX_NO_CALLS = 2;

    @Argument(shortName = "Q", doc = "Minimum base quality. Any barcode bases falling below this quality will be considered a mismatch even if the bases match.")
    public int MINIMUM_BASE_QUALITY = 0;

    @Argument(doc = "The minimum quality (after transforming 0s to 1s) expected from reads.  If qualities are lower than this value, an error is thrown. " +
            "The default of 2 is what the Illumina's spec describes as the minimum, but in practice the value has been observed lower.")
    public int MINIMUM_QUALITY = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

    @Argument(doc = "Lane number. This can be specified multiple times. Reads with the same index in multiple lanes" +
            " will be added to the same output file.", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public List<Integer> LANE;

    @Argument(doc = ReadStructure.PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Argument(shortName = "GZIP", doc = "Compress output FASTQ files using gzip and append a .gz extension to the file names.")
    public boolean COMPRESS_OUTPUTS = false;

    @Argument(doc = "The Illumina basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(doc = "Per-barcode and per-lane metrics written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, optional = true)
    public File METRICS_FILE;

    @Argument(doc = "The input file that defines parameters for the program. This is the BARCODE_FILE for" +
            " `ExtractIlluminaBarcodes` or the MULTIPLEX_PARAMS or LIBRARY_PARAMS file for `IlluminaBasecallsToFastq` " +
            " or `IlluminaBasecallsToSam`", optional = true)
    public File INPUT_PARAMS_FILE;
    /**
     * Column header for the first barcode sequence (preferred).
     */
    public static final String BARCODE_COLUMN = "barcode";
    public static final String BARCODE_SEQUENCE_COLUMN = "barcode_sequence";

    /**
     * Column header for the barcode name.
     */
    public static final String BARCODE_NAME_COLUMN = "barcode_name";
    /**
     * Column header for the library name.
     */
    public static final String LIBRARY_NAME_COLUMN = "library_name";
    public static final Set<String> BARCODE_PREFIXES = new HashSet<>(
            Arrays.asList(BARCODE_SEQUENCE_COLUMN, BARCODE_COLUMN)
    );



    protected Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<>();
    protected final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);
    protected BarcodeMetric noMatchMetric;
    private final NumberFormat tileNumberFormatter = NumberFormat.getNumberInstance();

    /**
     * The read structure of the actual Illumina Run, i.e. the readStructure of the input data
     */
    protected ReadStructure inputReadStructure;

    protected BarcodeExtractor createBarcodeExtractor() {
        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[inputReadStructure.sampleBarcodes.length()];
        int index = 0;
        for (final ReadDescriptor d : inputReadStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                noMatchBarcode[index++] = StringUtil.repeatCharNTimes('N', d.length);
            }
        }
        this.noMatchMetric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

        return new BarcodeExtractor(barcodeToMetrics,
                noMatchMetric,
                inputReadStructure,
                MAX_NO_CALLS,
                MAX_MISMATCHES,
                MIN_MISMATCH_DELTA,
                MINIMUM_BASE_QUALITY,
                DISTANCE_MODE);
    }

    /**
     * Parses all barcodes from input files and validates all barcodes are the same length and unique
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     * to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        inputReadStructure = new ReadStructure(READ_STRUCTURE);
        List<String> messages = new ArrayList<>();
        tileNumberFormatter.setMinimumIntegerDigits(4);
        tileNumberFormatter.setGroupingUsed(false);

        if (INPUT_PARAMS_FILE != null) {
            Tuple<Map<String, BarcodeMetric>, List<String>> test = parseInputFile(INPUT_PARAMS_FILE, inputReadStructure);
            barcodeToMetrics = test.a;
            messages = test.b;

            if (barcodeToMetrics.keySet().isEmpty()) {
                messages.add("No barcodes have been specified.");
            }
        }

       return messages.toArray(new String[0]);
    }

    protected String[] collectErrorMessages(List<String> messages, String[] superErrors) {
        if (superErrors != null && superErrors.length > 0) {
            messages.addAll(Arrays.asList(superErrors));
        }

        if (messages.isEmpty()) {
            return null;
        }
        return messages.toArray(new String[0]);
    }

    protected void outputMetrics() {
        final MetricsFile<BarcodeMetric, Integer> metrics = getMetricsFile();
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            metrics.addMetric(barcodeMetric);
        }
        metrics.addMetric(noMatchMetric);
        metrics.write(METRICS_FILE);
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
     * Parses any one of the following types of files:
     *
     * ExtractIlluminaBarcodes      BARCODE_FILE
     * IlluminaBasecallsToFastq     MULTIPLEX_PARAMS
     * IlluminaBasecallsToSam       LIBRARY_PARAMS
     *
     * This will validate to file format as well as populate a Map of barcodes to metrics.
     *
     * @param inputFile         The input file that is being parsed
     * @param readStructure     The read structure for the reads of the run
     */
    protected static Tuple<Map<String, BarcodeMetric>, List<String>> parseInputFile(final File inputFile,
                                                                                    final ReadStructure readStructure) {
        List<String> messages = new ArrayList<>();
        Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<>();
        try (final TabbedTextFileWithHeaderParser barcodesParser = new TabbedTextFileWithHeaderParser(inputFile)) {
            List<String> validBarcodeColumns = barcodesParser.columnLabels().stream().filter(name -> {
                boolean isValidPrefix = false;
                for (String columnPrefix : BARCODE_PREFIXES) {
                    isValidPrefix |= name.toUpperCase().startsWith(columnPrefix.toUpperCase()) && !name.equalsIgnoreCase(BARCODE_NAME_COLUMN);
                }
                return isValidPrefix;
            }).collect(Collectors.toList());

            if (readStructure.sampleBarcodes.length() != validBarcodeColumns.size()) {
                messages.add("Expected " + readStructure.sampleBarcodes.length() + " valid barcode columns, but found " +
                        String.join(",", validBarcodeColumns));
            }

            validBarcodeColumns.sort((s, t1) -> {
                int lengthDiff =  s.length() - t1.length();
                if (lengthDiff == 0) {
                    return s.compareTo(t1);
                }
                return lengthDiff;
            });

            Matcher matcher = Pattern.compile("^(.*)_\\d").matcher(validBarcodeColumns.get(0));

            final String sequenceColumn;
            boolean hasMultipleNumberedBarcodeColumns = matcher.matches();

            // If there are multiple numbered barcode columns we extract the base column name (ie for BARCODE_1, BARCODE_2)
            // the sequence column is BARCODE. In the case of BARCODE_PARAMS file the columns are BARCODE_SEQUENCE_1
            // and BARCODE_SEQUENCE_2.
            if (hasMultipleNumberedBarcodeColumns) {
                sequenceColumn = matcher.group(1);
            } else {
                sequenceColumn = validBarcodeColumns.get(0);
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
                    final String header = hasMultipleNumberedBarcodeColumns ? sequenceColumn + "_" + (1 + barcodeNum) : sequenceColumn;
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
                // if the barcode is all Ns don't add it to metrics (we add noCallMetric separately)
                if (bcStr.contains("N") || bcStr.contains("n")) {
                    continue;
                }
                if (barcodes.contains(bcStr)) {
                    messages.add("Barcode " + bcStr + " specified more than once in " + inputFile);
                }
                barcodes.add(bcStr);
                final String barcodeName = (hasBarcodeName ? row.getField(BARCODE_NAME_COLUMN) : "");
                final String libraryName = (hasLibraryName ? row.getField(LIBRARY_NAME_COLUMN) : "");
                final BarcodeMetric metric = new BarcodeMetric(barcodeName, libraryName, bcStr, bcStrings);
                barcodeToMetrics.put(StringUtil.join("", bcStrings), metric);
            }
        }
        return new Tuple<>(barcodeToMetrics, messages);
    }
}
