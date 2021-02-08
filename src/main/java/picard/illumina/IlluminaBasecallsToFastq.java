/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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

import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.BaseCallingProgramGroup;
import picard.fastq.Casava18ReadNameEncoder;
import picard.fastq.IlluminaReadNameEncoder;
import picard.fastq.ReadNameEncoder;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.ReadData;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.IlluminaUtil;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

@CommandLineProgramProperties(
        summary = IlluminaBasecallsToFastq.USAGE_SUMMARY + IlluminaBasecallsToFastq.USAGE_DETAILS,
        oneLineSummary = IlluminaBasecallsToFastq.USAGE_SUMMARY,
        programGroup = BaseCallingProgramGroup.class
)
@DocumentedFeature
public class IlluminaBasecallsToFastq extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Generate FASTQ file(s) from Illumina basecall read data.  ";
    static final String USAGE_DETAILS = "<p>This tool generates FASTQ files from data in an Illumina BaseCalls output directory.  " +
            "Separate FASTQ files are created for each template, barcode, and index (molecular barcode) read.  Briefly, the template reads " +
            "are the target sequence of your experiment, the barcode sequence reads facilitate sample demultiplexing, and the index reads " +
            "help mitigate instrument phasing errors.  For additional information on the read types, please see the following " +
            "reference <a href'=http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245947/'>here</a>.</p>" +
            "" +
            "<p>In the absence of sample pooling (multiplexing) and/or barcodes, then an OUTPUT_PREFIX (file directory) must be " +
            "provided as the sample identifier.  For multiplexed samples, a MULTIPLEX_PARAMS file must be specified.  " +
            "The MULTIPLEX_PARAMS file contains the list of sample barcodes used to sort template, barcode, and index reads.  " +
            "It is essentially the same as the BARCODE_FILE used in the" +
            "<a href='http://broadinstitute.github.io/picard/command-line-overview.html#ExtractIlluminaBarcodes'>ExtractIlluminaBarcodes</a> " +
            "tool.</p>     " +
            "" +
            "<p>Files from this tool use the following naming format: {prefix}.{type}_{number}.fastq with the {prefix} indicating the sample " +
            "barcode, the {type} indicating the types of reads e.g. index, barcode, or blank (if it contains a template read).  " +
            "The {number} indicates the read number, either first (1) or second (2) for paired-end sequencing. </p> " +

            "<h4>Usage examples:</h4>" +
            "<pre>" +
            "Example 1: Sample(s) with either no barcode or barcoded without multiplexing <br />" +
            "java -jar picard.jar IlluminaBasecallsToFastq \\<br />" +
            "      READ_STRUCTURE=25T8B25T \\<br />" +
            "      BASECALLS_DIR=basecallDirectory \\<br />" +
            "      LANE=001 \\<br />" +
            "      OUTPUT_PREFIX=noBarcode.1 \\<br />" +
            "      RUN_BARCODE=run15 \\<br />" +
            "      FLOWCELL_BARCODE=abcdeACXX <br /><br />" +

            "Example 2: Multiplexed samples <br />" +
            "java -jar picard.jar IlluminaBasecallsToFastq \\<br />" +
            "      READ_STRUCTURE=25T8B25T \\<br />" +
            "      BASECALLS_DIR=basecallDirectory \\<br />" +
            "      LANE=001 \\<br />" +
            "      MULTIPLEX_PARAMS=demultiplexed_output.txt \\<br />" +
            "      RUN_BARCODE=run15 \\<br />" +
            "      FLOWCELL_BARCODE=abcdeACXX <br />" +
            "</pre>" +
            "<p>The FLOWCELL_BARCODE is required if emitting Casava 1.8-style read name headers.</p>" +
            "<hr />";

    // The following attributes define the command-line arguments

    @Argument(doc = "The basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(doc = "The barcodes directory with _barcode.txt files (generated by ExtractIlluminaBarcodes). If not set, use BASECALLS_DIR. ", shortName = "BCD", optional = true)
    public File BARCODES_DIR;

    @Argument(doc = "Lane number. ", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Argument(doc = "The prefix for output FASTQs.  Extensions as described above are appended.  Use this option for a non-barcoded run, or" +
            " for a barcoded run in which it is not desired to demultiplex reads into separate files by barcode.",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            mutex = {"MULTIPLEX_PARAMS"})
    public File OUTPUT_PREFIX;

    @Argument(doc = "The barcode of the run.  Prefixed to read names.")
    public String RUN_BARCODE;

    @Argument(doc = "The name of the machine on which the run was sequenced; required if emitting Casava1.8-style read name headers", optional = true)
    public String MACHINE_NAME;

    @Argument(doc = "The barcode of the flowcell that was sequenced; required if emitting Casava1.8-style read name headers", optional = true)
    public String FLOWCELL_BARCODE;

    @Argument(doc = ReadStructure.PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Argument(doc = "Tab-separated file for creating all output FASTQs demultiplexed by barcode for a lane with single " +
            "IlluminaBasecallsToFastq invocation.  The columns are OUTPUT_PREFIX, and BARCODE_1, BARCODE_2 ... BARCODE_X " +
            "where X = number of barcodes per cluster (optional).  Row with BARCODE_1 set to 'N' is used to specify " +
            "an output_prefix for no barcode match.",
            mutex = {"OUTPUT_PREFIX"})
    public File MULTIPLEX_PARAMS;

    @Deprecated
    @Argument(doc = "Deprecated (No longer used). Which adapters to look for in the read.", optional = true)
    public List<IlluminaUtil.IlluminaAdapterPair> ADAPTERS_TO_CHECK = null;

    @Argument(doc = "The number of threads to run in parallel. If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0, then the number of cores used will" +
            " be the number available on the machine less NUM_PROCESSORS.")
    public Integer NUM_PROCESSORS = 0;

    @Argument(doc = "If set, this is the first tile to be processed (used for debugging).  Note that tiles are not processed" +
            " in numerical order.",
            optional = true)
    public Integer FIRST_TILE;

    @Argument(doc = "If set, process no more than this many tiles (used for debugging).", optional = true)
    public Integer TILE_LIMIT;

    @Argument(doc = "Apply EAMSS filtering to identify inappropriately quality scored bases towards the ends of reads" +
            " and convert their quality scores to Q2.")
    public boolean APPLY_EAMSS_FILTER = true;

    @Argument(doc = "If true, call System.gc() periodically.  This is useful in cases in which the -Xmx value passed " +
            "is larger than the available memory.")
    public Boolean FORCE_GC = true;

    @Argument(doc = "If true, the output records are sorted by read name. Otherwise they are output in the same order " +
            "that the data was produced on the sequencer (ordered by tile and position).")
    public Boolean SORT = true;

    @Deprecated
    @Argument(doc = "Configure SortingCollections to store this many records before spilling to disk. For an indexed" +
            " run, each SortingCollection gets this value/number of indices. Deprecated: use `MAX_RECORDS_IN_RAM`")
    public int MAX_READS_IN_RAM_PER_TILE = -1;

    @Argument(doc = "The minimum quality (after transforming 0s to 1s) expected from reads.  If qualities are lower than this value, an error is thrown." +
            "The default of 2 is what the Illumina's spec describes as the minimum, but in practice the value has been observed lower.")
    public int MINIMUM_QUALITY = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

    @Argument(doc = "Whether to include non-PF reads", shortName = "NONPF", optional = true)
    public boolean INCLUDE_NON_PF_READS = true;

    @Argument(doc = "Whether to ignore reads whose barcodes are not found in MULTIPLEX_PARAMS.  Useful when outputting " +
            "FASTQs for only a subset of the barcodes in a lane.", shortName = "INGORE_UNEXPECTED")
    public boolean IGNORE_UNEXPECTED_BARCODES = false;

    @Argument(doc = "The read name header formatting to emit.  Casava1.8 formatting has additional information beyond Illumina, including: " +
            "the passing-filter flag value for the read, the flowcell name, and the sequencer name.")
    public ReadNameFormat READ_NAME_FORMAT = ReadNameFormat.CASAVA_1_8;

    @Argument(shortName = "GZIP", doc = "Compress output FASTQ files using gzip and append a .gz extension to the file names.")
    public boolean COMPRESS_OUTPUTS = false;

    /**
     * Simple switch to control the read name format to emit.
     */
    public enum ReadNameFormat {
        CASAVA_1_8, ILLUMINA
    }
    private final Map<String, AsyncClusterWriter> sampleBarcodeClusterWriterMap = new HashMap<>(1, 0.5f);
    final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);
    private ReadStructure readStructure;
    private BasecallsConverter<?> basecallsConverter;
    private static final Log log = Log.getInstance(IlluminaBasecallsToFastq.class);
    private final FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
    private ReadNameEncoder readNameEncoder;

    @Override
    protected int doWork() {
        initialize();
        final Set<String> barcodes = sampleBarcodeClusterWriterMap.keySet();
        basecallsConverter.processTilesAndWritePerSampleOutputs(barcodes);
        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        final LinkedList<String> errors = new LinkedList<>();

        // Remove once deprecated parameter is deleted.
        if (MAX_READS_IN_RAM_PER_TILE != -1) {
            log.warn("Setting deprecated parameter `MAX_READS_IN_RAM_PER_TILE` use ` MAX_RECORDS_IN_RAM` instead");
            MAX_RECORDS_IN_RAM = MAX_READS_IN_RAM_PER_TILE;
        }

        if (READ_NAME_FORMAT == ReadNameFormat.CASAVA_1_8 && MACHINE_NAME == null) {
            errors.add("MACHINE_NAME is required when using Casava1.8-style read name headers.");
        }

        if (READ_NAME_FORMAT == ReadNameFormat.CASAVA_1_8 && FLOWCELL_BARCODE == null) {
            errors.add("FLOWCELL_BARCODE is required when using Casava1.8-style read name headers.");
        }

        if (ADAPTERS_TO_CHECK != null) {
            log.warn("ADAPTERS_TO_CHECK is not used");
        }

        if (errors.isEmpty()) {
            return null;
        } else {
            return errors.toArray(new String[errors.size()]);
        }
    }

    /**
     * Prepares loggers, initiates garbage collection thread, parses arguments and initialized variables appropriately/
     */
    private void initialize() {
        fastqWriterFactory.setCreateMd5(CREATE_MD5_FILE);

        switch (READ_NAME_FORMAT) {
            case CASAVA_1_8:
                readNameEncoder = new Casava18ReadNameEncoder(MACHINE_NAME, RUN_BARCODE, FLOWCELL_BARCODE);
                break;
            case ILLUMINA:
                readNameEncoder = new IlluminaReadNameEncoder(RUN_BARCODE);
                break;
        }

        readStructure = new ReadStructure(READ_STRUCTURE);
        if (MULTIPLEX_PARAMS != null) {
            IOUtil.assertFileIsReadable(MULTIPLEX_PARAMS);
        }
        final boolean demultiplex;
        if (OUTPUT_PREFIX != null) {
            sampleBarcodeClusterWriterMap.put(null, buildWriter(OUTPUT_PREFIX));
            demultiplex = false;
        } else {
            populateWritersFromMultiplexParams();
            demultiplex = true;
        }

        BasecallsConverterBuilder<ClusterData> converterBuilder = new BasecallsConverterBuilder<>(BASECALLS_DIR, LANE, readStructure, sampleBarcodeClusterWriterMap)
                .barcodesDir(BARCODES_DIR)
                .withDemultiplex(demultiplex)
                .numProcessors(NUM_PROCESSORS)
                .firstTile(FIRST_TILE)
                .tileLimit(TILE_LIMIT)
                .withMaxRecordsInRam(MAX_RECORDS_IN_RAM)
                .withApplyEamssFiltering(APPLY_EAMSS_FILTER)
                .withIncludeNonPfReads(INCLUDE_NON_PF_READS)
                .withIgnoreUnexpectedBarcodes(IGNORE_UNEXPECTED_BARCODES)
                .withBclQualityEvaluationStrategy(bclQualityEvaluationStrategy);

        if (SORT) {
            Comparator<ClusterData> queryNameComparator = new ClusterDataQueryNameComparator(readNameEncoder);
            converterBuilder = converterBuilder.withSorting(
                    queryNameComparator,
                    new ClusterDataCodec(),
                    ClusterData.class,
                    TMP_DIR);
        }


        final BasecallsConverter<ClusterData> converter = converterBuilder.build();
        converter.setConverter(new NoOpClusterConverter());
        this.basecallsConverter = converter;

        log.info("READ STRUCTURE IS " + readStructure.toString());
    }


    /**
     * Assert that expectedCols are present
     *
     * @param actualCols   The columns present in the MULTIPLEX_PARAMS file
     * @param expectedCols The columns that are REQUIRED
     */
    private void assertExpectedColumns(final Set<String> actualCols, final Set<String> expectedCols) {
        final Set<String> missingColumns = new HashSet<>(expectedCols);
        missingColumns.removeAll(actualCols);

        if (!missingColumns.isEmpty()) {
            throw new PicardException(String.format(
                    "MULTIPLEX_PARAMS file %s is missing the following columns: %s.",
                    MULTIPLEX_PARAMS.getAbsolutePath(), StringUtil.join(", ", missingColumns
                    )));
        }
    }

    /**
     * For each line in the MULTIPLEX_PARAMS file create a FastqRecordsWriter and put it in the sampleBarcodeFastqWriterMap map,
     * where the key to the map is the concatenation of all sampleBarcodes in order for the given line.
     */
    private void populateWritersFromMultiplexParams() {
        final TabbedTextFileWithHeaderParser libraryParamsParser = new TabbedTextFileWithHeaderParser(MULTIPLEX_PARAMS);

        final Set<String> expectedColumnLabels = CollectionUtil.makeSet("OUTPUT_PREFIX");
        final List<String> sampleBarcodeColumnLabels = new ArrayList<>();
        for (int i = 1; i <= readStructure.sampleBarcodes.length(); i++) {
            sampleBarcodeColumnLabels.add("BARCODE_" + i);
        }

        expectedColumnLabels.addAll(sampleBarcodeColumnLabels);
        assertExpectedColumns(libraryParamsParser.columnLabels(), expectedColumnLabels);

        final List<TabbedTextFileWithHeaderParser.Row> rows = libraryParamsParser.iterator().toList();
        final Set<String> seenBarcodes = new HashSet<>();

        for (final TabbedTextFileWithHeaderParser.Row row : rows) {
            List<String> sampleBarcodeValues = null;

            if (!sampleBarcodeColumnLabels.isEmpty()) {
                sampleBarcodeValues = new ArrayList<>();
                for (final String sampleBarcodeLabel : sampleBarcodeColumnLabels) {
                    sampleBarcodeValues.add(row.getField(sampleBarcodeLabel));
                }
            }

            final String key = (sampleBarcodeValues == null || sampleBarcodeValues.contains("N")) ? null : StringUtil.join("", sampleBarcodeValues);
            if (seenBarcodes.contains(key)) {    //This will catch the case of having more than 1 line in a non-barcoded MULTIPLEX_PARAMS file
                throw new PicardException("Row for barcode " + key + " appears more than once in MULTIPLEX_PARAMS file " + MULTIPLEX_PARAMS);
            } else {
                seenBarcodes.add(key);
            }

            sampleBarcodeClusterWriterMap.put(key, buildWriter(new File(row.getField("OUTPUT_PREFIX"))));
        }

        if (seenBarcodes.isEmpty()) {
            throw new PicardException("MULTIPLEX_PARAMS file " + MULTIPLEX_PARAMS + " does have any data rows.");
        }
        libraryParamsParser.close();
    }

    /**
     * Builds an asynchronous writer used to output all fastq records.
     *
     * @return AsyncClusterWriter that contains one or more ClusterWriters (amount depends on read structure), all using
     * outputPrefix to determine the filename(s).
     */
    private AsyncClusterWriter buildWriter(final File outputPrefix) {
        final File outputDir = outputPrefix.getAbsoluteFile().getParentFile();
        IOUtil.assertDirectoryIsWritable(outputDir);
        final String prefixString = outputPrefix.getName();
        final String suffixString = COMPRESS_OUTPUTS ? "fastq.gz" : "fastq";

        final File[] templateFiles = new File[readStructure.templates.length()];
        final File[] sampleBarcodeFiles = new File[readStructure.sampleBarcodes.length()];
        final File[] molecularBarcodeFiles = new File[readStructure.molecularBarcode.length()];

        for (int i = 0; i < templateFiles.length; ++i) {
            templateFiles[i] = new File(outputDir, String.format("%s.%d.%s", prefixString, i + 1, suffixString));
        }

        for (int i = 0; i < sampleBarcodeFiles.length; ++i) {
            sampleBarcodeFiles[i] = new File(outputDir, String.format("%s.barcode_%d.%s", prefixString, i + 1, suffixString));
        }

        for (int i = 0; i < molecularBarcodeFiles.length; ++i) {
            molecularBarcodeFiles[i] = new File(outputDir, String.format("%s.index_%d.%s", prefixString, i + 1, suffixString));
        }


        return new AsyncClusterWriter(new ClusterToFastqWriter(templateFiles, sampleBarcodeFiles, molecularBarcodeFiles), 1024);
    }

    /**
     * Trivial class to avoid converting ClusterData to another type when not sorting outputs.
     */
    private static final class NoOpClusterConverter implements BasecallsConverter.ClusterDataConverter<ClusterData> {
        @Override
        public ClusterData convertClusterToOutputRecord(ClusterData cluster) {
            return cluster;
        }
    }

    private final class AsyncClusterWriter extends AbstractAsyncWriter<ClusterData> implements BasecallsConverter.ConvertedClusterDataWriter<ClusterData>  {
        private final ClusterToFastqWriter writer;

        public AsyncClusterWriter(final ClusterToFastqWriter out, final int queueSize) {
            super(queueSize);
            this.writer = out;
        }

        @Override protected String getThreadNamePrefix() { return "FastqWriterThread-"; }
        @Override protected void synchronouslyWrite(final ClusterData item) { this.writer.write(item); }
        @Override protected void synchronouslyClose() { this.writer.close(); }
    }
    /**
     * An optimized writer for writing ClusterData directly to a set of Fastq files.
     */
    private final class ClusterToFastqWriter implements BasecallsConverter.ConvertedClusterDataWriter<ClusterData> {
        public static final char NEW_LINE = '\n';
        public static final char AT_SYMBOL = '@';
        public static final char PLUS = '+';
        private final OutputStream[] templateOut;
        private final OutputStream[] sampleBarcodeOut;
        private final OutputStream[] molecularBarcodeOut;
        private final boolean appendTemplateNumber;
        private final boolean appendMolecularBarcodeNumber;
        private final int numReads;

        public ClusterToFastqWriter(final File[] templateFiles,
                                    final File[] sampleBarcodeFiles,
                                    final File[] molecularBarcodeFiles) {

            this.templateOut = Arrays.stream(templateFiles).map(this::makeWriter).toArray(OutputStream[]::new);
            this.sampleBarcodeOut = Arrays.stream(sampleBarcodeFiles).map(this::makeWriter).toArray(OutputStream[]::new);
            this.molecularBarcodeOut = Arrays.stream(molecularBarcodeFiles).map(this::makeWriter).toArray(OutputStream[]::new);
            this.appendTemplateNumber = this.templateOut.length > 1;
            this.appendMolecularBarcodeNumber = this.molecularBarcodeOut.length > 1;
            this.numReads = templateOut.length + sampleBarcodeOut.length + molecularBarcodeOut.length;
        }

        private OutputStream makeWriter(final File file) {
            Path outputPath = file.toPath();
            try {
                OutputStream os = Files.newOutputStream(outputPath);
                if (IOUtil.hasGzipFileExtension(outputPath)) {
                    os = new BlockCompressedOutputStream(os, (File) null, COMPRESSION_LEVEL);
                } else {
                    os = IOUtil.maybeBufferOutputStream(os);
                }
                if (Defaults.CREATE_MD5) os = new Md5CalculatingOutputStream(os, IOUtil.addExtension(outputPath, ".md5"));
                return os;
            } catch (final IOException ioe) {
                throw new RuntimeIOException("Error opening file: " + outputPath.toUri(), ioe);
            }
        }

        @Override
        public void write(final ClusterData rec) {
            int templateIndex = 0;
            int sampleBarcodeIndex = 0;
            int molecularBarcodeIndex = 0;

            for (int i = 0; i < this.numReads; ++i) {
                final ReadData read = rec.getRead(i);
                final OutputStream out;
                final String name;

                switch (read.getReadType()) {
                    case T:
                        out = templateOut[templateIndex++];
                        name = readNameEncoder.generateReadName(rec, appendTemplateNumber ? templateIndex : null);
                        break;
                    case B:
                        out = sampleBarcodeOut[sampleBarcodeIndex++];
                        name = readNameEncoder.generateReadName(rec, null);
                        break;
                    case M:
                        out = molecularBarcodeOut[molecularBarcodeIndex++];
                        name = readNameEncoder.generateReadName(rec, appendMolecularBarcodeNumber ? molecularBarcodeIndex : null);
                        break;
                    default:
                        throw new IllegalStateException("Read type other than T/B/M encountered.");
                }

                writeSingle(out, name, read);
            }
        }

        /**
         * Writes out a single read to a single FASTQ's output stream.
         */
        private void writeSingle(final OutputStream out, final String name, final ReadData read) {
            try {
                final byte[] bases = read.getBases();
                final byte[] quals = read.getQualities();
                final int len = bases.length;
                for (int i = 0; i < len; ++i) {
                    quals[i] = (byte) SAMUtils.phredToFastq(quals[i]);
                }

                out.write(AT_SYMBOL);
                out.write(name.getBytes(StandardCharsets.UTF_8));
                out.write(NEW_LINE);
                out.write(bases);
                out.write(NEW_LINE);
                out.write(PLUS);
                out.write(NEW_LINE);
                out.write(quals);
                out.write(NEW_LINE);
            } catch (IOException ioe) {
                throw new RuntimeIOException(ioe);
            }
        }

        @Override
        public void close() {
            try {
                for (final OutputStream out : templateOut) out.close();
                for (final OutputStream out : sampleBarcodeOut) out.close();
                for (final OutputStream out : molecularBarcodeOut) out.close();
            } catch (IOException ioe) {
                throw new RuntimeIOException(ioe);
            }
        }
    }


    /**
     * Binary codec used to encode/decode `ClusterData`
     */
    private static class ClusterDataCodec implements SortingCollection.Codec<ClusterData> {
        private final BinaryCodec binaryCodec = new BinaryCodec();

        @Override
        public void setOutputStream(OutputStream os) {
            this.binaryCodec.setOutputStream(os);
        }

        @Override
        public void setInputStream(InputStream is) {
            this.binaryCodec.setInputStream(is);
        }

        @Override
        public void encode(ClusterData clusterData) {
            this.binaryCodec.writeInt(clusterData.getNumReads());
            this.binaryCodec.writeInt(clusterData.getLane());
            this.binaryCodec.writeInt(clusterData.getTile());
            this.binaryCodec.writeInt(clusterData.getX());
            this.binaryCodec.writeInt(clusterData.getY());
            this.binaryCodec.writeBoolean(clusterData.isPf());

            if (clusterData.getMatchedBarcode() != null) {
                this.binaryCodec.writeString(clusterData.getMatchedBarcode(), true, true);
            } else {
                this.binaryCodec.writeString("", true, true);
            }

            for (int i = 0; i < clusterData.getNumReads(); i++) {
                ReadData read = clusterData.getRead(i);
                byte[] bases = read.getBases();
                byte[] quals = read.getQualities();

                binaryCodec.writeInt(bases.length);
                binaryCodec.writeString(read.getReadType().name(), false, false);

                for (int j = 0; j < bases.length; j++) {
                    this.binaryCodec.writeByte(bases[j]);
                    this.binaryCodec.writeByte(quals[j]);
                }
            }
        }

        @Override
        public ClusterData decode() {
            int numReads;
            try { numReads = this.binaryCodec.readInt(); }
            catch (final RuntimeEOFException e) { return null; }

            ReadData[] readData = new ReadData[numReads];
            ClusterData clusterData = new ClusterData(readData);
            clusterData.setLane(this.binaryCodec.readInt());
            clusterData.setTile(this.binaryCodec.readInt());
            clusterData.setX(this.binaryCodec.readInt());
            clusterData.setY(this.binaryCodec.readInt());
            clusterData.setPf(this.binaryCodec.readBoolean());
            String matchedBarcode = this.binaryCodec.readLengthAndString(true);
            if (matchedBarcode.length() == 0) {
                clusterData.setMatchedBarcode(null);
            } else {
                clusterData.setMatchedBarcode(matchedBarcode);
            }

            for (int i = 0; i < numReads; i++) {
                ReadData read = new ReadData();
                int numBases = this.binaryCodec.readInt();
                read.setReadType(ReadType.valueOf(this.binaryCodec.readString(1)));

                byte[] bases = new byte[numBases];
                byte[] quals = new byte[numBases];

                for (int j = 0; j < numBases; j++) {
                    bases[j] = this.binaryCodec.readByte();
                    quals[j] = this.binaryCodec.readByte();
                }

                read.setBases(bases);
                read.setQualities(quals);

                readData[i] = read;
            }

            return clusterData;
        }

        @Override
        public SortingCollection.Codec<ClusterData> clone() {
            return new ClusterDataCodec();
        }
    }

    /**
     * Comparator used for sorting `ClusterData` objects by read name.
     */
    private static class ClusterDataQueryNameComparator implements Comparator<ClusterData> {
        private final ReadNameEncoder readNameEncoder;

        /**
         * Creates a `ClusterData` comparator used to sort by query name given an read name encoder.
         *
         * @param readNameEncoder The read name encoder used to generate the read names for the clusters.
         */
        public ClusterDataQueryNameComparator(ReadNameEncoder readNameEncoder) {
            this.readNameEncoder = readNameEncoder;
        }

        @Override
        public int compare(ClusterData o1, ClusterData o2) {
            return readNameEncoder.generateShortName(o1).compareTo(readNameEncoder.generateShortName(o2));
        }
    }
}
