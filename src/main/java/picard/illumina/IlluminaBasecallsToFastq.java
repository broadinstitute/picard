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

import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.Illumina;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Illumina;
import picard.fastq.Casava18ReadNameEncoder;
import picard.fastq.IlluminaReadNameEncoder;
import picard.fastq.ReadNameEncoder;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaFileUtil;
import picard.illumina.parser.ReadData;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.IlluminaUtil;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

@CommandLineProgramProperties(
        summary = IlluminaBasecallsToFastq.USAGE_SUMMARY + IlluminaBasecallsToFastq.USAGE_DETAILS,
        oneLineSummary = IlluminaBasecallsToFastq.USAGE_SUMMARY,
        programGroup = Illumina.class
)
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

    @Argument(doc = "Configure SortingCollections to store this many records before spilling to disk. For an indexed" +
            " run, each SortingCollection gets this value/number of indices.")
    public int MAX_READS_IN_RAM_PER_TILE = 1200000;

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

    private final Map<String, FastqRecordsWriter> sampleBarcodeFastqWriterMap = new HashMap<>();
    private ReadStructure readStructure;
    private BasecallsConverter<FastqRecordsForCluster> basecallsConverter;
    private static final Log log = Log.getInstance(IlluminaBasecallsToFastq.class);
    private final FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
    private ReadNameEncoder readNameEncoder;
    private static final Comparator<FastqRecordsForCluster> queryNameComparator = (r1, r2) -> SAMRecordQueryNameComparator.compareReadNames(r1.templateRecords[0].getReadHeader(),
            r2.templateRecords[0].getReadHeader());

    @Override
    protected int doWork() {
        initialize();
        basecallsConverter.doTileProcessing();
        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        final LinkedList<String> errors = new LinkedList<>();
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

        final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);
        readStructure = new ReadStructure(READ_STRUCTURE);
        if (MULTIPLEX_PARAMS != null) {
            IOUtil.assertFileIsReadable(MULTIPLEX_PARAMS);
        }
        final boolean demultiplex;
        if (OUTPUT_PREFIX != null) {
            sampleBarcodeFastqWriterMap.put(null, buildWriter(OUTPUT_PREFIX));
            demultiplex = false;
        } else {
            populateWritersFromMultiplexParams();
            demultiplex = true;
        }
        final int readsPerCluster = readStructure.templates.length() + readStructure.sampleBarcodes.length();
        if (IlluminaFileUtil.hasCbcls(BASECALLS_DIR, LANE)) {
            if (BARCODES_DIR == null) BARCODES_DIR = BASECALLS_DIR;
            basecallsConverter = new NewIlluminaBasecallsConverter<>(BASECALLS_DIR, BARCODES_DIR, LANE, readStructure,
                    sampleBarcodeFastqWriterMap, demultiplex, Math.max(1, MAX_READS_IN_RAM_PER_TILE / readsPerCluster),
                    TMP_DIR, NUM_PROCESSORS,
                    FIRST_TILE, TILE_LIMIT, queryNameComparator,
                    new FastqRecordsForClusterCodec(readStructure.templates.length(),
                            readStructure.sampleBarcodes.length(), readStructure.molecularBarcode.length()),
                    FastqRecordsForCluster.class, bclQualityEvaluationStrategy, IGNORE_UNEXPECTED_BARCODES);
        } else {
            basecallsConverter = new IlluminaBasecallsConverter<>(BASECALLS_DIR, BARCODES_DIR, LANE, readStructure,
                    sampleBarcodeFastqWriterMap, demultiplex, Math.max(1, MAX_READS_IN_RAM_PER_TILE / readsPerCluster), TMP_DIR, NUM_PROCESSORS,
                    FORCE_GC, FIRST_TILE, TILE_LIMIT, queryNameComparator,
                    new FastqRecordsForClusterCodec(readStructure.templates.length(),
                            readStructure.sampleBarcodes.length(), readStructure.molecularBarcode.length()), FastqRecordsForCluster.class, bclQualityEvaluationStrategy,
                    this.APPLY_EAMSS_FILTER, INCLUDE_NON_PF_READS, IGNORE_UNEXPECTED_BARCODES);
        }

        basecallsConverter.setConverter(
                new ClusterToFastqRecordsForClusterConverter(
                        basecallsConverter.getFactory().getOutputReadStructure()));

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

        for (final TabbedTextFileWithHeaderParser.Row row : libraryParamsParser) {
            List<String> sampleBarcodeValues = null;

            if (!sampleBarcodeColumnLabels.isEmpty()) {
                sampleBarcodeValues = new ArrayList<>();
                for (final String sampleBarcodeLabel : sampleBarcodeColumnLabels) {
                    sampleBarcodeValues.add(row.getField(sampleBarcodeLabel));
                }
            }

            final String key = (sampleBarcodeValues == null || sampleBarcodeValues.contains("N")) ? null : StringUtil.join("", sampleBarcodeValues);
            if (sampleBarcodeFastqWriterMap.containsKey(key)) {    //This will catch the case of having more than 1 line in a non-barcoded MULTIPLEX_PARAMS file
                throw new PicardException("Row for barcode " + key + " appears more than once in MULTIPLEX_PARAMS file " +
                        MULTIPLEX_PARAMS);
            }

            final FastqRecordsWriter writer = buildWriter(new File(row.getField("OUTPUT_PREFIX")));
            sampleBarcodeFastqWriterMap.put(key, writer);
        }
        if (sampleBarcodeFastqWriterMap.isEmpty()) {
            throw new PicardException("MULTIPLEX_PARAMS file " + MULTIPLEX_PARAMS + " does have any data rows.");
        }
        libraryParamsParser.close();
    }

    /**
     * @return FastqRecordsWriter that contains one or more FastqWriters (amount depends on read structure), all using
     * outputPrefix to determine the filename(s).
     */
    private FastqRecordsWriter buildWriter(final File outputPrefix) {
        final File outputDir = outputPrefix.getAbsoluteFile().getParentFile();
        IOUtil.assertDirectoryIsWritable(outputDir);
        final String prefixString = outputPrefix.getName();
        final String suffixString = COMPRESS_OUTPUTS ? "fastq.gz" : "fastq";
        final FastqWriter[] templateWriters = new FastqWriter[readStructure.templates.length()];
        final FastqWriter[] sampleBarcodeWriters = new FastqWriter[readStructure.sampleBarcodes.length()];
        final FastqWriter[] molecularBarcodeWriters = new FastqWriter[readStructure.molecularBarcode.length()];

        for (int i = 0; i < templateWriters.length; ++i) {
            final String filename = String.format("%s.%d.%s", prefixString, i + 1, suffixString);
            templateWriters[i] = fastqWriterFactory.newWriter(new File(outputDir, filename));
        }

        for (int i = 0; i < sampleBarcodeWriters.length; ++i) {
            final String filename = String.format("%s.barcode_%d.%s", prefixString, i + 1, suffixString);
            sampleBarcodeWriters[i] = fastqWriterFactory.newWriter(new File(outputDir, filename));
        }

        for (int i = 0; i < molecularBarcodeWriters.length; ++i) {
            final String filename = String.format("%s.index_%d.%s", prefixString, i + 1, suffixString);
            molecularBarcodeWriters[i] = fastqWriterFactory.newWriter(new File(outputDir, filename));
        }
        return new FastqRecordsWriter(templateWriters, sampleBarcodeWriters, molecularBarcodeWriters);
    }

    public static void main(final String[] args) {
        new IlluminaBasecallsToFastq().instanceMainWithExit(args);
    }

    /**
     * Container for various FastqWriters, one for each template read, one for each sample barcode read,
     * and one for each molecular barcode read.
     */
    private static final class FastqRecordsWriter implements BasecallsConverter.ConvertedClusterDataWriter<FastqRecordsForCluster> {
        final FastqWriter[] templateWriters;
        final FastqWriter[] sampleBarcodeWriters;
        final FastqWriter[] molecularBarcodeWriters;

        /**
         * @param templateWriters         Writers for template reads in order, e,g. 0th element is for template read 1.
         * @param sampleBarcodeWriters    Writers for sample barcode reads in order, e,g. 0th element is for sample barcode read 1.
         * @param molecularBarcodeWriters Writers for molecular barcode reads in order, e,g. 0th element is for molecualr barcode read 1.
         */
        private FastqRecordsWriter(final FastqWriter[] templateWriters, final FastqWriter[] sampleBarcodeWriters, final FastqWriter[] molecularBarcodeWriters) {
            this.templateWriters = templateWriters;
            this.sampleBarcodeWriters = sampleBarcodeWriters;
            this.molecularBarcodeWriters = molecularBarcodeWriters;
        }

        @Override
        public void write(final FastqRecordsForCluster records) {
            write(templateWriters, records.templateRecords);
            write(sampleBarcodeWriters, records.sampleBarcodeRecords);
            write(molecularBarcodeWriters, records.molecularBarcodeRecords);
        }

        private void write(final FastqWriter[] writers, final FastqRecord[] records) {
            for (int i = 0; i < writers.length; ++i) {
                writers[i].write(records[i]);
            }
        }

        @Override
        public void close() {
            for (final FastqWriter writer : templateWriters) {
                writer.close();
            }
            for (final FastqWriter writer : sampleBarcodeWriters) {
                writer.close();
            }
            for (final FastqWriter writer : molecularBarcodeWriters) {
                writer.close();
            }
        }
    }

    /**
     * Contains the results of transforming one cluster into the record(s) to be written to output file(s).
     */
    static class FastqRecordsForCluster {
        // These are accessed directly by converter and writer rather than through getters and setters.
        final FastqRecord[] templateRecords;
        final FastqRecord[] sampleBarcodeRecords;
        final FastqRecord[] molecularBarcodeRecords;

        FastqRecordsForCluster(final int numTemplates, final int numSampleBarcodes, final int numMolecularBarcodes) {
            templateRecords = new FastqRecord[numTemplates];
            sampleBarcodeRecords = new FastqRecord[numSampleBarcodes];
            molecularBarcodeRecords = new FastqRecord[numMolecularBarcodes];
        }
    }

    /**
     * Passed to IlluminaBaseCallsConverter to do the conversion from input format to output format.
     */
    class ClusterToFastqRecordsForClusterConverter
            implements IlluminaBasecallsConverter.ClusterDataConverter<FastqRecordsForCluster> {

        private final int[] templateIndices;
        private final int[] sampleBarcodeIndicies;
        private final int[] molecularBarcodeIndicies;

        ClusterToFastqRecordsForClusterConverter(final ReadStructure outputReadStructure) {
            this.templateIndices = outputReadStructure.templates.getIndices();
            this.sampleBarcodeIndicies = outputReadStructure.sampleBarcodes.getIndices();
            this.molecularBarcodeIndicies = outputReadStructure.molecularBarcode.getIndices();
        }

        @Override
        public FastqRecordsForCluster convertClusterToOutputRecord(final ClusterData cluster) {
            final FastqRecordsForCluster ret = new FastqRecordsForCluster(readStructure.templates.length(), readStructure.sampleBarcodes.length(), readStructure.molecularBarcode.length());
            final boolean appendTemplateNumberSuffix = ret.templateRecords.length > 1;
            final boolean appendMolecularBarcodeNumber = ret.molecularBarcodeRecords.length > 1;

            makeFastqRecords(ret.templateRecords, templateIndices, cluster, appendTemplateNumberSuffix);
            makeFastqRecords(ret.sampleBarcodeRecords, sampleBarcodeIndicies, cluster, false);
            makeFastqRecords(ret.molecularBarcodeRecords, molecularBarcodeIndicies, cluster, appendMolecularBarcodeNumber);

            return ret;
        }

        private void makeFastqRecords(final FastqRecord[] recs, final int[] indices,
                                      final ClusterData cluster, final boolean appendReadNumberSuffix) {
            for (short i = 0; i < indices.length; ++i) {
                final ReadData readData = cluster.getRead(indices[i]);
                final String readBases = StringUtil.bytesToString(readData.getBases()).replace('.', 'N');
                final String readName = readNameEncoder.generateReadName(cluster, appendReadNumberSuffix ? i + 1 : null);
                recs[i] = new FastqRecord(
                        readName,
                        readBases,
                        null,
                        SAMUtils.phredToFastq(readData.getQualities())
                );
            }
        }
    }

    /**
     * Codec passed to IlluminaBasecallsConverter for use in SortingCollections of output records.
     */
    static class FastqRecordsForClusterCodec implements SortingCollection.Codec<FastqRecordsForCluster> {
        private final int numTemplates;
        private final int numSampleBarcodes;
        private final int numMolecularBarcodes;

        private BasicFastqWriter writer = null;
        private FastqReader reader = null;

        FastqRecordsForClusterCodec(final int numTemplates, final int numSampleBarcodes, final int numMolecularBarcodes) {
            this.numTemplates = numTemplates;
            this.numSampleBarcodes = numSampleBarcodes;
            this.numMolecularBarcodes = numMolecularBarcodes;
        }

        @Override
        public void setOutputStream(final OutputStream os) {
            writer = new BasicFastqWriter(new PrintStream(os));
        }

        @Override
        public void setInputStream(final InputStream is) {
            reader = new FastqReader(new BufferedReader(new InputStreamReader(is)));
        }

        //TODO: add tests to encode and decode
        @Override
        public void encode(final FastqRecordsForCluster val) {
            if (numTemplates != val.templateRecords.length) throw new IllegalStateException();
            if (numSampleBarcodes != val.sampleBarcodeRecords.length) throw new IllegalStateException();
            encodeArray(val.templateRecords);
            encodeArray(val.sampleBarcodeRecords);
            encodeArray(val.molecularBarcodeRecords);
            writer.flush();
        }

        private void encodeArray(final FastqRecord[] recs) {
            for (final FastqRecord rec : recs) {
                writer.write(rec);
            }
        }

        @Override
        public FastqRecordsForCluster decode() {
            if (!reader.hasNext()) return null;
            final FastqRecordsForCluster ret = new FastqRecordsForCluster(numTemplates, numSampleBarcodes, numMolecularBarcodes);
            decodeArray(ret.templateRecords);
            decodeArray(ret.sampleBarcodeRecords);
            decodeArray(ret.molecularBarcodeRecords);
            return ret;
        }

        private void decodeArray(final FastqRecord[] recs) {
            for (int i = 0; i < recs.length; ++i) {
                recs[i] = reader.next();
            }
        }

        @Override
        public SortingCollection.Codec<FastqRecordsForCluster> clone() {
            return new FastqRecordsForClusterCodec(numTemplates, numSampleBarcodes, numMolecularBarcodes);
        }
    }
}
