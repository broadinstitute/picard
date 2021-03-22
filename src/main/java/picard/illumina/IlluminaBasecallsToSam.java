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

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.BaseCallingProgramGroup;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.util.AdapterPair;
import picard.util.IlluminaUtil;
import picard.util.IlluminaUtil.IlluminaAdapterPair;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.*;

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
 * previous paragraph.  As writes complete, the TileReadAggregator re-evaluates the state of reads/writes and may queue
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
@CommandLineProgramProperties(
        summary = IlluminaBasecallsToSam.USAGE_SUMMARY + IlluminaBasecallsToSam.USAGE_DETAILS,
        oneLineSummary = IlluminaBasecallsToSam.USAGE_SUMMARY,
        programGroup = BaseCallingProgramGroup.class
)
@DocumentedFeature
public class IlluminaBasecallsToSam extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Transforms raw Illumina sequencing data into an unmapped SAM or BAM file.";
    static final String USAGE_DETAILS = "<p>The IlluminaBaseCallsToSam program collects, demultiplexes, and sorts reads across all " +
            "of the tiles of a lane via barcode to produce an unmapped SAM/BAM file.  An unmapped BAM file is often referred to as a uBAM.  " +
            "All barcode, sample, and library data is provided in the LIBRARY_PARAMS file.  Note, this LIBRARY_PARAMS file " +
            "should be formatted according to the specifications indicated below.  The following is an example of a properly" +
            " formatted LIBRARY_PARAMS file:</p>" +
            "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\n" +
            "AAAAAAAA\tSA_AAAAAAAA.bam\tSA_AAAAAAAA\tLN_AAAAAAAA\n" +
            "AAAAGAAG\tSA_AAAAGAAG.bam\tSA_AAAAGAAG\tLN_AAAAGAAG\n" +
            "AACAATGG\tSA_AACAATGG.bam\tSA_AACAATGG\tLN_AACAATGG\n" +
            "N\tSA_non_indexed.bam\tSA_non_indexed\tLN_NNNNNNNN\n " +
            "" +
            "<p>The BARCODES_DIR file is produced by the " +
            "<a href='http://broadinstitute.github.io/picard/command-line-overview.html#ExtractIlluminaBarcodes'>ExtractIlluminaBarcodes</a> " +
            "tool for each lane of a flow cell.</p>  " +

            "<h4>Usage example:</h4>" +
            "<pre>" +
            "" +
            "java -jar picard.jar IlluminaBasecallsToSam \\<br />" +
            "      BASECALLS_DIR=/BaseCalls/ \\<br />" +
            "      LANE=001 \\<br />" +
            "      READ_STRUCTURE=25T8B25T \\<br />" +
            "      RUN_BARCODE=run15 \\<br />" +
            "      IGNORE_UNEXPECTED_BARCODES=true \\<br />" +
            "      LIBRARY_PARAMS=library.params " +
            "</pre>" +
            "<hr />";


    // The following attributes define the command-line arguments

    public static final String USAGE = "Generate a SAM or BAM file from data in an Illumina basecalls output directory";

    @Argument(doc = "The basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(doc = "The barcodes directory with _barcode.txt files (generated by ExtractIlluminaBarcodes). If not set, use BASECALLS_DIR. ", shortName = "BCD", optional = true)
    public File BARCODES_DIR;

    @Argument(doc = "Lane number. ", shortName = StandardOptionDefinitions.LANE_SHORT_NAME)
    public Integer LANE;

    @Argument(doc = "Deprecated (use LIBRARY_PARAMS).  The output SAM or BAM file. Format is determined by extension.",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            mutex = {"BARCODE_PARAMS", "LIBRARY_PARAMS"})
    public File OUTPUT;

    @Argument(doc = "The barcode of the run.  Prefixed to read names.")
    public String RUN_BARCODE;

    @Argument(doc = "Deprecated (use LIBRARY_PARAMS).  The name of the sequenced sample",
            shortName = StandardOptionDefinitions.SAMPLE_ALIAS_SHORT_NAME,
            mutex = {"BARCODE_PARAMS", "LIBRARY_PARAMS"})
    public String SAMPLE_ALIAS;

    @Argument(doc = "ID used to link RG header record with RG tag in SAM record.  " +
            "If these are unique in SAM files that get merged, merge performance is better.  " +
            "If not specified, READ_GROUP_ID will be set to <first 5 chars of RUN_BARCODE>.<LANE> .",
            shortName = StandardOptionDefinitions.READ_GROUP_ID_SHORT_NAME, optional = true)
    public String READ_GROUP_ID;

    @Argument(doc = "Deprecated (use LIBRARY_PARAMS).  The name of the sequenced library",
            shortName = StandardOptionDefinitions.LIBRARY_NAME_SHORT_NAME,
            optional = true,
            mutex = {"BARCODE_PARAMS", "LIBRARY_PARAMS"})
    public String LIBRARY_NAME;

    @Argument(doc = "The name of the sequencing center that produced the reads.  Used to set the @RG->CN header tag.")
    public String SEQUENCING_CENTER;

    @Argument(doc = "The start date of the run.", optional = true)
    public Date RUN_START_DATE;

    @Argument(doc = "The name of the sequencing technology that produced the read.", optional = true)
    public String PLATFORM = "ILLUMINA";

    @Argument(doc = "Whether to include the barcode information in the @RG->BC header tag. Defaults to false until included in the SAM spec.")
    public boolean INCLUDE_BC_IN_RG_TAG = false;

    @Argument(doc = ReadStructure.PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Argument(doc = "Deprecated (use LIBRARY_PARAMS).  Tab-separated file for creating all output BAMs for barcoded run " +
            "with single IlluminaBasecallsToSam invocation.  Columns are BARCODE, OUTPUT, SAMPLE_ALIAS, and " +
            "LIBRARY_NAME.  Row with BARCODE=N is used to specify a file for no barcode match",
            mutex = {"OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "LIBRARY_PARAMS"})
    public File BARCODE_PARAMS;

    @Argument(doc = "Tab-separated file for creating all output BAMs for a lane with single IlluminaBasecallsToSam " +
            "invocation.  The columns are OUTPUT, SAMPLE_ALIAS, and LIBRARY_NAME, BARCODE_1, BARCODE_2 ... BARCODE_X " +
            "where X = number of barcodes per cluster (optional).  Row with BARCODE_1 set to 'N' is used to specify a file " +
            "for no barcode match.  You may also provide any 2 letter RG header attributes (excluding PU, CN, PL, and" +
            " DT)  as columns in this file and the values for those columns will be inserted into the RG tag for the" +
            " BAM file created for a given row.",
            mutex = {"OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "BARCODE_PARAMS"})
    public File LIBRARY_PARAMS;

    @Argument(doc = "Which adapters to look for in the read.")
    public List<IlluminaAdapterPair> ADAPTERS_TO_CHECK = new ArrayList<>(
            Arrays.asList(IlluminaAdapterPair.INDEXED,
                    IlluminaAdapterPair.DUAL_INDEXED,
                    IlluminaAdapterPair.NEXTERA_V2,
                    IlluminaAdapterPair.FLUIDIGM));

    @Argument(doc = "For specifying adapters other than standard Illumina", optional = true)
    public String FIVE_PRIME_ADAPTER;

    @Argument(doc = "For specifying adapters other than standard Illumina", optional = true)
    public String THREE_PRIME_ADAPTER;

    @Argument(doc = "The number of threads to run in parallel. If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0, then the number of cores used will" +
            " be the number available on the machine less NUM_PROCESSORS.")
    public Integer NUM_PROCESSORS = 0;

    @Argument(doc = "If set, this is the first tile to be processed (used for debugging).  Note that tiles are not processed" +
            " in numerical order.",
            mutex = "PROCESS_SINGLE_TILE",
            optional = true)
    public Integer FIRST_TILE;

    @Argument(doc = "If set, process no more than this many tiles (used for debugging).", optional = true)
    public Integer TILE_LIMIT;

    @Argument(doc = "If set, process only the tile number given and prepend the tile number to the output file name.",
            mutex = "FIRST_TILE",
            optional = true)
    public Integer PROCESS_SINGLE_TILE;

    @Argument(doc = "Apply EAMSS filtering to identify inappropriately quality scored bases towards the ends of reads" +
            " and convert their quality scores to Q2.")
    public boolean APPLY_EAMSS_FILTER = true;

    @Argument(doc = "Configure SortingCollections to store this many records before spilling to disk. For an indexed" +
            " run, each SortingCollection gets this value/number of indices. Deprecated: use `MAX_RECORDS_IN_RAM`")
    public int MAX_READS_IN_RAM_PER_TILE = -1;

    @Argument(doc = "The minimum quality (after transforming 0s to 1s) expected from reads.  If qualities are lower than this value, an error is thrown." +
            "The default of 2 is what the Illumina's spec describes as the minimum, but in practice the value has been observed lower.")
    public int MINIMUM_QUALITY = BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

    @Argument(doc = "Whether to include non-PF reads", shortName = "NONPF", optional = true)
    public boolean INCLUDE_NON_PF_READS = true;

    @Argument(doc = "Whether to ignore reads whose barcodes are not found in LIBRARY_PARAMS.  Useful when outputting " +
            "BAMs for only a subset of the barcodes in a lane.", shortName = "IGNORE_UNEXPECTED")
    public boolean IGNORE_UNEXPECTED_BARCODES = false;

    @Argument(doc = "The tag to use to store any molecular indexes.  If more than one molecular index is found, they will be concatenated and stored here.", optional = true)
    public String MOLECULAR_INDEX_TAG = "RX";

    @Argument(doc = "The tag to use to store any molecular index base qualities.  If more than one molecular index is found, their qualities will be concatenated and stored here " +
            "(.i.e. the number of \"M\" operators in the READ_STRUCTURE)", optional = true)
    public String MOLECULAR_INDEX_BASE_QUALITY_TAG = "QX";

    @Argument(doc = "The list of tags to store each molecular index.  The number of tags should match the number of molecular indexes.", optional = true)
    public List<String> TAG_PER_MOLECULAR_INDEX;

    @Argument(doc = "When should the sample barcode (as read by the sequencer) be placed on the reads in the BC tag?")
    public ClusterDataToSamConverter.PopulateBarcode BARCODE_POPULATION_STRATEGY = ClusterDataToSamConverter.PopulateBarcode.ORPHANS_ONLY;

    @Argument(doc = "Should the barcode quality be included when the sample barcode is included?")
    public boolean INCLUDE_BARCODE_QUALITY = false;

    @Argument(doc = "If true, the output records are sorted by read name. Otherwise they are unsorted.")
    public Boolean SORT = true;

    private Map<String, SAMFileWriterWrapper> barcodeSamWriterMap;
    private ReadStructure readStructure;
    private BasecallsConverter<SAMRecordsForCluster> basecallsConverter;
    private static final Log log = Log.getInstance(IlluminaBasecallsToSam.class);
    private final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);

    @Override
    protected int doWork() {
        initialize();
        try {
            basecallsConverter.processTilesAndWritePerSampleOutputs(barcodeSamWriterMap.keySet());
        } catch (IOException e) {
            throw new PicardException("Error converting basecalls to SAM.", e);
        }
        return 0;
    }

    /**
     * Prepares loggers, initiates garbage collection thread, parses arguments and initialized variables appropriately/
     */
    private void initialize() {
        if (OUTPUT != null) {
            IOUtil.assertFileIsWritable(OUTPUT);
        }

        if (LIBRARY_PARAMS != null) {
            IOUtil.assertFileIsReadable(LIBRARY_PARAMS);
        }

        if (OUTPUT != null) {
            barcodeSamWriterMap = new HashMap<>(1, 1.0f);
            barcodeSamWriterMap.put(null, buildSamFileWriter(OUTPUT, SAMPLE_ALIAS, LIBRARY_NAME, buildSamHeaderParameters(null), SORT));
        } else {
            populateWritersFromLibraryParams();
        }

        final int numOutputRecords = readStructure.templates.length();
        // Combine any adapters and custom adapter pairs from the command line into an array for use in clipping
        final List<AdapterPair> adapters = new ArrayList<>(ADAPTERS_TO_CHECK);

        if (FIVE_PRIME_ADAPTER != null && THREE_PRIME_ADAPTER != null) {
            adapters.add(new CustomAdapterPair(FIVE_PRIME_ADAPTER, THREE_PRIME_ADAPTER));
        }

        final boolean demultiplex = readStructure.hasSampleBarcode();
        BasecallsConverterBuilder<SAMRecordsForCluster> converterBuilder = new BasecallsConverterBuilder<>(BASECALLS_DIR, LANE, readStructure, barcodeSamWriterMap)
                .barcodesDir(BARCODES_DIR)
                .withDemultiplex(demultiplex)
                .numProcessors(NUM_PROCESSORS)
                .firstTile(FIRST_TILE)
                .tileLimit(TILE_LIMIT)
                .withApplyEamssFiltering(APPLY_EAMSS_FILTER)
                .withIncludeNonPfReads(INCLUDE_NON_PF_READS)
                .withIgnoreUnexpectedBarcodes(IGNORE_UNEXPECTED_BARCODES)
                .withBclQualityEvaluationStrategy(bclQualityEvaluationStrategy);

        if (SORT) {
            converterBuilder = converterBuilder
                    .withSorting(
                            new QueryNameComparator(),
                            new Codec(numOutputRecords),
                            SAMRecordsForCluster.class,
                            TMP_DIR);
        }

        basecallsConverter = converterBuilder.build();
        /*
         * Be sure to pass the outputReadStructure to ClusterDataToSamConverter, which reflects the structure of the output cluster
         * data which may be different from the input read structure (specifically if there are skips).
         */
        final ClusterDataToSamConverter converter = new ClusterDataToSamConverter(RUN_BARCODE, READ_GROUP_ID,
                basecallsConverter.getFactory().getOutputReadStructure(), adapters, BARCODE_POPULATION_STRATEGY, INCLUDE_BARCODE_QUALITY)
                .withMolecularIndexTag(MOLECULAR_INDEX_TAG)
                .withMolecularIndexQualityTag(MOLECULAR_INDEX_BASE_QUALITY_TAG)
                .withTagPerMolecularIndex(TAG_PER_MOLECULAR_INDEX);
        basecallsConverter.setConverter(converter);
        log.info("DONE_READING STRUCTURE IS " + readStructure.toString());
    }

    /**
     * Assert that expectedCols are present and return actualCols - expectedCols
     *
     * @param actualCols   The columns present in the LIBRARY_PARAMS file
     * @param expectedCols The columns that are REQUIRED
     * @return actualCols - expectedCols
     */
    private Set<String> findAndFilterExpectedColumns(final Set<String> actualCols, final Set<String> expectedCols) {
        final Set<String> missingColumns = new HashSet<>(expectedCols);
        missingColumns.removeAll(actualCols);

        if (!missingColumns.isEmpty()) {
            throw new PicardException(String.format(
                    "LIBRARY_PARAMS file %s is missing the following columns: %s.",
                    LIBRARY_PARAMS.getAbsolutePath(), StringUtil.join(", ", missingColumns
                    )));
        }

        final Set<String> remainingColumns = new HashSet<>(actualCols);
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

        if (!forbiddenHeaders.isEmpty()) {
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
     * the key to the map is the concatenation of all sampleBarcodes in order for the given line
     */
    private void populateWritersFromLibraryParams() {
        final TabbedTextFileWithHeaderParser libraryParamsParser = new TabbedTextFileWithHeaderParser(LIBRARY_PARAMS);

        final Set<String> expectedColumnLabels = CollectionUtil.makeSet("OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME");
        final List<String> barcodeColumnLabels = new ArrayList<>();
        if (readStructure.sampleBarcodes.length() == 1) {
            //For the single barcode read case, the barcode label name can either by BARCODE or BARCODE_1
            if (libraryParamsParser.hasColumn("BARCODE")) {
                barcodeColumnLabels.add("BARCODE");
            } else if (libraryParamsParser.hasColumn("BARCODE_1")) {
                barcodeColumnLabels.add("BARCODE_1");
            } else {
                throw new PicardException("LIBRARY_PARAMS(BARCODE_PARAMS) file " + LIBRARY_PARAMS + " does not have column BARCODE or BARCODE_1.");
            }
        } else {
            for (int i = 1; i <= readStructure.sampleBarcodes.length(); i++) {
                barcodeColumnLabels.add("BARCODE_" + i);
            }
        }

        expectedColumnLabels.addAll(barcodeColumnLabels);
        final Set<String> rgTagColumns = findAndFilterExpectedColumns(libraryParamsParser.columnLabels(), expectedColumnLabels);
        checkRgTagColumns(rgTagColumns);

        final List<TabbedTextFileWithHeaderParser.Row> rows = libraryParamsParser.iterator().toList();
        barcodeSamWriterMap = new HashMap<>(rows.size(), 1);

        for (final TabbedTextFileWithHeaderParser.Row row : rows) {
            List<String> barcodeValues = null;

            if (!barcodeColumnLabels.isEmpty()) {
                barcodeValues = new ArrayList<>();
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

            File outputFile = new File(row.getField("OUTPUT"));

            // If we are processing a single tile we want to append the tile number to the output file name. This is
            // done to avoid file overwrites if you are running tile based processing in parallel.
            if (PROCESS_SINGLE_TILE != null) {
                outputFile = new File(outputFile.getParentFile(),
                        PROCESS_SINGLE_TILE + "." + outputFile.getName());
            }

            final SAMFileWriterWrapper writer = buildSamFileWriter(outputFile,
                    row.getField("SAMPLE_ALIAS"), row.getField("LIBRARY_NAME"), samHeaderParams, SORT);
            barcodeSamWriterMap.put(key, writer);
        }
        if (barcodeSamWriterMap.isEmpty()) {
            throw new PicardException("LIBRARY_PARAMS(BARCODE_PARAMS) file " + LIBRARY_PARAMS + " does have any data rows.");
        }
        libraryParamsParser.close();
    }

    /**
     * Create the list of headers that will be added to the SAMFileHeader for a library with the given sampleBarcodes (or
     * the entire run if sampleBarcodes == NULL).  Note that any value that is null will NOT be added via buildSamFileWriter
     * but is placed in the map in order to be able to query the tags that we automatically add.
     *
     * @param barcodes The list of sampleBarcodes that uniquely identify the read group we are building parameters for
     * @return A Map of ReadGroupHeaderTags -> Values
     */
    private Map<String, String> buildSamHeaderParameters(final List<String> barcodes) {
        final Map<String, String> params = new LinkedHashMap<>();

        String platformUnit = RUN_BARCODE + "." + LANE;
        if (barcodes != null) {
            final String barcodeString = IlluminaUtil.barcodeSeqsToString(barcodes);
            platformUnit += "." + barcodeString;
            if (INCLUDE_BC_IN_RG_TAG) {
                params.put("BC", barcodeString);
            }
        }

        if (PLATFORM != null) {
            params.put(SAMReadGroupRecord.PLATFORM_TAG, PLATFORM);
        }

        params.put(SAMReadGroupRecord.PLATFORM_UNIT_TAG, platformUnit);
        if (SEQUENCING_CENTER != null) {
            params.put(SAMReadGroupRecord.SEQUENCING_CENTER_TAG, SEQUENCING_CENTER);
        }
        params.put(SAMReadGroupRecord.DATE_RUN_PRODUCED_TAG, RUN_START_DATE == null ? null : new Iso8601Date(RUN_START_DATE).toString());

        return params;
    }

    /**
     * Build a SamFileWriter that will write its contents to the output file.
     *
     * @param output           The file to which to write
     * @param sampleAlias      The sample alias set in the read group header
     * @param libraryName      The name of the library to which this read group belongs
     * @param headerParameters Header parameters that will be added to the RG header for this SamFile
     * @return A SAMFileWriter
     */
    private SAMFileWriterWrapper buildSamFileWriter(final File output, final String sampleAlias,
                                                    final String libraryName, final Map<String, String> headerParameters,
                                                    final boolean presorted) {
        IOUtil.assertFileIsWritable(output);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(READ_GROUP_ID);
        rg.setSample(sampleAlias);

        if (libraryName != null) rg.setLibrary(libraryName);
        for (final Map.Entry<String, String> tagNameToValue : headerParameters.entrySet()) {
            if (tagNameToValue.getValue() != null) {
                rg.setAttribute(tagNameToValue.getKey(), tagNameToValue.getValue());
            }
        }

        final SAMFileHeader header = new SAMFileHeader();

        if (presorted) {
            header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        } else {
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        }
        header.addReadGroup(rg);
        return new SAMFileWriterWrapper(new SAMFileWriterFactory().makeSAMOrBAMWriter(header, presorted, output));
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access args.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     * to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {

        if (NUM_PROCESSORS == 0) {
            NUM_PROCESSORS = Runtime.getRuntime().availableProcessors();
        } else if (NUM_PROCESSORS < 0) {
            NUM_PROCESSORS = Runtime.getRuntime().availableProcessors() + NUM_PROCESSORS;
        }

        if (BARCODE_PARAMS != null) {
            LIBRARY_PARAMS = BARCODE_PARAMS;
        }

        // Remove once deprecated parameter is deleted.
        if (MAX_READS_IN_RAM_PER_TILE != -1) {
            log.warn("Setting deprecated parameter `MAX_READS_IN_RAM_PER_TILE` use ` MAX_RECORDS_IN_RAM` instead");
            MAX_RECORDS_IN_RAM = MAX_READS_IN_RAM_PER_TILE;
        }

        final ArrayList<String> messages = new ArrayList<>();

        readStructure = new ReadStructure(READ_STRUCTURE);
        if (readStructure.hasSampleBarcode() && LIBRARY_PARAMS == null) {
            messages.add("BARCODE_PARAMS or LIBRARY_PARAMS is missing.  If READ_STRUCTURE contains a B (barcode)" +
                    " then either LIBRARY_PARAMS or BARCODE_PARAMS(deprecated) must be provided!");
        }

        if (READ_GROUP_ID == null) {
            READ_GROUP_ID = RUN_BARCODE.substring(0, Math.min(RUN_BARCODE.length(), 5)) + "." + LANE;
        }

        if (!TAG_PER_MOLECULAR_INDEX.isEmpty() && TAG_PER_MOLECULAR_INDEX.size() != readStructure.molecularBarcode.length()) {
            messages.add("The number of tags given in TAG_PER_MOLECULAR_INDEX does not match the number of molecular indexes in READ_STRUCTURE");
        }

        if ((FIVE_PRIME_ADAPTER == null) != (THREE_PRIME_ADAPTER == null)) {
            messages.add("THREE_PRIME_ADAPTER and FIVE_PRIME_ADAPTER must either both be null or both be set.");
        }

        // If we are processing a single tile we need to set TILE_LIMIT and FIRST_TILE for the underlying
        // basecalls converter.
        if (PROCESS_SINGLE_TILE != null) {
            TILE_LIMIT = 1;
            FIRST_TILE = PROCESS_SINGLE_TILE;
        }

        if (messages.isEmpty()) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }

    private static final class SAMFileWriterWrapper
            implements BasecallsConverter.ConvertedClusterDataWriter<SAMRecordsForCluster> {
        public final SAMFileWriter writer;

        private SAMFileWriterWrapper(final SAMFileWriter writer) {
            this.writer = writer;
        }

        @Override
        public void write(final SAMRecordsForCluster records) {
            for (final SAMRecord rec : records.records) {
                writer.addAlignment(rec);
            }
        }

        @Override
        public void close() {
            writer.close();
        }
    }

    static class SAMRecordsForCluster {
        final SAMRecord[] records;

        SAMRecordsForCluster(final int numRecords) {
            records = new SAMRecord[numRecords];
        }
    }

    static class QueryNameComparator implements Comparator<SAMRecordsForCluster> {
        private final SAMRecordQueryNameComparator comparator = new SAMRecordQueryNameComparator();

        @Override
        public int compare(final SAMRecordsForCluster s1, final SAMRecordsForCluster s2) {
            return comparator.compare(s1.records[0], s2.records[0]);
        }
    }

    static class Codec implements SortingCollection.Codec<SAMRecordsForCluster> {
        private final BAMRecordCodec bamCodec;
        private final int numRecords;

        Codec(final int numRecords, final BAMRecordCodec bamCodec) {
            this.numRecords = numRecords;
            this.bamCodec = bamCodec;
        }

        Codec(final int numRecords) {
            this(numRecords, new BAMRecordCodec(null));
        }

        @Override
        public void setOutputStream(final OutputStream os) {
            bamCodec.setOutputStream(os);
        }

        @Override
        public void setInputStream(final InputStream is) {
            bamCodec.setInputStream(is);
        }

        @Override
        public void encode(final SAMRecordsForCluster val) {
            if (val.records.length != numRecords) {
                throw new IllegalStateException(String.format("Expected number of clusters %d != actual %d",
                        numRecords, val.records.length));
            }
            for (final SAMRecord rec : val.records) {
                bamCodec.encode(rec);
            }
        }

        @Override
        public SAMRecordsForCluster decode() {
            final SAMRecord zerothRecord = bamCodec.decode();
            if (zerothRecord == null) return null;
            final SAMRecordsForCluster ret = new SAMRecordsForCluster(numRecords);
            ret.records[0] = zerothRecord;
            for (int i = 1; i < numRecords; ++i) {
                ret.records[i] = bamCodec.decode();
                if (ret.records[i] == null) {
                    throw new IllegalStateException(String.format("Expected to read %d records but read only %d", numRecords, i));
                }
            }
            return ret;
        }

        @Override
        public SortingCollection.Codec<SAMRecordsForCluster> clone() {
            return new Codec(numRecords, bamCodec.clone());
        }
    }
}
