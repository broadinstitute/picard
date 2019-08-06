/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

package picard.sam.SamErrorMetric;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.reference.SamLocusAndReferenceIterator;
import htsjdk.samtools.reference.SamLocusAndReferenceIterator.SAMLocusAndReference;
import htsjdk.samtools.util.*;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Program to collect error metrics on bases stratified in various ways.
 */

@CommandLineProgramProperties(
        summary = "Program to collect error metrics on bases stratified in various ways.\n" +
                "<p>" +
                "Sequencing errors come in different 'flavors'. For example, some occur during sequencing while " +
                "others happen during library construction, prior to the sequencing. They may be correlated with " +
                "various aspect of the sequencing experiment: position in the read, base context, length of insert and so on.\n " +
                "<p>" +
                "This program collects two different kinds of error metrics (one which attempts to distinguish between pre- and " +
                "post- sequencer errors, and on which doesn't) and a collation of 'stratifiers' " +
                "each of which assigns bases into various bins. The stratifiers can be used together to generate a composite " +
                "stratification. " +
                "<p>" +
                "For example:" +
                "<p>" +
                "The BASE_QUALITY stratifier will place bases in bins according to their declared base quality. " +
                "The READ_ORDINALITY stratifier will place bases in one of two bins depending on whether their read " +
                "is 'first' or 'second'. One could generate a composite stratifier BASE_QUALITY:READ_ORDINALITY which will " +
                "do both stratifications as the same time. \n" +
                "<p>" +
                "The resulting metric file will be named according to a provided prefix and a suffix which is generated " +
                " automatically according to the error metric. " +
                "The tool can collect multiple metrics in a single pass and there should be hardly any " +
                "performance loss when specifying multiple metrics at the same time; the default includes a " +
                "large collection of metrics. \n" +
                "<p>" +
                "To estimate the error rate the tool assumes that all differences from the reference are errors. For this to be " +
                "a reasonable assumption the tool needs to know the sites at which the sample is actually polymorphic and a " +
                "confidence interval where the user is relatively certain that the polymorphic sites are known and accurate. " +
                "These two inputs are provided as a VCF and INTERVALS. The program will only process sites that are in the " +
                "intersection of the interval lists in the INTERVALS argument as long as they are not polymorphic in the " +
                "VCF.\n\n",
        oneLineSummary = "Program to collect error metrics on bases stratified in various ways.",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CollectSamErrorMetrics extends CommandLineProgram {

    /**
     * Classifies whether the given event is a match, insertion, or deletion. This is deliberately not expressed
     * as CIGAR operators, since there is no knowledge of the CIGAR string at the time that this is determined.
     */
    public enum BaseOperation {
        Match,
        Insertion,
        Deletion,
    }

    private static final int MAX_DIRECTIVES = ReadBaseStratification.Stratifier.values().length + 1;
    private static final Log log = Log.getInstance(CollectSamErrorMetrics.class);

    private final List<SamErrorReadFilter> samErrorReadFilters = new ArrayList<>();
    private SAMFileWriter filterOutputWriter;
    private final HashSet<String> outputReadOccurrence = new HashSet<>();

    // =====================================================================

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Base name for output files. Actual file names will be " +
            "generated from the basename and suffixes from the ERROR and STRATIFIER by adding a '.' and then " +
            "error_by_stratifier[_and_stratifier]* where 'error' is ERROR's extension, and 'stratifier' is STRATIFIER's suffix. " +
            "For example, an ERROR_METRIC of ERROR:BASE_QUALITY:GC_CONTENT will produce an extension '.error_by_base_quality_and_gc'. " +
            "The suffixes can be found in the documentation for ERROR_VALUE and SUFFIX_VALUE.")
    public File OUTPUT;

    @Argument(doc = "Errors to collect in the form of \"ERROR(:STRATIFIER)*\". " +
            "To see the values available for ERROR and STRATIFIER look at the documentation for the arguments " +
            "ERROR_VALUE and STRATIFIER_VALUE.")
    public List<String> ERROR_METRICS = CollectionUtil.makeList(
            "ERROR",
            "ERROR:BASE_QUALITY",
            "ERROR:INSERT_LENGTH",
            "ERROR:GC_CONTENT",
            "ERROR:READ_DIRECTION",
            "ERROR:PAIR_ORIENTATION",
            "ERROR:HOMOPOLYMER",
            "ERROR:BINNED_HOMOPOLYMER",
            "ERROR:CYCLE",
            "ERROR:INSERTIONS_IN_READ",
            "ERROR:DELETIONS_IN_READ",
            "ERROR:INDELS_IN_READ",
            "ERROR:READ_ORDINALITY",
            "ERROR:READ_ORDINALITY:CYCLE",
            "ERROR:READ_ORDINALITY:HOMOPOLYMER",
            "ERROR:READ_ORDINALITY:GC_CONTENT",
            "ERROR:READ_ORDINALITY:PRE_DINUC",
            "ERROR:MAPPING_QUALITY",
            "ERROR:READ_GROUP",
            "ERROR:MISMATCHES_IN_READ",
            "ERROR:ONE_BASE_PADDED_CONTEXT",
            "OVERLAPPING_ERROR",
            "OVERLAPPING_ERROR:BASE_QUALITY",
            "OVERLAPPING_ERROR:INSERT_LENGTH",
            "OVERLAPPING_ERROR:READ_ORDINALITY",
            "OVERLAPPING_ERROR:READ_ORDINALITY:CYCLE",
            "OVERLAPPING_ERROR:READ_ORDINALITY:HOMOPOLYMER",
            "OVERLAPPING_ERROR:READ_ORDINALITY:GC_CONTENT",
            "INDEL_ERROR",
            "INDEL_ERROR:INDEL_LENGTH"
    );

    @Argument(doc = "A fake argument used to show the options of ERROR (in ERROR_METRICS).", optional = true)
    public ErrorType ERROR_VALUE;

    @Argument(doc = "A fake argument used to show the options of STRATIFIER (in ERROR_METRICS).", optional = true)
    public ReadBaseStratification.Stratifier STRATIFIER_VALUE;

    @Argument(shortName = "V", doc = "VCF of known variation for sample. program will skip over polymorphic sites in this VCF and " +
            "avoid collecting data on these loci.")
    public File VCF;

    @Argument(shortName = "L", doc = "Region(s) to limit analysis to. Supported formats are VCF or interval_list. Will intersect inputs if multiple are given. ", optional = true)
    public List<File> INTERVALS;

    @Argument(shortName = "F", doc = "Filters for writing reads to output which match certain criteria in order to " +
            "being able to analyze them subsequently. Within each file, the first line defines the name of the " +
            "filter, which determines the name of the output file which lists reads that match this filter's criteria. " +
            "Each subsequent line represents one criterion. Each criterion is represented by 3 tab-separated fields: " +
            "\"suffix(TAB)datatype(TAB)operator(TAB)value\", where \"suffix\" is the suffix of the stratifier used for " +
            "this criterion. \"datatype\" can be any of \"boolean, int\" and \"operator\" can be one of " +
            "\"=, !=, <, <=, >, >=\" (depending on the datatype). These criteria are conjunctive, i.e. all criteria must " +
            "be met for the read to be included in the output. Multiple files can be provided which are disjuncitve to each " +
            "other (i.e. reads are included if they meet either of the filters).", optional = true)
    public List<File> FILTERS;

    @Argument(shortName = StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "Minimum mapping quality to include read.")
    public int MIN_MAPPING_Q = 20;

    @Argument(shortName = "BQ", doc = "Minimum base quality to include base.")
    public int MIN_BASE_Q = 20;

    @Argument(shortName = "PE", doc = "The prior error, in phred-scale (used for calculating empirical error rates).",
            optional = true)
    public int PRIOR_Q = 30;

    @Argument(shortName = "MAX", doc = "Maximum number of loci to process (or unlimited if 0).", optional = true)
    public long MAX_LOCI;

    @Argument(shortName = "LH", doc = "Shortest homopolymer which is considered long.  Used by the BINNED_HOMOPOLYMER stratifier.", optional = true)
    public int LONG_HOMOPOLYMER = 6;

    @Argument(shortName = "P", doc = "The probability of selecting a locus for analysis (for downsampling).", optional = true)
    public double PROBABILITY = 1;

    @Argument(
            fullName = "progressStepInterval",
            doc = "The interval between which progress will be displayed.",
            optional = true
    )
    public int progressStepInterval = 100000;

    // =====================================================================

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected String[] customCommandLineValidation() {
        List<String> errors = new ArrayList<>();

        if (ERROR_VALUE != null) {
            errors.add("ERROR_VALUE is a fake argument that is only there to show what are the different Error aggregation options. Please use it within the ERROR_METRICS argument.");
        }

        if (STRATIFIER_VALUE != null) {
            errors.add("STRATIFIER_VALUE is a fake argument that is only there to show what are the different Stratification options. Please use it within the STRATIFIER_VALUE argument.");
        }

        if (MIN_MAPPING_Q < 0) {
            errors.add("MIN_MAPPING_Q must be non-negative. found value: " + MIN_MAPPING_Q);
        }

        if (MIN_BASE_Q < 0) {
            errors.add("MIN_BASE_Q must be non-negative. found value: " + MIN_BASE_Q);
        }
        if (PRIOR_Q < 0) {
            errors.add("PRIOR_Q must be 2 or more. found value: " + PRIOR_Q);
        }

        if (MAX_LOCI < 0) {
            errors.add("MAX_LOCI must be non-negative. found value: " + MAX_LOCI);
        }

        if (LONG_HOMOPOLYMER < 0) {
            errors.add("LONG_HOMOPOLYMER must be non-negative. found value: " + LONG_HOMOPOLYMER);
        }

        if (PROBABILITY < 0 || PROBABILITY > 1) {
            errors.add("PROBABILITY must be between 0 and 1. found value: " + PROBABILITY);
        }

        final String[] superValidation = super.customCommandLineValidation();
        if (superValidation != null) {
            errors.addAll(Arrays.asList(superValidation));
        }
        if (!errors.isEmpty()) {
            return errors.toArray(new String[0]);
        } else {
            return null;
        }
    }

    @Override
    protected int doWork() {

        final Random random = new Random(42);

        final ProgressLogger progressLogger = new ProgressLogger(log, progressStepInterval);
        long nTotalLoci = 0;
        long nSkippedLoci = 0;
        long nProcessedLoci = 0;

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFilesAreReadable(INTERVALS);

        final Collection<BaseErrorAggregation> aggregatorList = getAggregatorList();
        // Open up the input resources:
        try (
                final SamReader sam = SamReaderFactory.makeDefault()
                        .referenceSequence(REFERENCE_SEQUENCE)
                        .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                        .open(INPUT);
                final ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
                final VCFFileReader vcfFileReader = new VCFFileReader(VCF, true)
        ) {

            // Make sure we can query our file:
            if (!vcfFileReader.isQueryable()) {
                throw new PicardException("Cannot query VCF File!  VCF Files must be queryable!");
            }

            final SAMSequenceDictionary sequenceDictionary = referenceSequenceFileWalker.getSequenceDictionary();
            if (sam.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new PicardException("Input BAM must be sorted by coordinate");
            }

            sequenceDictionary.assertSameDictionary(sam.getFileHeader().getSequenceDictionary());

            final IntervalList regionOfInterest = getIntervals(sequenceDictionary);

            log.info("Reading output read filters from files");
            readFiltersFromFiles();

            log.info("Getting SamLocusIterator");

            final SamLocusIterator samLocusIterator = new SamLocusIterator(sam, regionOfInterest);

            // We want to know about indels:
            samLocusIterator.setIncludeIndels(true);

            // Now go through and compare information at each of these sites
            samLocusIterator.setEmitUncoveredLoci(false);
            samLocusIterator.setMappingQualityScoreCutoff(MIN_MAPPING_Q);
            samLocusIterator.setQualityScoreCutoff(MIN_BASE_Q);

            log.info("Using " + aggregatorList.size() + " aggregators.");

            aggregatorList.forEach(la ->
                    IOUtil.assertFileIsWritable(new File(OUTPUT + la.getSuffix())));

            // iterate over loci
            log.info("Starting iteration over loci");

            final SamLocusAndReferenceIterator iterator = new SamLocusAndReferenceIterator(referenceSequenceFileWalker, samLocusIterator);

            if (samErrorReadFilters.size() > 0) {
                initializeFilterOutput(sam.getFileHeader());
            }

            // This hasNext() call has side-effects. It loads up the index and makes sure that the iterator is really ready for
            // action. Calling this allows for the logging to be more accurate.
            iterator.hasNext();
            log.info("Really starting iteration now.");

            for (final SAMLocusAndReference info : iterator) {

                if (random.nextDouble() > PROBABILITY) {
                    continue;
                }
                nTotalLoci++;

                // while there is a next (non-filtered) variant and it is before the locus, advance the pointer.
                if (checkLocus(vcfFileReader, info.getLocus(), sequenceDictionary)) {
                    log.debug("Locus does not overlap any variants: " + locusToInterval(info.getLocus(), sequenceDictionary));
                    nSkippedLoci++;
                    continue;
                }

                addLocusBases(aggregatorList, info);

                nProcessedLoci++;
                progressLogger.record(info.getLocus().getSequenceName(), info.getLocus().getPosition());

                if (MAX_LOCI != 0 && nProcessedLoci >= MAX_LOCI) {
                    log.warn("Early stopping due to having processed MAX_LOCI loci.");
                    break;
                }
            }

            log.info("Writing filtered reads");

            finalizeFilterOutput();

        } catch (final IOException e) {
            log.error(e, "A problem occurred:");
            return 1;
        }
        log.info("Iteration complete, generating metric files");

        aggregatorList.forEach(this::writeMetricsFileForAggregator);

        log.info(String.format("Examined %d loci, Processed %d loci, Skipped %d loci.\n" +
                "Computation took %d seconds.", nTotalLoci, nProcessedLoci, nSkippedLoci, progressLogger.getElapsedSeconds()));

        return 0;
    }

    /**
     * Initializes the writer to generate a BAM file for reads that match filter criteria and should be written to output
     *
     * @param samFileHeader The header to be used for the output BAM
     */
    private void initializeFilterOutput(final SAMFileHeader samFileHeader) {
        filterOutputWriter = new SAMFileWriterFactory()
                .setCompressionLevel(2)
                .makeWriter(samFileHeader, false, new File(OUTPUT + ".outputreads.bam"), REFERENCE_SEQUENCE);
    }

    /**
     * Closes the BAM file writer for the filtered output reads. In addition to that, for each filter specified in the
     * FILTER argument, a separate file is created which lists all reads that matched that specific filter's criteria.
     * Within each file, each line consists of the unique read ID and the number of occurrences that this read has
     * matched this filter's criteria. The reads in the output BAM file are annotated with this unique ID in the
     * CO:Z:uniquefilterid tag and can be queried that way.
     */
    private void finalizeFilterOutput() {
        if (filterOutputWriter == null) {
            return;
        }

        filterOutputWriter.close();

        for (final SamErrorReadFilter filter : samErrorReadFilters) {
            try (final BufferedWriter out = IOUtil.openFileForBufferedWriting(new File(OUTPUT + "." + filter.getName() + ".filteredreads"))) {
                for (final Map.Entry<String, Integer> entry : filter.getFilteredReads().entrySet()) {
                    out.write(entry.getKey() + "\t" + entry.getValue() + "\n");
                }
            } catch (final IOException e) {
                throw new SAMException("Could not write output files.", e);
            }
        }
    }

    /**
     * Given a "full" BaseCalculator, get all of its metrics and put them in the appropriate Metrics file.
     *
     * @param locusAggregator a BaseCalculator that has been "loaded up" with bases.
     */
    private void writeMetricsFileForAggregator(final BaseErrorAggregation locusAggregator) {
        final MetricsFile<ErrorMetric, Integer> file = getMetricsFile();

        ErrorMetric.setPriorError(QualityUtil.getErrorProbabilityFromPhredScore(PRIOR_Q));

        for (final ErrorMetric metric : locusAggregator.getMetrics()) {
            metric.calculateDerivedFields();
            file.addMetric(metric);
        }

        file.write(new File(OUTPUT + "." + locusAggregator.getSuffix()));
    }

    /**
     * If there's a variant at the locus return true, otherwise false.
     * <p>
     * HAS SIDE EFFECTS!!! Queries the vcfFileReader
     *
     * @param vcfFileReader      a {@link VCFFileReader} to query for the given locus
     * @param locus              a {@link Locus} at which to examine the variants
     * @param sequenceDictionary a dictionary with which to compare the Locatable to the Locus...
     * @return true if there's a variant over the locus, false otherwise.
     */
    private static boolean checkLocus(final VCFFileReader vcfFileReader, final Locus locus, final SAMSequenceDictionary sequenceDictionary) {

        boolean overlaps = false;
        final Interval queryInterval = locusToInterval(locus, sequenceDictionary);

        if (queryInterval != null) {

            try (final CloseableIterator<VariantContext> vcfIterator = vcfFileReader.query(queryInterval)) {

                overlaps = true;

                boolean allFiltered = true;

                while (vcfIterator.hasNext()) {
                    final VariantContext vcf = vcfIterator.next();
                    if (vcf.isFiltered()) {
                        continue;
                    }
                    allFiltered = false;
                }

                if (allFiltered) {
                    overlaps = false;
                }
            }
        }

        return overlaps;
    }

    /**
     * Converts the given locus into an interval using the given sequenceDictionary.
     *
     * @param locus              The {@link Locus} to convert.
     * @param sequenceDictionary The {@link SAMSequenceDictionary} to use to convert the given {@code locus}.
     * @return An {@link Interval} representing the given {@code locus} or {@code null} if it cannot be converted.
     */
    private static Interval locusToInterval(final Locus locus, final SAMSequenceDictionary sequenceDictionary) {
        final SAMSequenceRecord samSequenceRecord = sequenceDictionary.getSequence(locus.getSequenceIndex());
        if (samSequenceRecord == null) {
            return null;
        }

        return new Interval(samSequenceRecord.getSequenceName(), locus.getPosition(), locus.getPosition());
    }

    /**
     * Map of previously seen deletion records, associated with the locus they have been last seen
     */
    private final HashMap<SAMRecord, SamLocusIterator.LocusInfo> previouslySeenDeletions = new HashMap<>();

    /**
     * Current locus for taking care of deleting deletion records from the above map
     */
    private SamLocusIterator.LocusInfo currentLocus = null;

    /**
     * Checks if the same record has been seen at the previous locus already, thereby determining
     * whether or not a deletion has already been processed. Note that calling this method will
     * signal that the record is processed at this location.
     *
     * @param deletionRaO The RecordAndOffset to be checked
     * @param locusInfo   The LocusInfo to determine the current locus
     * @return True, if the record has been seen at the previous locus. False, if it has not yet been seen.
     */
    @VisibleForTesting
    protected boolean processDeletionLocus(final SamLocusIterator.RecordAndOffset deletionRaO, final SamLocusIterator.LocusInfo locusInfo) {
        if (currentLocus == null) {
            currentLocus = locusInfo;
        }

        // Check if we have moved to a new locus
        else if (!currentLocus.withinDistanceOf(locusInfo, 0)) {
            // If yes, remove all entries that have not been seen in the previous locus
            currentLocus = locusInfo;
            previouslySeenDeletions.entrySet().removeIf(entry -> !entry.getValue().withinDistanceOf(currentLocus, 1));
        }
        if (previouslySeenDeletions.containsKey(deletionRaO.getRecord())) {
            previouslySeenDeletions.put(deletionRaO.getRecord(), currentLocus);
            return true;
        }
        previouslySeenDeletions.put(deletionRaO.getRecord(), currentLocus);
        return false;
    }

    /**
     * Method for applying filters to stratifiers recursively, in case they are made up of multiple stratifiers
     *
     * @param stratifier Stratifier(s) to apply filters to. Can also be a CollectionStratifier or PairStratifier
     * @param rao        The RecordAndOffset for applying the filter
     * @param info       The SAMLocusAndReference for applying the filter
     * @param operation  Classifies whether the current event is a match, insertion, or deletion
     */
    private void applyFiltersIfFilterable(final ReadBaseStratification.RecordAndOffsetStratifier stratifier, final SamLocusIterator.RecordAndOffset rao, final SAMLocusAndReference info, final CollectSamErrorMetrics.BaseOperation operation) {
        if (stratifier instanceof ReadBaseStratification.FilterableRecordAndOffsetStratifier) {
            ((ReadBaseStratification.FilterableRecordAndOffsetStratifier) stratifier).stratifyAndApplyFilters(rao, info, operation, samErrorReadFilters);
        } else if (stratifier instanceof ReadBaseStratification.CollectionStratifier) {
            applyFiltersIfFilterable(((ReadBaseStratification.CollectionStratifier) stratifier).getStratifier(), rao, info, operation);
        } else if (stratifier instanceof ReadBaseStratification.PairStratifier) {
            applyFiltersIfFilterable(((ReadBaseStratification.PairStratifier) stratifier).getFirstStratifier(), rao, info, operation);
            applyFiltersIfFilterable(((ReadBaseStratification.PairStratifier) stratifier).getSecondStratifier(), rao, info, operation);
        }
    }

    /**
     * Stratifies the current RecordAndOffset and applies the filters to it. In case isDeletionRecord is true, the record
     * is checked if this deletion has already been processed, as it will be populated for each locus in the reference
     *
     * @param aggregatorList The aggregators to add the bases to
     * @param rao            The ReadAndOffset object
     * @param info           The SAMLocusAndReference object
     * @param operation      Classifies whether the current event is a match, insertion, or deletion
     */
    private void addAndFilterRecordAndOffset(final Collection<BaseErrorAggregation> aggregatorList, final SamLocusIterator.RecordAndOffset rao, final SAMLocusAndReference info, final BaseOperation operation) {
        // If deletion has been processed already, skip it

        if (operation == BaseOperation.Deletion && processDeletionLocus(rao, info.getLocus()))
            return;

        samErrorReadFilters.forEach(SamErrorReadFilter::reset);
        for (final BaseErrorAggregation aggregation : aggregatorList) {
            Object stratus = aggregation.addBase(rao, info, operation);

            applyFiltersIfFilterable(aggregation.getStratifier(), rao, info, operation);
        }

        for (final SamErrorReadFilter filter : samErrorReadFilters) {
            if (filter.isSatisfied()) {
                // Add the read to the filter for it to keep track of how many reads it matched
                filter.addReadById(getUniqueReadId(rao.getRecord()));

                // And add it to the output list
                addOutputRead(rao.getRecord());
            }
        }
    }

    /**
     * Iterate over the different records in the locus and add bases to aggregators
     */
    private void addLocusBases(final Collection<BaseErrorAggregation> aggregatorList, final SAMLocusAndReference info) {
        // Matching bases
        for (final SamLocusIterator.RecordAndOffset rao : info.getRecordAndOffsets()) {
            addAndFilterRecordAndOffset(aggregatorList, rao, info, BaseOperation.Match);
        }

        // Deleted bases
        for (final SamLocusIterator.RecordAndOffset deletionRao : info.getLocus().getDeletedInRecord()) {
            addAndFilterRecordAndOffset(aggregatorList, deletionRao, info, BaseOperation.Deletion);
        }

        // Inserted bases
        for (final SamLocusIterator.RecordAndOffset insertionRao : info.getLocus().getInsertedInRecord()) {
            addAndFilterRecordAndOffset(aggregatorList, insertionRao, info, BaseOperation.Insertion);
        }
    }

    /**
     * Gets a unique identifier for a read. In BAM files, this identifier is guaranteed to be unique, however, in SAM
     * files, this identifier is not available and the read name is used instead. This is likely to cause collisions
     * and as a consequence not all matching reads might be written to the output. To avoid this, consider using a BAM file.
     *
     * @param read The read to calculate the unique Id for
     * @return A String representing the position of the read in the BAM file, thereby acting as a unique identifier.
     * For SAM files, the read name is used instead.
     */
    private String getUniqueReadId(final SAMRecord read) {
        if (read.getFileSource() == null || read.getFileSource().getFilePointer() == null) {
            log.warn("There is no supported way to generate a unique read ID for SAM files. Instead, the " +
                    "read name is used, which is likely to cause collisions. As a consequence, not all reads matching " +
                    "the filter criteria might be written to the output. To avoid this, consider using a BAM file.");
            return read.getReadName();
        }
        return String.valueOf(((BAMFileSpan) read.getFileSource().getFilePointer()).getFirstOffset());
    }

    /**
     * Whenever a read matches filter criteria, write it to the output BAM file and remember that it has been written,
     * so we don't need to output it again. The unique identifier for that read is written to the SAM "CO" tag
     * (prefixed with "uniquefilterid:"), which can be used to search for reads specified in each filter's output file
     * (see {@link #finalizeFilterOutput()}).
     *
     * @param read The read that matched filter criteria and should be written to the output BAM file
     */
    private void addOutputRead(final SAMRecord read) {
        if (filterOutputWriter == null)
            return;

        String readId = getUniqueReadId(read);
        if (!outputReadOccurrence.contains(readId)) {
            // If it hasn't been seen yet, write that read to the output BAM file
            outputReadOccurrence.add(readId);
            read.setAttribute("CO", "uniquefilterid:" + readId);
            filterOutputWriter.addAlignment(read);
        }
    }

    /**
     * Parses the list of "directives" and creates a collection of {@link BaseErrorAggregation}s based on them.
     *
     * @return Collection of {@link BaseErrorAggregation}
     * @throws IllegalArgumentException if directives cannot be parsed or objects cannot be created.
     */

    private Collection<BaseErrorAggregation> getAggregatorList() {
        final List<BaseErrorAggregation> aggregatorList = new ArrayList<>();
        Set<String> suffixes = new HashSet<>();

        ReadBaseStratification.setLongHomopolymer(LONG_HOMOPOLYMER);
        for (final String directive : ERROR_METRICS) {
            final BaseErrorAggregation aggregator;
            aggregator = parseDirective(directive);
            aggregatorList.add(aggregator);
            if (!suffixes.add(aggregator.getSuffix())) {
                throw new IllegalArgumentException(String.format("Duplicated suffix (%s) found in aggregator %s.",
                        aggregator.getSuffix(), aggregator.getClass()));
            }
        }
        return aggregatorList;
    }

    /**
     * Interprets intervals from the input INTERVALS, if there's a file with that name, it opens the file, otherwise
     * it tries to parse it, checks that their dictionaries all agree with the input SequenceDictionary and returns the
     * intersection of all the lists.
     *
     * @param sequenceDictionary a {@link SAMSequenceDictionary} to parse the intervals against.
     * @return the intersection of the interval lists pointed to by the input parameter INTERVALS
     */
    private IntervalList getIntervals(final SAMSequenceDictionary sequenceDictionary) {

        IntervalList regionOfInterest = null;

        for (final File intervalListFile : INTERVALS) {

            if (!intervalListFile.exists()) {
                throw new IllegalArgumentException("Input file " + intervalListFile + " doesn't seem to exist. ");
            }

            log.info("Reading IntervalList ", intervalListFile, ".");
            final IntervalList temp = IntervalList.fromFile(intervalListFile);
            sequenceDictionary.assertSameDictionary(temp.getHeader().getSequenceDictionary());

            if (regionOfInterest == null) {
                regionOfInterest = temp;
            } else {
                log.info("Intersecting interval_list: ", intervalListFile, ".");
                regionOfInterest = IntervalList.intersection(regionOfInterest, temp);
            }
        }
        return regionOfInterest;
    }

    /**
     * Reads SamErrorReadFilters from the files specified as arguments. See argument documentation for detailed
     * format description
     */
    private void readFiltersFromFiles() {
        for (final File filterFile : FILTERS) {
            samErrorReadFilters.add(SamErrorReadFilter.fromFile(filterFile));
        }
    }

    /**
     * Compares a VariantContext to a Locus providing information regarding possible overlap, or relative location
     *
     * @param dictionary     The {@link SAMSequenceDictionary} to use for ordering the sequences
     * @param variantContext the {@link VariantContext} to compare
     * @param locus          the {@link Locus} to compare
     * @return negative if variantContext comes before locus (with no overlap)
     * zero if variantContext and locus overlap
     * positive if variantContext comes after locus (with no overlap)
     * <p/>
     * if variantContext and locus are in the same contig the return value will be the number of bases apart they are,
     * otherwise it will be MIN_INT/MAX_INT
     */
    public static int CompareVariantContextToLocus(final SAMSequenceDictionary dictionary, final VariantContext variantContext, final Locus locus) {

        log.debug("Comparing variant (" + variantContext.toStringWithoutGenotypes() + ") to locus (" + locus.toString() + ")");

        final int indexDiff = dictionary.getSequenceIndex(variantContext.getContig()) - locus.getSequenceIndex();
        if (indexDiff != 0) {
            return indexDiff < 0 ? Integer.MIN_VALUE : Integer.MAX_VALUE;
        }

        //same SequenceIndex, can compare by genomic position
        if (locus.getPosition() < variantContext.getStart())
            return variantContext.getStart() - locus.getPosition();
        if (locus.getPosition() > variantContext.getEnd())
            return variantContext.getEnd() - locus.getPosition();
        return 0;
    }

    /**
     * Parses a "Directive" of the form "ERROR,STRATIFIER,STRATIFIER...etc." into a {@link BaseErrorAggregation} consisting of the appropriate
     * {@link BaseCalculator} and the {@link ReadBaseStratification.CollectionStratifier CollectionStratifier} constructed from the various
     * individual stratifiers.
     * The conversion from string to object is performed by the enums
     * {@link ErrorType Errors} and {@link ReadBaseStratification.Stratifier Stratifier}.
     *
     * @param stratifierDirective The string directive describing the error type and collection of stratifiers to use
     * @return The appropriate {@link BaseErrorAggregation}.
     */
    protected static BaseErrorAggregation parseDirective(final String stratifierDirective) {
        final String[] directiveUnits = new String[MAX_DIRECTIVES + 1];
        final char directiveSeparator = ':';
        final int numberOfTerms = ParsingUtils.split(stratifierDirective, directiveUnits, directiveSeparator, false);

        if (numberOfTerms > MAX_DIRECTIVES) {
            throw new IllegalArgumentException(String.format("Cannot parse more than the number of different stratifiers plus one (%d) terms in a single directive. (What are you trying to do?)", MAX_DIRECTIVES));
        }
        if (numberOfTerms == 0) {
            throw new IllegalArgumentException("Found no directives at all. Cannot process.");
        }

        // make a linkedList due to removal and addition operations below
        final List<ReadBaseStratification.RecordAndOffsetStratifier<?>> stratifiers = Arrays.stream(directiveUnits, 1, numberOfTerms)
                .map(String::trim)
                .map(ReadBaseStratification.Stratifier::valueOf)
                .map(ReadBaseStratification.Stratifier::makeStratifier)
                .collect(Collectors.toList());

        final ReadBaseStratification.RecordAndOffsetStratifier jointStratifier;

        if (stratifiers.isEmpty()) {
            jointStratifier = ReadBaseStratification.nonStratifier;
        } else {
            jointStratifier = new ReadBaseStratification.CollectionStratifier(stratifiers);
        }

        // build an error supplier from the first term
        final Supplier<? extends BaseCalculator> supplier =
                ErrorType.valueOf(directiveUnits[0].trim()).getErrorSupplier();

        // return the aggregator made from the stratifier and the error supplier
        return new BaseErrorAggregation<>(supplier, jointStratifier);
    }
}

