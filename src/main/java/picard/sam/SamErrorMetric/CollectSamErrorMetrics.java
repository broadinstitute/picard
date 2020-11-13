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
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.reference.SamLocusAndReferenceIterator;
import htsjdk.samtools.reference.SamLocusAndReferenceIterator.SAMLocusAndReference;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locus;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SamLocusIterator;
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
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

    // =====================================================================

    private static final int MAX_DIRECTIVES = ReadBaseStratification.Stratifier.values().length + 1;
    private static final Log log = Log.getInstance(CollectSamErrorMetrics.class);

    // =====================================================================

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public String INPUT;

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
            "INDEL_ERROR"
    );

    @Argument(doc = "A fake argument used to show the options of ERROR (in ERROR_METRICS).", optional = true)
    public ErrorType ERROR_VALUE;

    @Argument(doc = "A fake argument used to show the options of STRATIFIER (in ERROR_METRICS).", optional = true)
    public ReadBaseStratification.Stratifier STRATIFIER_VALUE;

    @Argument(shortName = "V", doc = "VCF of known variation for sample. program will skip over polymorphic sites in this VCF and " +
            "avoid collecting data on these loci.")
    public String VCF;

    @Argument(shortName = "L", doc = "Region(s) to limit analysis to. Supported formats are VCF or interval_list. Will *intersect* inputs if multiple are given. " +
            "When this argument is supplied, the VCF provided must be *indexed*.", optional = true)
    public List<File> INTERVALS;

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

    @Argument(shortName = "LBS", doc = "Size of location bins. Used by the FLOWCELL_X and FLOWCELL_Y stratifiers", optional = true)
    public int LOCATION_BIN_SIZE = 2500;

    @Argument(shortName = "P", doc = "The probability of selecting a locus for analysis (for downsampling).", optional = true)
    public double PROBABILITY = 1;

    @Argument(
            fullName = "PROGRESS_STEP_INTERVAL",
            doc = "The interval between which progress will be displayed.",
            optional = true
    )
    public int PROGRESS_STEP_INTERVAL = 100000;

    @Argument(
            fullName = "INTERVAL_ITERATOR",
            doc = "Iterate through the file assuming it consists of a pre-created subset interval of the full genome.  " +
                    "This enables fast processing of files with reads at disperate parts of the genome.  " +
                    "Requires that the provided VCF file is indexed. ",
            optional = true
    )
    public boolean INTERVAL_ITERATOR = false;

    // =====================================================================

    /** Random object from which to pull pseudo-random numbers.  Initialized in {@link #initializeAggregationState()}.*/
    private Random random;

    /** Aggregator list to store errors.  Initialized in {@link #initializeAggregationState()}.*/
    private Collection<BaseErrorAggregation> aggregatorList;

    /** Logger with which to keep the user apprised of our progress.  Initialized in {@link #initializeAggregationState()}.*/
    private ProgressLogger progressLogger;

    /** Count of the total number of loci visited.  Initialized in {@link #initializeAggregationState()}.*/
    private long nTotalLoci;

    /** Count of the number of skipped loci.  Initialized in {@link #initializeAggregationState()}.*/
    private long nSkippedLoci;

    /** Count of the number of processed loci.  Initialized in {@link #initializeAggregationState()}.*/
    private long nProcessedLoci;

    /** The {@link SAMSequenceDictionary} for the reference, variants, and reads being processed. */
    private SAMSequenceDictionary sequenceDictionary;

    /** A {@link VCFFileReader} to read in variants when running in {@link #INTERVAL_ITERATOR} mode. */
    private VCFFileReader vcfFileReader;

    /** A {@link PeekableIterator<VariantContext>} to read in variants when running in default mode. */
    private PeekableIterator<VariantContext> vcfIterator;

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

    /**
     * Initialize the source for {@link VariantContext} information.
     * If running in {@link #INTERVAL_ITERATOR} mode, will initialize the {@link #vcfFileReader}.
     * Otherwise, will initialize the {@link #vcfIterator}.
     */
    private void initializeVcfDataSource() throws IOException {
        if ( INTERVAL_ITERATOR ) {
            vcfFileReader = new VCFFileReader(IOUtil.getPath(VCF), false);
            // Make sure we can query our file for interval mode:
            if (!vcfFileReader.isQueryable()) {
                throw new PicardException("Cannot query VCF File!  VCF Files must be queryable!  Please index input VCF and re-run.");
            }
        }
        else {
            vcfIterator = new PeekableIterator<>(
                    VCF == null ? Collections.emptyIterator() : new VCFFileReader(IOUtil.getPath(VCF), false).iterator());
        }
    }

    /**
     * Check if the given locus overlaps a variant from our variant source.
     *
     * If running in {@link #INTERVAL_ITERATOR} mode:
     * This will query the input variant file for overlapping variants, rather than just iterating through them.
     * It will be a little slower than the normal iteration, and much faster if a subset of a genome is used.
     *
     * If running in default mode:
     * Process our input data by iterating through input variants and reads.
     * Useful for going through whole genomes / exomes, but slow for sets of reads not beginning at the start of the
     * given variant file.
     * This slowness is because {@link CollectSamErrorMetrics} by default iterates through both the reads and the
     * variants sequentially from start to finish.  In the case your region of interest begins much later in the genome
     * (e.g. chr20:145000) then you will have to wait for the variant file to be traversed from the start to your
     * locus of interest.
     *
     * NOTE: This method HAS SIDE EFFECTS for both the {@link #vcfFileReader} and {@link #vcfIterator}.
     *
     * @return {@code true} if a variant overlaps the given locus.  {@code false} otherwise.
     */
    private boolean checkLocusForVariantOverlap(final SamLocusIterator.LocusInfo locusInfo) {
        final boolean returnValue;
        if (INTERVAL_ITERATOR) {
            returnValue = checkLocus(vcfFileReader, locusInfo);
            if ( returnValue ) {
                log.debug("Locus overlaps a known variant: " + locusInfo);
            }
        }
        else {
            returnValue = advanceIteratorAndCheckLocus(vcfIterator, locusInfo, sequenceDictionary);
            if ( returnValue ) {
                log.debug(String.format("Locus overlaps a known variant from VCF: %s -> %s", locusInfo.toString(),
                        vcfIterator.peek().toStringWithoutGenotypes()));
            }
        }
        return returnValue;
    }

    /**
     * Close the variant data source.
     * If running in {@link #INTERVAL_ITERATOR} mode will close the {@link #vcfFileReader};
     * If running in default mode will close the {@link #vcfIterator};
     */
    private void closeVcfDataSource() {
        if (INTERVAL_ITERATOR) {
            vcfFileReader.close();
        }
        else {
            vcfIterator.close();
        }
    }

    private int processData() {
        try (
                final SamReader sam = SamReaderFactory.makeDefault()
                        .referenceSequence(REFERENCE_SEQUENCE)
                        .open(IOUtil.getPath(INPUT));
                final ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE)
        ) {
            // Initialize our variants:
            initializeVcfDataSource();

            // Go through our reads and variants now:
            final SamLocusAndReferenceIterator iterator = createSamLocusAndReferenceIterator(sam, referenceSequenceFileWalker);

            log.info("Really starting iteration now.");
            for (final SAMLocusAndReference info : iterator) {
                if (random.nextDouble() > PROBABILITY) {
                    continue;
                }
                nTotalLoci++;

                // Check if a variant overlaps the current locus:
                if ( checkLocusForVariantOverlap(info.getLocus()) ) {
                    log.debug("Skipping overlapping locus.");
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

        } catch (final IOException e) {
            log.error(e, "A problem occurred:", e.getMessage());
            return 1;
        }
        finally {
            // Close our data source:
            closeVcfDataSource();
        }

        return 0;
    }

    /**
     * Creates the {@link SamLocusAndReferenceIterator} object to use to traverse the input files.
     * Also initializes the {@link #sequenceDictionary} and performs some checks on the output
     * files in {@link #aggregatorList} to make sure we can write our output.
     * @param sam an open {@link SamReader} from which to pull reads.
     * @param referenceSequenceFileWalker an open {@link ReferenceSequenceFileWalker} from which to get reference sequence information.
     * @return An open {@link SamLocusAndReferenceIterator} ready to iterate.
     */
    private SamLocusAndReferenceIterator createSamLocusAndReferenceIterator(final SamReader sam, final ReferenceSequenceFileWalker referenceSequenceFileWalker) {

        sequenceDictionary = referenceSequenceFileWalker.getSequenceDictionary();
        if (sam.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException("Input BAM must be sorted by coordinate");
        }

        // Make sure our reference and reads have the same sequence dictionary:
        sequenceDictionary.assertSameDictionary(sam.getFileHeader().getSequenceDictionary());

        final IntervalList regionOfInterest = getIntervals(sequenceDictionary);

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

        // This hasNext() call has side-effects. It loads up the index and makes sure that
        // the iterator is really ready for
        // action. Calling this allows for the logging to be more accurate.
        iterator.hasNext();
        return iterator;
    }

    /**
     * Initializes the counts, random object, and iterators for this class.
     *
     * Specifically, the following fields are initialized:
     *
     *   random
     *   aggregatorList
     *   progressLogger
     *   nTotalLoci
     *   nSkippedLoci
     *   nProcessedLoci
     */
    private void initializeAggregationState() {
        random = new Random(42);

        aggregatorList = getAggregatorList();

        progressLogger = new ProgressLogger(log, PROGRESS_STEP_INTERVAL);
        nTotalLoci = 0;
        nSkippedLoci = 0;
        nProcessedLoci = 0;
    }

    @Override
    protected int doWork() {

        // Initialize our iteration:
        initializeAggregationState();

        // Make sure we can read our files:
        try {
            IOUtil.assertFileIsReadable(IOUtil.getPath(INPUT));
            IOUtil.assertFileIsReadable(IOUtil.getPath(VCF));
            IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
            IOUtil.assertFilesAreReadable(INTERVALS);
        }
        catch (final IOException e) {
            log.error(e, "A problem occurred:", e.getMessage());
            return 1;
        }

        // Process our data based on how we will be iterating:
        final int returnValue = processData();

        // Check if we had an error and if so, immediately return
        // (to preserve old functionality):
        if (returnValue == 0) {
            log.info("Iteration complete, generating metric files");

            aggregatorList.forEach(this::writeMetricsFileForAggregator);

            log.info(String.format("Examined %d loci, Processed %d loci, Skipped %d loci.\n" +
                    "Computation took %d seconds.", nTotalLoci, nProcessedLoci, nSkippedLoci, progressLogger.getElapsedSeconds()));
        }
        return returnValue;
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
     * Advance the iterator until at or ahead of locus. if there's a variant at the locus return true, otherwise false.
     * <p>
     * HAS SIDE EFFECTS!!! draws from vcfIterator!!!
     *
     * @param vcfIterator        a {@link VariantContext} iterator to advance (assumed sorted)
     * @param locus              a {@link Locus} at which to examine the variants
     * @param sequenceDictionary a dictionary with which to compare the Locatable to the Locus...
     * @return true if there's a variant over the locus, false otherwise.
     */
    private static boolean advanceIteratorAndCheckLocus(final PeekableIterator<VariantContext> vcfIterator, final Locus locus, final SAMSequenceDictionary sequenceDictionary) {
        while (vcfIterator.hasNext() && (vcfIterator.peek().isFiltered() ||
                CompareVariantContextToLocus(sequenceDictionary, vcfIterator.peek(), locus) < 0)) {
            vcfIterator.next();
        }
        return vcfIterator.hasNext() && CompareVariantContextToLocus(sequenceDictionary, vcfIterator.peek(), locus) == 0;
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
     * have the side effect of signaling that the record is processed at this location.
     *
     * @param deletionRao The RecordAndOffset to be checked
     * @param locusInfo   The LocusInfo to determine the current locus
     * @return True, if the record has been seen at the previous locus. False, if it has not yet been seen.
     */
    @VisibleForTesting
    protected boolean processDeletionLocus(final SamLocusIterator.RecordAndOffset deletionRao, final SamLocusIterator.LocusInfo locusInfo) {
        if (currentLocus == null) {
            currentLocus = locusInfo;
        }

        // Check if we have moved to a new locus
        else if (!currentLocus.withinDistanceOf(locusInfo, 0)) {
            // If yes, remove all entries that have not been seen in the previous locus
            currentLocus = locusInfo;
            previouslySeenDeletions.entrySet().removeIf(entry -> !entry.getValue().withinDistanceOf(currentLocus, 1));
        }
        if (previouslySeenDeletions.containsKey(deletionRao.getRecord())) {
            previouslySeenDeletions.put(deletionRao.getRecord(), currentLocus);
            return true;
        }
        previouslySeenDeletions.put(deletionRao.getRecord(), currentLocus);
        return false;
    }

    /**
     * Stratifies the current RecordAndOffset. In case isDeletionRecord is true, the record is checked for whether or not
     * this deletion has already been processed, as it will be populated for each locus in the reference
     *
     * @param aggregatorList The aggregators to add the bases to
     * @param rao            The ReadAndOffset object
     * @param info           The SAMLocusAndReference object
     */
    private void addRecordAndOffset(final Collection<BaseErrorAggregation> aggregatorList, final SamLocusIterator.RecordAndOffset rao, final SAMLocusAndReference info) {
        // If deletion has been processed already, skip it
        if (rao.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Deletion && processDeletionLocus(rao, info.getLocus()))
            return;

        for (final BaseErrorAggregation aggregation : aggregatorList) {
            aggregation.addBase(rao, info);
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

        final int indexDiff = dictionary.getSequenceIndex(variantContext.getContig()) - locus.getSequenceIndex();
        if (indexDiff != 0) {
            return indexDiff < 0 ? Integer.MIN_VALUE : Integer.MAX_VALUE;
        }

        //same SequenceIndex, can compare by genomic position
        if (locus.getPosition() < variantContext.getStart()) {
            return variantContext.getStart() - locus.getPosition();
        }
        if (locus.getPosition() > variantContext.getEnd()) {
            return variantContext.getEnd() - locus.getPosition();
        }
        return 0;
    }

    /**
     * If there's a variant at the locus return true, otherwise false.
     * <p>
     * HAS SIDE EFFECTS!!! Queries the vcfFileReader
     *
     * @param vcfFileReader      a {@link VCFFileReader} to query for the given locus
     * @param locusInfo          a {@link SamLocusIterator.LocusInfo} at which to examine the variants
     * @return true if there's a variant over the locus, false otherwise.
     */
    private static boolean checkLocus(final VCFFileReader vcfFileReader, final SamLocusIterator.LocusInfo locusInfo) {
        boolean overlaps = false;

        if (locusInfo != null) {
            try (final CloseableIterator<VariantContext> vcfIterator = vcfFileReader.query(locusInfo)) {

                while (vcfIterator.hasNext()) {
                    if (vcfIterator.next().isFiltered()) {
                        continue;
                    }
                    overlaps = true;
                    break;
                }
            }
        }
        return overlaps;
    }

    /**
     * Iterate over the different records in the locus and add bases to aggregators
     */
    private void addLocusBases(final Collection<BaseErrorAggregation> aggregatorList, final SAMLocusAndReference info) {
        // Matching bases
        for (final SamLocusIterator.RecordAndOffset rao : info.getRecordAndOffsets()) {
            addRecordAndOffset(aggregatorList, rao, info);
        }

        // Deleted bases
        for (final SamLocusIterator.RecordAndOffset deletionRao : info.getLocus().getDeletedInRecord()) {
            addRecordAndOffset(aggregatorList, deletionRao, info);
        }

        // Inserted bases
        for (final SamLocusIterator.RecordAndOffset insertionRao : info.getLocus().getInsertedInRecord()) {
            addRecordAndOffset(aggregatorList, insertionRao, info);
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
        ReadBaseStratification.setLocationBinSize(LOCATION_BIN_SIZE);
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

