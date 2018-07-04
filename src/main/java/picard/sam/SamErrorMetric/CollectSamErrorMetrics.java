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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
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
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
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
                "The tool cal collect multiple metrics in a single pass and there should be hardly any " +
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
public class CollectSamErrorMetrics extends CommandLineProgram {
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
            "OVERLAPPING_ERROR:READ_ORDINALITY:GC_CONTENT");

    @Argument(doc = "A fake argument used to show the options of ERROR (in ERROR_METRICS).", optional = true)
    public ErrorType ERROR_VALUE;

    @Argument(doc = "A fake argument used to show the options of STRATIFIER (in ERROR_METRICS).", optional = true)
    public ReadBaseStratification.Stratifier STRATIFIER_VALUE;

    @Argument(shortName = "V", doc = "VCF of known variation for sample. program will skip over polymorphic sites in this VCF and " +
            "avoid collecting data on these loci.")
    public File VCF;

    @Argument(shortName = "L", doc = "Region(s) to limit analysis to. Supported formats are VCF or interval_list. Will intersect inputs if multiple are given. ", optional = true)
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

    @Argument(shortName = "P", doc = "The probability of selecting a locus for analysis (for downsampling).", optional = true)
    public double PROBABILITY = 1;

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

    private static final int MAX_DIRECTIVES = ReadBaseStratification.Stratifier.values().length + 1;
    private static final Log log = Log.getInstance(CollectSamErrorMetrics.class);

    @Override
    protected int doWork() {

        final Random random = new Random(42);

        final ProgressLogger progressLogger = new ProgressLogger(log, 100000);
        long nTotalLoci = 0;
        long nSkippedLoci = 0;
        long nProcessedLoci = 0;

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFilesAreReadable(INTERVALS);

        final Collection<BaseErrorAggregation> aggregatorList = getAggregatorList();
        // Open up the input resources:
        try (
                final SamReader sam = SamReaderFactory.makeDefault().open(INPUT);
                final ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
                final PeekableIterator<VariantContext> vcfIterator = new PeekableIterator<>(
                        VCF == null ? Collections.emptyIterator() : new VCFFileReader(VCF, true).iterator())
        ) {

            final SAMSequenceDictionary sequenceDictionary = referenceSequenceFileWalker.getSequenceDictionary();
            if (sam.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new PicardException("Input BAM must be sorted by coordinate");
            }

            sequenceDictionary.assertSameDictionary(sam.getFileHeader().getSequenceDictionary());

            final IntervalList regionOfInterest = getIntervals(sequenceDictionary);

            log.info("Getting SamLocusIterator");

            final SamLocusIterator samLocusIterator = new SamLocusIterator(sam, regionOfInterest);

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
                if (advanceIteratorAndCheckLocus(vcfIterator, info.getLocus(), sequenceDictionary)) {
                    log.debug(String.format("Skipping locus from VCF: %s", vcfIterator.peek().toStringWithoutGenotypes()));
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

        } catch (IOException e) {
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
     * Iterate over the different records in the locus and add bases to aggregators
     */
    private void addLocusBases(final Collection<BaseErrorAggregation> aggregatorList, final SAMLocusAndReference info) {
        info.getRecordAndOffsets()
                .forEach(rao -> aggregatorList
                        .forEach(l -> l.addBase(rao, info)));
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
     * Reads Interprets intervals from the input INTERVALS, if there's a file with that name, it opens the file, otherwise
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

