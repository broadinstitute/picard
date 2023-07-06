/*
 * The MIT License
 *
 * Copyright (c) 2011-2016 The Broad Institute
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
 * 
 */

/**
 * $Id$
 */
package picard.sam;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.nio.PicardHtsPath;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <h3>Summary</h3>
 * Subsets a SAM file by either selecting or excluding certain reads
 * <p>
 * <h3>Details</h3>
 * Subsets a SAM or BAM file by either excluding or selecting reads as specified by FILTER.
 * Other parameters influence the behavior of the FILTER algorithm as described below.
 * <p>
 * <h3>Usage examples</h3>
 * <h4>Filter by queryname:</h4>
 * <pre>
 * java -jar picard.jar FilterSamReads \
 *       I=input.bam \
 *       O=output.bam \
 *       READ_LIST_FILE=read_names.txt \
 *       FILTER=includeReadList
 * </pre>
 * <h4>Filter by interval:</h4>
 * <pre>
 * java -jar picard.jar FilterSamReads \
 *       I=input.bam \
 *       O=output.bam \
 *       INTERVAL_LIST=regions.interval_list \
 *       FILTER=includePairedIntervals
 * </pre>
 * <h4>Filter reads having a (2-base or more) soft clip on the beginning of the read:</h4>
 * <pre>
 * cat <<EOF > script.js
 * // reads having a soft clip larger than 2 bases in start of read
 * function accept(rec) {
 *     if (rec.getReadUnmappedFlag()) return false;
 *     var cigar = rec.getCigar();
 *     if (cigar == null) return false;
 *     var ce = cigar.getCigarElement(0);
 *     return ce.getOperator().name() == "S" && ce.length() > 2;
 * }
 *
 * accept(record);
 * EOF
 *
 * java -jar picard.jar FilterSamReads \
 *       I=input.bam \
 *       O=output.bam \
 *       JAVASCRIPT_FILE=script.js \
 *       FILTER=includeJavascript
 * </pre>
 */
@CommandLineProgramProperties(
        summary = FilterSamReads.USAGE_SUMMARY + FilterSamReads.USAGE_DETAILS,
        oneLineSummary = FilterSamReads.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class FilterSamReads extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Subsets reads from a SAM/BAM/CRAM file by applying one of several filters.";
    static final String USAGE_DETAILS = "\nTakes a SAM/BAM/CRAM file and subsets it by either excluding or " +
            "only including certain reads such as aligned or unaligned reads, specific reads based on a list of reads names, " +
            "an interval list, by Tag Values (type Z / String values only), or using a JavaScript script.\n" +
            "<br />" +
            "<h3>Usage example:</h3>" +
            "<h4>Filter by queryname</h4>" +
            "<pre>" +
            "java -jar picard.jar FilterSamReads \\<br /> " +
            "      I=input.bam \\ <br /> " +
            "      O=output.bam \\ <br /> " +
            "      READ_LIST_FILE=read_names.txt \\ <br />" +
            "      FILTER=includeReadList" +
            "</pre> " +
            "<h4>Filter by interval</h4>" +
            "<pre>" +
            "java -jar picard.jar FilterSamReads \\ <br /> " +
            "      I=input.bam \\ <br /> " +
            "      O=output.bam \\ <br /> " +
            "      INTERVAL_LIST=regions.interval_list \\ <br/>" +
            "      FILTER=includePairedIntervals" +
            "</pre> " +
            "<h4>Filter by Tag Value (type Z / String values only)</h4>" +
            "<pre>" +
            "java -jar picard.jar FilterSamReads \\ <br /> " +
            "      I=input.bam \\ <br /> " +
            "      O=output.bam \\ <br /> " +
            "      TAG=CR \\ <br/>" +
            "      TAG_VALUE=TTTGTCATCTCGAGTA \\ <br/>" +
            "      FILTER=includeTagValues" +
            "</pre> " +
            "<h4>Filter reads having a soft clip on the beginning of the read larger than 2 bases with a JavaScript script</h4>" +
            "<pre>" +
            "cat <<EOF > script.js <br/>" +
            "/** reads having a soft clip larger than 2 bases in beginning of read*/ <br/>" +
            "function accept(rec) {   <br/>" +
            "    if (rec.getReadUnmappedFlag()) return false; <br/>" +
            "    var cigar = rec.getCigar(); <br/>" +
            "    if (cigar == null) return false; <br/>" +
            "    var ce = cigar.getCigarElement(0); <br/>" +
            "    return ce.getOperator().name() == \"S\" && ce.length() > 2; <br/>" +
            "} <br />" +
            "<br />" +
            "accept(record); <br/>" +
            "EOF <br/>" +
            "<br/>" +
            "java -jar picard.jar FilterSamReads \\ <br /> " +
            "      I=input.bam \\ <br /> " +
            "      O=output.bam \\ <br /> " +
            "      JAVASCRIPT_FILE=script.js \\ <br/>" +
            "      FILTER=includeJavascript" +
            "</pre> ";
    private static final Log log = Log.getInstance(FilterSamReads.class);

    @VisibleForTesting
    protected enum Filter implements CommandLineParser.ClpEnum {
        includeAligned("Output aligned reads only. INPUT SAM/BAM/CRAM must be in queryname SortOrder. (Note: first and second of paired reads must both be aligned to be included in OUTPUT.)"),
        excludeAligned("Output Unmapped reads only. INPUT SAM/BAM/CRAM must be in queryname SortOrder. (Note: first and second of pair must both be aligned to be excluded from OUTPUT.)"),
        includeReadList("Output reads with names contained in READ_LIST_FILE. See READ_LIST_FILE for more detail."),
        excludeReadList("Output reads with names *not* contained in READ_LIST_FILE. See READ_LIST_FILE for more detail."),
        includeJavascript("Output reads that have been accepted by the JAVASCRIPT_FILE script, that is, reads for which the value of the script is true. " +
                "See the JAVASCRIPT_FILE argument for more detail. "),
        includePairedIntervals("Output reads that overlap with an interval from INTERVAL_LIST (and their mate). INPUT must be coordinate sorted."),
        includeTagValues("Output reads that have a value of tag TAG that is contained in the values for TAG_VALUES"),
        excludeTagValues("Output reads that do not have a value of tag TAG that is contained in the values for TAG_VALUES");
       private final String description;

        Filter(final String description) {
            this.description = description;
        }

        @Override
        public String getHelpDoc() {
            return description;
        }
    }

    @Argument(doc = "The SAM/BAM/CRAM file that will be filtered.",
            shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public PicardHtsPath INPUT;

    @Argument(doc = "Which filter to use.")
    public Filter FILTER = null;

    @Argument(doc = "File containing reads that will be included in or excluded from the OUTPUT SAM/BAM/CRAM file, when using FILTER=includeReadList or FILTER=excludeReadList.",
            optional = true,
            shortName = "RLF")
    public File READ_LIST_FILE;

   @Argument(doc = "Interval List File containing intervals that will be included in the OUTPUT when using FILTER=includePairedIntervals",
            optional = true,
            shortName = "IL")
    public PicardHtsPath INTERVAL_LIST;

    @Argument(doc = "The tag to select from input SAM/BAM",
            optional = true,
            shortName = "T")
    public String TAG;

    @Argument(doc = "The tag value(s) to filter by",
            optional = true,
            shortName = "TV")
    public List<String> TAG_VALUE;

    @Argument(
            doc = "SortOrder of the OUTPUT file, otherwise use the SortOrder of the INPUT file.",
            optional = true,
            shortName = "SO")
    public SAMFileHeader.SortOrder SORT_ORDER;

    @Argument(doc = "SAM/BAM/CRAM file for resulting reads.",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public PicardHtsPath OUTPUT;

    @Argument(shortName = "JS",
            doc = "Filters the INPUT with a javascript expression using the java javascript-engine, when using FILTER=includeJavascript. "
                    + " The script puts the following variables in the script context: \n"
                    + " 'record' a SamRecord ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html ) and \n "
                    + " 'header' a SAMFileHeader ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html ).\n"
                    + " all the public members of SamRecord and SAMFileHeader are accessible. "
                    + "A record is accepted if the last value of the script evaluates to true.",
            optional = true)
    public File JAVASCRIPT_FILE = null; // tsato: update this too?

    /**
     * This method should be used only for testing needs
     * @param referenceSequence
     */
	void setReferenceSequence(File referenceSequence) {
        this.REFERENCE_SEQUENCE = referenceSequence;
    }

    @Argument(
            doc = "Create <OUTPUT>.reads file containing names of reads from INPUT and OUTPUT (for debugging purposes.)",
            optional = true)
    public boolean WRITE_READS_FILES = false;

    private void filterReads(final FilteringSamIterator filteringIterator) {

        // get OUTPUT header from INPUT and overwrite it if necessary
        final SAMFileHeader fileHeader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(INPUT.toPath());
        final SAMFileHeader.SortOrder inputSortOrder = fileHeader.getSortOrder();
        if (SORT_ORDER != null) {
            fileHeader.setSortOrder(SORT_ORDER); // tsato: why would one change the sort order?
        }

        if (FILTER == Filter.includePairedIntervals && fileHeader.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new UnsupportedOperationException("Input must be coordinate sorted to use includePairedIntervals");
        }

        final boolean presorted = inputSortOrder.equals(fileHeader.getSortOrder());
        log.info("Filtering [presorted=" + presorted + "] " + INPUT.getURIString() + " -> OUTPUT=" + // tsato: URI string the right thing?
                OUTPUT.getURIString() + " [sortorder=" + fileHeader.getSortOrder().name() + "]");
        // tsato: getURIString() or getRawInputString()
        // create OUTPUT file // tsato: REFERENCE -> PicardHTSPath? Would have to update all of Picard really...
        final SAMFileWriter outputWriter = new SAMFileWriterFactory().makeWriter(fileHeader, presorted, OUTPUT.toPath(), REFERENCE_SEQUENCE);
        // ^ makeWriter is deprecated with REFERENCE type being File. Will update when Chris's reference PR is merged
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Written");

        while (filteringIterator.hasNext()) {
            final SAMRecord rec = filteringIterator.next();
            outputWriter.addAlignment(rec);
            progress.record(rec);
        }

        filteringIterator.close();
        outputWriter.close();
        log.info(new DecimalFormat("#,###").format(progress.getCount()) + " SAMRecords written to " + OUTPUT.getRawInputString());
    }

    /**
     * Write out a file of read names for debugging purposes.
     *
     * @param samOrBamFile The SAM/BAM/CRAM file for which we are going to write out a file of its
     *                     containing read names
     */
    private void writeReadsFile(final Path samOrBamFile) throws IOException {
        final PicardHtsPath readsFile =new PicardHtsPath(OUTPUT.getURIString() + ".reads"); // tsato: probably not exacdtly right
                // new File(OUTPUT.getParentFile(), IOUtil.basename(samOrBamFile) + ".reads");
        // IOUtil.assertFileIsWritable(readsFile); // tsato: what to do about it, probably can implement a method in HTSJDK
        try (final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(samOrBamFile);
             final BufferedWriter bw = IOUtil.openFileForBufferedWriting(readsFile.toPath(), null)) { // tsato: what should theo OpenOption be?

            for (final SAMRecord rec : reader) {
                bw.write(rec.toString() + "\n");
            }
        }
        IOUtil.assertPathsAreReadable(Collections.singletonList(readsFile.toPath()));
    }

    private List<Interval> getIntervalList(final Path intervalFile) throws IOException {
        IOUtil.assertFileIsReadable(intervalFile);
        return IntervalList.fromPath(intervalFile).getIntervals();
    }

    @Override
    protected int doWork() {

        try {
            IOUtil.assertFileIsReadable(INPUT.toPath());
            // IOUtil.assertFileIsWritable(OUTPUT);
            if (WRITE_READS_FILES) writeReadsFile(INPUT.toPath());

            final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT.toPath());
            final FilteringSamIterator filteringIterator;

            // Used for exclude/include tag filter which expects a List<Object> input so casting here
            // This is done to get around poor constructors of TagFilter that should be addressed in
            // https://github.com/samtools/htsjdk/issues/1082
            List<Object> tagList = (List) TAG_VALUE; // tsato: wtf is this?

            switch (FILTER) {
                case includeAligned: // tsato: does this option assume that the reads are query-name sorted?
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new AlignedFilter(true), true);
                    break;
                case excludeAligned:
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new AlignedFilter(false), true);
                    break;
                case includeReadList:
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new ReadNameFilter(READ_LIST_FILE, true)); // tsato: need to update htsjdk here...
                    break;
                case excludeReadList:
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new ReadNameFilter(READ_LIST_FILE, false));
                    break;
                case includeJavascript:
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new JavascriptSamRecordFilter(
                                    JAVASCRIPT_FILE,
                                    samReader.getFileHeader()));
                    break;
                case includePairedIntervals:
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new IntervalKeepPairFilter(getIntervalList(INTERVAL_LIST.toPath())));
                    break;
                case includeTagValues:
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new TagFilter(TAG, tagList, true));
                    break;
                case excludeTagValues:
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new TagFilter(TAG, tagList, false));
                    break;
                default:
                    throw new UnsupportedOperationException(FILTER.name() + " has not been implemented!");
            }

            filterReads(filteringIterator);

            IOUtil.assertFileIsReadable(OUTPUT.toPath());
            if (WRITE_READS_FILES) writeReadsFile(OUTPUT.toPath());
            return 0;

        } catch (Exception e) {
            log.error(e, "Failed to filter " + INPUT.getRawInputString());

            if (Files.exists(OUTPUT.toPath())) {
                try {
                    Files.delete(OUTPUT.toPath());
                } catch (IOException ex) {
                    log.warn("Failed to delete possibly incomplete output file:" + OUTPUT.getRawInputString());
                    // tsato: good enough to first order.
                }
            }

            return 1;
        }
    }

    @Override
    protected String[] customCommandLineValidation() {

        List<String> errors = new ArrayList<>();

        if (INPUT.equals(OUTPUT)) errors.add("INPUT file and OUTPUT file must differ!");

        List<Filter> tagFilters = Arrays.asList(Filter.includeTagValues, Filter.excludeTagValues);

        checkInputs(Arrays.asList(Filter.includeReadList, Filter.excludeReadList), READ_LIST_FILE, "READ_LIST_FILE").ifPresent(errors::add);
        checkInputs(Collections.singletonList(Filter.includePairedIntervals), INTERVAL_LIST, "INTERVAL_LIST").ifPresent(errors::add);
        checkInputs(Collections.singletonList(Filter.includeJavascript), JAVASCRIPT_FILE, "JAVASCRIPT_FILE").ifPresent(errors::add);
        checkInputs(tagFilters, TAG, "TAG").ifPresent(errors::add);

        if (tagFilters.contains(FILTER) && TAG_VALUE.isEmpty()) {
            log.warn("Running FilterSamReads with a Tag Filter but no TAG_VALUE argument provided.  This " +
                    "will recreate the original input file i.e. not filter anything");
        }

        if (!errors.isEmpty()) return errors.toArray(new String[errors.size()]);

        return super.customCommandLineValidation();
    }

    private Optional<String> checkInputs(final List<Filter> filters, final Object inputObject, final String inputFileVariable) {
        if (filters.contains(FILTER) && inputObject == null)
            return Optional.of(String.format("%s must be specified when using FILTER=%s, but it was null.", inputFileVariable, FILTER));
        if (!filters.contains(FILTER) && inputObject != null)
            return Optional.of(String.format("%s may only be specified when using FILTER from %s, FILTER value: %s, %s value: %s",
                    inputFileVariable, String.join(", ", filters.stream().map(Enum::toString).collect(Collectors.toList())),
                    FILTER, inputFileVariable, inputObject));
        return Optional.empty();
    }
}
