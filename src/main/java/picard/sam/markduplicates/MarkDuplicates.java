/*
 * The MIT License
 *
 * Copyright (c) 2009-2016 The Broad Institute
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

package picard.sam.markduplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.*;
import picard.sam.util.ParallelBAMFileWriter;
import picard.sam.util.ParallelBAMRecordIterator;
import picard.sam.util.RawBlockBAMRecord;
import picard.sam.util.ReadAheadIterator;
import picard.sam.util.ReadNameParser;
import picard.sam.util.RepresentativeReadIndexer;

import java.io.File;
import java.io.IOException;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Objects;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.LockSupport;
import java.util.function.Consumer;
import java.util.function.LongConsumer;
import java.util.function.Supplier;

/**
 * A better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = MarkDuplicates.USAGE_SUMMARY + MarkDuplicates.USAGE_DETAILS,
        oneLineSummary = MarkDuplicates.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class MarkDuplicates extends AbstractMarkDuplicatesCommandLineProgram implements MarkDuplicatesHelper {
    static final String USAGE_SUMMARY = "Identifies duplicate reads.  ";
    static final String USAGE_DETAILS = "<p>This tool locates and tags duplicate reads in a SAM, BAM or CRAM file, where duplicate reads are " +
            "defined as originating from a single fragment of DNA.  Duplicates can arise during sample preparation e.g. library " +
            "construction using PCR.  See also " +
            "<a href='https://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity'>EstimateLibraryComplexity</a>" +
            " for additional notes on PCR duplication artifacts.  Duplicate reads can also result from a single amplification cluster, " +
            "incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.  These duplication artifacts are " +
            "referred to as optical duplicates.</p>" +
            "" +
            "<p>The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file.  " +
            "A BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes.  After duplicate reads are" +
            " collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums " +
            "of their base-quality scores (default method). Note that this is different from directly checking if the sequences match, which " +
            "MarkDuplicates does not do.</p>  " +
            "" +
            "<p>The tool's main output is a new SAM, BAM or CRAM file, in which duplicates have been identified in the SAM flags field for each" +
            " read.  Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024.  " +
            "If you are not familiar with this type of annotation, please see the following " +
            "<a href='https://www.broadinstitute.org/gatk/blog?id=7019'>blog post</a> for additional information.</p>" +
            "" +
            "<p>Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type of " +
            "duplicate.  To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in  " +
            "the 'optional field' section of a SAM/BAM/CRAM file.  Invoking the TAGGING_POLICY option," +
            " you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no " +
            "duplicates (DontTag).  The records within the output of a SAM/BAM/CRAM file will have values for the 'DT' tag (depending on the invoked " +
            "TAGGING_POLICY), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ).  " +
            "This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify " +
            "and differentiate duplicate types.  Set READ_NAME_REGEX to null to skip optical duplicate detection, e.g. for RNA-seq " +
            "or other data where duplicate sets are extremely large and estimating library complexity is not an aim.  " +
            "Note that without optical duplicate counts, library size estimation will be inaccurate.</p> " +

            "<p>MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.</p>  " +

            "<p>The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different.  " +
            "When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not " +
            "marked as duplicates.  However, when the input is query-sorted (actually query-grouped), " +
            "then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be" +
            " marked as duplicate reads.</p>  " +

            "<p>If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.</p>" +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar MarkDuplicates \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=marked_duplicates.bam \\<br />" +
            "      M=marked_dup_metrics.txt" +
            "</pre>" +
            "" +
            "Please see " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics'>MarkDuplicates</a> " +
            "for detailed explanations of the output metrics." +
            "<hr />";

    /**
     * Enum used to control how duplicates are flagged in the DT optional tag on each read.
     */
    public enum DuplicateTaggingPolicy {
        DontTag, OpticalOnly, All
    }

    /**
     * The optional attribute in SAM/BAM/CRAM files used to store the duplicate type.
     */
    public static final String DUPLICATE_TYPE_TAG = "DT";
    /**
     * The duplicate type tag value for duplicate type: library.
     */
    public static final String DUPLICATE_TYPE_LIBRARY = "LB";
    /**
     * The duplicate type tag value for duplicate type: sequencing (optical & pad-hopping, or "co-localized").
     */
    public static final String DUPLICATE_TYPE_SEQUENCING = "SQ";
    /**
     * The attribute in the SAM/BAM file used to store which read was selected as representative out of a duplicate set
     */
    public static final String DUPLICATE_SET_INDEX_TAG = "DI";
    /**
     * The attribute in the SAM/BAM file used to store the size of a duplicate set
     */
    public static final String DUPLICATE_SET_SIZE_TAG = "DS";

    /**
     * Enum for the possible values that a duplicate read can be tagged with in the DT attribute.
     */
    public enum DuplicateType {
        LIBRARY(DUPLICATE_TYPE_LIBRARY),
        SEQUENCING(DUPLICATE_TYPE_SEQUENCING);

        private final String code;

        DuplicateType(final String code) {
            this.code = code;
        }

        public String code() {
            return this.code;
        }
    }

    private final Log log = Log.getInstance(MarkDuplicates.class);

    /**
     * If more than this many sequences in SAM file, don't spill to disk because there will not
     * be enough file handles.
     */
    @Argument(shortName = "MAX_SEQS",
            doc = "This option is obsolete. ReadEnds will always be spilled to disk.")
    public int MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP = 50000;

    @Argument(shortName = "MAX_FILE_HANDLES",
            doc = "Maximum number of file handles to keep open when spilling read ends to disk. " +
                    "Set this number a little lower than the per-process maximum number of file that may be open. " +
                    "This number can be found by executing the 'ulimit -n' command on a Unix system.")
    public int MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 8000;

    @Argument(doc = "This number, plus the maximum RAM available to the JVM, determine the memory footprint used by " +
            "some of the sorting collections.  If you are running out of memory, try reducing this number.")
    public double SORTING_COLLECTION_SIZE_RATIO = 0.25;

    @Argument(doc = "Barcode SAM tag (ex. BC for 10X Genomics)", optional = true)
    public String BARCODE_TAG = null;

    @Argument(doc = "Read one barcode SAM tag (ex. BX for 10X Genomics)", optional = true)
    public String READ_ONE_BARCODE_TAG = null;

    @Argument(doc = "Read two barcode SAM tag (ex. BX for 10X Genomics)", optional = true)
    public String READ_TWO_BARCODE_TAG = null;

    @Argument(doc = "If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG (DS), " +
            "indicates the size of the duplicate set. The smallest possible DS value is 2 which occurs when two " +
            "reads map to the same portion of the reference only one of which is marked as duplicate. The second " +
            "tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier for the duplicate set to which the " +
            "record belongs. This identifier is the index-in-file of the representative read that was selected out " +
            "of the duplicate set.", optional = true)
    public boolean TAG_DUPLICATE_SET_MEMBERS = false;

    @Argument(doc = "If true remove 'optical' duplicates and other duplicates that appear to have arisen from the " +
            "sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false. " +
            "If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored.")
    public boolean REMOVE_SEQUENCING_DUPLICATES = false;

    @Argument(doc = "Determines how duplicate types are recorded in the DT optional attribute.")
    public DuplicateTaggingPolicy TAGGING_POLICY = DuplicateTaggingPolicy.DontTag;

    @Argument(doc = "Clear DT tag from input SAM records. Should be set to false if input SAM doesn't have this tag.  Default true")
    public boolean CLEAR_DT = true;

    @Argument(doc = "Treat UMIs as being duplex stranded.  This option requires that the UMI consist of two equal length " +
            "strings that are separated by a hyphen (e.g. 'ATC-GTC'). Reads are considered duplicates if, in addition to standard " +
            "definition, have identical normalized UMIs.  A UMI from the 'bottom' strand is normalized by swapping its content " +
            "around the hyphen (eg. ATC-GTC becomes GTC-ATC).  A UMI from the 'top' strand is already normalized as it is. " +
            "Both reads from a read pair considered top strand if the read 1 unclipped 5' coordinate is less than the read " +
            "2 unclipped 5' coordinate. All chimeric reads and read fragments are treated as having come from the top strand. " +
            "With this option is it required that the BARCODE_TAG hold non-normalized UMIs. Default false.")
    public boolean DUPLEX_UMI = false;

    @Argument(doc = "SAM tag to uniquely identify the molecule from which a read was derived.  Use of this option requires that " +
            "the BARCODE_TAG option be set to a non null value.  Default null.", optional = true)
    public String MOLECULAR_IDENTIFIER_TAG = null;

    @ArgumentCollection
    public MarkDuplicatesForFlowArgumentCollection flowBasedArguments = new MarkDuplicatesForFlowArgumentCollection();

    protected SortingCollection<ReadEndsForMarkDuplicates> pairSort;
    // Completed pairs as rows rather than objects; null when pairSort holds them.
    private PairedReadEndRows pairedEndRows;
    // The pair pass running beside pass 1; null when it is not run there.
    private PairPassDuringRead pairPassDuringRead;
    // Pass 2 running beside pass 1; null when it is not run there.
    private SecondPassDuringRead secondPassDuringRead;
    protected SortingCollection<ReadEndsForMarkDuplicates> fragSort;
    // Keys of the paired ends left out of fragSort; null when fragSort holds every end (barcodes, flow mode).
    private PairedReadEndKeys pairedEndKeys;
    protected SortingLongCollection duplicateIndexes;
    protected SortingLongCollection opticalDuplicateIndexes;
    protected SortingCollection<RepresentativeReadIndexer> representativeReadIndicesForDuplicates;

    // some calculations are performed using a helper class, which can be parameter specific
    // by default, this instance is the helper
    private MarkDuplicatesHelper calcHelper = this;

    private int numDuplicateIndices = 0;
    // Records pass 1 read; no more duplicate indices than this are ever added.
    private long numRecordsInPass1 = Long.MAX_VALUE;
    // Records pass 1 took from its input, each validated and eagerly decoded there at the default stringency.
    private long numRecordsValidatedInPass1 = 0;
    static private final long NO_SUCH_INDEX = Long.MAX_VALUE; // needs to be large so that that >= test fails for query-sorted traversal
    // Outside the short range of library IDs.
    static private final int NO_LIBRARY_ID = Integer.MIN_VALUE;
    // Records handed from the read-ahead thread to pass 1 at a time.
    private static final int READ_AHEAD_BATCH = 1024;
    // The smallest heap with which pass 2 runs beside pass 1.
    private static final long SECOND_PASS_MIN_HEAP = 1L << 30;

    /**
     * Decoding threads for one BAM input in pass 2; 0 below 4 processors, where pass 2 reads on the main thread.
     */
    private static int decodeThreads() {
        return Runtime.getRuntime().availableProcessors() >= 4 ? 2 : 0;
    }

    /**
     * Decoding threads for one BAM input in pass 1, where they also build read ends and pair-map keys; 0 below 4
     * processors, where pass 1 reads through the read-ahead thread.
     */
    private static int pass1DecodeThreads() {
        final int processors = Runtime.getRuntime().availableProcessors();
        return processors >= 4 ? Math.min(4, processors - 2) : 0;
    }

    protected LibraryIdGenerator libraryIdGenerator = null; // this is initialized in buildSortedReadEndLists

    private int getReadOneBarcodeValue(final SAMRecord record) {
        return EstimateLibraryComplexity.getReadBarcodeValue(record, READ_ONE_BARCODE_TAG);
    }

    private int getReadTwoBarcodeValue(final SAMRecord record) {
        return EstimateLibraryComplexity.getReadBarcodeValue(record, READ_TWO_BARCODE_TAG);
    }

    public MarkDuplicates() {
        DUPLICATE_SCORING_STRATEGY = ScoringStrategy.SUM_OF_BASE_QUALITIES;
    }

    /**
     * Main work method.  Reads the SAM file once and collects sorted information about
     * the 5' ends of both ends of each read (or just one end in the case of pairs).
     * Then makes a pass through those determining duplicates before re-reading the
     * input file and writing it out with duplication flags set correctly.
     */
    protected int doWork() {
        IOUtil.assertInputsAreValid(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(METRICS_FILE);

        final boolean useBarcodes = (null != BARCODE_TAG || null != READ_ONE_BARCODE_TAG || null != READ_TWO_BARCODE_TAG);

        // use flow based calculation helper?
        if ( flowBasedArguments.FLOW_MODE ) {
            calcHelper = new MarkDuplicatesForFlowHelper(this);
        }

        reportMemoryStats("Start of doWork");
        log.info("Reading input file and constructing read end information.");
        try {
            buildSortedReadEndLists(useBarcodes);
        } catch (final RuntimeException | Error e) {
            if (this.pairPassDuringRead != null) {
                this.pairPassDuringRead.abandon();
                this.pairPassDuringRead = null;
            }
            if (this.secondPassDuringRead != null) {
                if (this.secondPassDuringRead.isStarted()) {
                    log.debug("Pass 2 beside pass 1 stopped: pass 1 failed.");
                }
                this.secondPassDuringRead.discard();
                this.secondPassDuringRead = null;
            }
            throw e;
        }
        reportMemoryStats("After buildSortedReadEndLists");
        final SecondPassDuringRead secondPass = this.secondPassDuringRead;
        this.secondPassDuringRead = null;
        if (secondPass != null && writtenDuringRead(secondPass)) {
            log.info("Closed outputs. Getting more Memory Stats.");
            reportMemoryStats("After output close");
            finalizeAndWriteMetrics(libraryIdGenerator, getMetricsFile(), METRICS_FILE);
            return 0;
        }
        calcHelper.generateDuplicateIndexes(useBarcodes, this.REMOVE_SEQUENCING_DUPLICATES || this.TAGGING_POLICY != DuplicateTaggingPolicy.DontTag);

        reportMemoryStats("After generateDuplicateIndexes");
        log.info("Marking " + this.numDuplicateIndices + " records as duplicates.");

        if (this.READ_NAME_REGEX == null) {
            log.warn("Skipped optical duplicate cluster discovery; library size estimation may be inaccurate!");
        } else {
            log.info("Found " + (this.libraryIdGenerator.getNumberOfOpticalDuplicateClusters()) + " optical duplicate clusters.");
        }

        // RawBlockBAMRecords let ParallelBAMFileWriter copy each record's name, cigar, bases and qualities. The records
        // pass 1 took are not validated again.
        final SamHeaderAndIterator headerAndIterator = openInputs(false, false, decodeThreads(), RawBlockBAMRecord.FACTORY,
                null, this.numRecordsValidatedInPass1);
        final SAMFileHeader header = headerAndIterator.header;
        final SAMFileHeader.SortOrder sortOrder = header.getSortOrder();

        final SAMFileHeader outputHeader = header.clone();

        log.info("Reads are assumed to be ordered by: " + sortOrder);

        // Setting the ASSUME_SORT_ORDER to equal queryname is understood to mean that the input is
        // queryname **grouped**. So that's what we set the output order to be, so that the validation will pass
        if (ASSUME_SORT_ORDER == SAMFileHeader.SortOrder.queryname) {
            outputHeader.setGroupOrder(SAMFileHeader.GroupOrder.query);
            outputHeader.setSortOrder(SAMFileHeader.SortOrder.unknown);
            log.info("Output will not be re-sorted. Output header will state SO:unknown GO:query");
        }

        if (ASSUME_SORT_ORDER == null && sortOrder != SAMFileHeader.SortOrder.coordinate && sortOrder != SAMFileHeader.SortOrder.queryname ||
                ASSUME_SORT_ORDER != null && ASSUME_SORT_ORDER != SAMFileHeader.SortOrder.coordinate && ASSUME_SORT_ORDER != SAMFileHeader.SortOrder.queryname) {
            throw new PicardException("This program requires input that are either coordinate or query sorted (according to the header, or at least ASSUME_SORT_ORDER and the content.) " +
                    "Found ASSUME_SORT_ORDER=" + ASSUME_SORT_ORDER + " and header sortorder=" + sortOrder);
        }

        COMMENT.forEach(outputHeader::addComment);

        // Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
        final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);

        // Records are encoded and compressed on other threads, so none may be changed after addAlignment.
        try (SAMFileWriter out = ParallelBAMFileWriter.canWrite(OUTPUT)
                ? ParallelBAMFileWriter.makeAsyncWriter(outputHeader, true, OUTPUT)
                : new SAMFileWriterFactory().setUseAsyncIo(true).makeWriter(outputHeader, true, OUTPUT,REFERENCE_SEQUENCE)){

            // Now copy over the file while marking all the necessary indexes as duplicates
            long recordInFileIndex = 0;
            long nextOpticalDuplicateIndex = this.opticalDuplicateIndexes != null && this.opticalDuplicateIndexes.hasNext() ? this.opticalDuplicateIndexes.next() : NO_SUCH_INDEX;
            long nextDuplicateIndex = (this.duplicateIndexes.hasNext() ? this.duplicateIndexes.next() : NO_SUCH_INDEX);

            // initialize variables for optional representative read tagging
            CloseableIterator<RepresentativeReadIndexer> representativeReadIterator = null;
            RepresentativeReadIndexer rri = null;
            int representativeReadIndexInFile = -1;
            int duplicateSetSize = -1;
            int nextReadInDuplicateSetIndex = -1;
            if (TAG_DUPLICATE_SET_MEMBERS) {
                representativeReadIterator = this.representativeReadIndicesForDuplicates.iterator();
                if (representativeReadIterator.hasNext()) {
                    rri = representativeReadIterator.next();
                    nextReadInDuplicateSetIndex = rri.readIndexInFile;
                    representativeReadIndexInFile = rri.representativeReadIndexInFile;
                    duplicateSetSize = rri.setSize;
                }
            }

            final ProgressLoggerInterface progress = new SingleThreadProgressLogger(log, (int) 1e7, "Written");
            final CloseableIterator<SAMRecord> iterator = headerAndIterator.iterator;
            String duplicateQueryName = null;
            String representativeQueryName = null;

            while (iterator.hasNext()) {
                final SAMRecord rec = iterator.next();

                DuplicationMetrics metrics = AbstractMarkDuplicatesCommandLineProgram.addReadToLibraryMetrics(rec, header, libraryIdGenerator, flowBasedArguments.FLOW_MODE);

                // Now try and figure out the next duplicate index (if going by coordinate. if going by query name, only do this
                // if the query name has changed.
                nextDuplicateIndex = nextIndexIfNeeded(sortOrder, recordInFileIndex, nextDuplicateIndex, duplicateQueryName, rec, this.duplicateIndexes);

                final boolean isDuplicate = recordInFileIndex == nextDuplicateIndex ||
                        (sortOrder == SAMFileHeader.SortOrder.queryname &&
                                recordInFileIndex > nextDuplicateIndex && rec.getReadName().equals(duplicateQueryName));

                if (isDuplicate) {
                    rec.setDuplicateReadFlag(true);

                    metrics.addDuplicateReadToMetrics(rec);
                } else {
                    rec.setDuplicateReadFlag(false);
                }
                nextOpticalDuplicateIndex = nextIndexIfNeeded(sortOrder, recordInFileIndex, nextOpticalDuplicateIndex, duplicateQueryName, rec, this.opticalDuplicateIndexes);

                final boolean isOpticalDuplicate = sortOrder == SAMFileHeader.SortOrder.queryname &&
                        recordInFileIndex > nextOpticalDuplicateIndex &&
                        rec.getReadName().equals(duplicateQueryName) ||
                        recordInFileIndex == nextOpticalDuplicateIndex;

                if (CLEAR_DT) {
                    rec.setAttribute(DUPLICATE_TYPE_TAG, null);
                }

                if (this.TAGGING_POLICY != DuplicateTaggingPolicy.DontTag && rec.getDuplicateReadFlag()) {
                    if (isOpticalDuplicate) {
                        rec.setAttribute(DUPLICATE_TYPE_TAG, DuplicateType.SEQUENCING.code());
                    } else if (this.TAGGING_POLICY == DuplicateTaggingPolicy.All) {
                        rec.setAttribute(DUPLICATE_TYPE_TAG, DuplicateType.LIBRARY.code());
                    }
                }

                // Tag any read pair that was in a duplicate set with the duplicate set size and a representative read name
                if (TAG_DUPLICATE_SET_MEMBERS) {
                    final boolean needNextRepresentativeIndex = recordInFileIndex > nextReadInDuplicateSetIndex &&
                            (sortOrder == SAMFileHeader.SortOrder.coordinate || !rec.getReadName().equals(representativeQueryName));
                    if (needNextRepresentativeIndex && representativeReadIterator.hasNext()) {
                        rri = representativeReadIterator.next();
                        nextReadInDuplicateSetIndex = rri.readIndexInFile;
                        representativeReadIndexInFile = rri.representativeReadIndexInFile;
                        duplicateSetSize = rri.setSize;
                    }

                    /* If this record's index is readInDuplicateSetIndex, then it is in a duplicateset.
                    For queryname sorted data, we only have one representativeReadIndex entry per read name, so we need
                    to also look for additional reads with the same name.
                     */
                    final boolean isInDuplicateSet = recordInFileIndex == nextReadInDuplicateSetIndex ||
                            (sortOrder == SAMFileHeader.SortOrder.queryname &&
                                    recordInFileIndex > nextReadInDuplicateSetIndex && rec.getReadName().equals(representativeQueryName));
                    if (isInDuplicateSet) {
                        if (!rec.isSecondaryOrSupplementary() && !rec.getReadUnmappedFlag()) {
                            rec.setAttribute(DUPLICATE_SET_INDEX_TAG, representativeReadIndexInFile);
                            rec.setAttribute(DUPLICATE_SET_SIZE_TAG, duplicateSetSize);
                            representativeQueryName = rec.getReadName();
                        }
                    }
                }

                // Set MOLECULAR_IDENTIFIER_TAG for SAMRecord rec
                if (BARCODE_TAG != null) {
                    UmiUtil.setMolecularIdentifier(rec, "", MOLECULAR_IDENTIFIER_TAG, DUPLEX_UMI);
                }


                // Note, duplicateQueryName must be incremented after we have marked both optical and sequencing duplicates for queryname sorted files.
                if (isDuplicate) {
                    duplicateQueryName = rec.getReadName();
                }

                // Output the record if desired and bump the record index
                recordInFileIndex++;
                if (this.REMOVE_DUPLICATES && rec.getDuplicateReadFlag()) {
                    continue;
                }
                if (this.REMOVE_SEQUENCING_DUPLICATES && isOpticalDuplicate) {
                    continue;
                }
                if (PROGRAM_RECORD_ID != null && pgTagArgumentCollection.ADD_PG_TAG_TO_READS) {
                    rec.setAttribute(SAMTag.PG.name(), chainedPgIds.get(rec.getStringAttribute(SAMTag.PG.name())));
                }
                out.addAlignment(rec);
                progress.record(rec);
            }

            // remember to close the inputs
            log.info("Writing complete. Closing input iterator.");
            iterator.close();

            log.info("Duplicate Index cleanup.");
            this.duplicateIndexes.cleanup();
            if (TAG_DUPLICATE_SET_MEMBERS) {
                log.info("Representative read Index cleanup.");
                this.representativeReadIndicesForDuplicates.cleanup();
            }

            log.info("Getting Memory Stats.");
            reportMemoryStats("Before output close");
        }
        log.info("Closed outputs. Getting more Memory Stats.");
        reportMemoryStats("After output close");

        // Write out the metrics
        finalizeAndWriteMetrics(libraryIdGenerator, getMetricsFile(), METRICS_FILE);

        return 0;
    }

    /**
     * package-visible for testing
     */
    long numOpticalDuplicates() {
        return ((long) this.libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap().getSumOfValues());
    } // cast as long due to returning a double

    /**
     * Print out some quick JVM memory stats.
     */
    private void reportMemoryStats(final String stage) {
        System.gc();
        final Runtime runtime = Runtime.getRuntime();
        log.info(stage + " freeMemory: " + runtime.freeMemory() + "; totalMemory: " + runtime.totalMemory() +
                "; maxMemory: " + runtime.maxMemory());
    }

    /**
     * After pass 1, whether pass 2 run beside it wrote the output: the pair pass and pass 2 then both finished
     * without being abandoned, and the output was moved into place. The duplicate counts, histograms and library
     * metrics are then those generateDuplicateIndexes and pass 2 would have made, and the log has the lines they would
     * have logged, in their order: pass 2 holds its own back until the lines before it are logged, and reads the
     * records pass 1 left undecided and the unplaced reads only after that. Otherwise pass 2 is stopped and its output
     * deleted, and nothing else is changed, with one exception: when pass 2 fails on a record after its lines were
     * released, as pass 2 after pass 1 would on the same record, its output is moved into place as that one's would be
     * left, and the failure is thrown.
     */
    private boolean writtenDuringRead(final SecondPassDuringRead secondPass) {
        final PairPassDuringRead pairPass = this.pairPassDuringRead;
        if (pairPass == null || secondPass.isAbandoned() || !pairPass.finish() || secondPass.isAbandoned()) {
            if (secondPass.isStarted()) {
                log.debug("Pass 2 beside pass 1 abandoned during pass 1; pass 2 runs after it.");
            }
            secondPass.discard();
            return false;
        }
        final boolean indexOpticalDuplicates = this.REMOVE_SEQUENCING_DUPLICATES || this.TAGGING_POLICY != DuplicateTaggingPolicy.DontTag;
        final long duplicates = pairPass.numDuplicates + secondPass.fragmentDuplicates;
        // As generateDuplicateIndexes, doWork and the start of pass 2 log them.
        log.info("Will retain up to " + duplicateIndicesInMemory(indexOpticalDuplicates) + " duplicate indices before spilling to disk.");
        log.info("Traversing read pair information and detecting duplicates.");
        log.info("Traversing fragment information and detecting duplicates.");
        log.info("Sorting list of duplicate records.");
        reportMemoryStats("After generateDuplicateIndexes");
        log.info("Marking " + duplicates + " records as duplicates.");
        if (this.READ_NAME_REGEX == null) {
            log.warn("Skipped optical duplicate cluster discovery; library size estimation may be inaccurate!");
        } else {
            log.info("Found " + pairPass.numberOfOpticalDuplicateClusters() + " optical duplicate clusters.");
        }
        log.info("Reads are assumed to be ordered by: " + secondPass.sortOrder());
        secondPass.releaseLog();

        final boolean written = secondPass.awaitWritten();
        final RuntimeException failure = secondPass.failureAfterRelease();
        if ((written || failure != null) && secondPass.moveTo(OUTPUT)) {
            if (failure != null) {
                log.debug("Pass 2 beside pass 1 failed after pass 1; its output is left as pass 2's would be.");
                throw failure;
            }
            log.debug("Pass 2 beside pass 1 wrote the output.");
            this.pairPassDuringRead = null;
            this.pairedEndRows = null;
            pairPass.addHistograms();
            secondPass.addMetricsTo(this.libraryIdGenerator);
            this.numDuplicateIndices += (int) duplicates;
            this.fragSort.cleanup();
            this.fragSort = null;
            this.pairedEndKeys = null;
            return true;
        }
        log.debug("Pass 2 beside pass 1 abandoned after pass 1; pass 2 runs again.");
        secondPass.discard();
        return false;
    }

    /**
     * Goes through all the records in a file and generates a set of ReadEndsForMarkDuplicates objects that
     * hold the necessary information (reference sequence, 5' read coordinate) to do
     * duplication, caching to disk as necessary to sort them.
     */
    private void buildSortedReadEndLists(final boolean useBarcodes) {
        final int sizeInBytes;
        if (useBarcodes) {
            sizeInBytes = ReadEndsForMarkDuplicatesWithBarcodes.getSizeOf();
        } else {
            sizeInBytes = ReadEndsForMarkDuplicates.getSizeOf();
        }
        MAX_RECORDS_IN_RAM = Math.min((int) (Runtime.getRuntime().maxMemory() / sizeInBytes) / 2, Integer.MAX_VALUE - 32);
        final int maxInMemory = Math.min((int) ((Runtime.getRuntime().maxMemory() * SORTING_COLLECTION_SIZE_RATIO) / sizeInBytes), Integer.MAX_VALUE - 32);
        log.info("Will retain up to " + maxInMemory + " data points before spilling to disk.");

        final ReadEndsForMarkDuplicatesCodec fragCodec, pairCodec, diskCodec;
        if (useBarcodes) {
            fragCodec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
            pairCodec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
            diskCodec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        } else {
            fragCodec = new ReadEndsForMarkDuplicatesCodec();
            pairCodec = new ReadEndsForMarkDuplicatesCodec();
            diskCodec = new ReadEndsForMarkDuplicatesCodec();
        }

        //In Illumina the two ends of the fragment are saved only if the read is paired
        // in flow based sequencing we are interested in both ends of the read
        if (flowBasedArguments.FLOW_USE_END_IN_UNPAIRED_READS){
            fragCodec.setForceSerialization();
            pairCodec.setForceSerialization();
            diskCodec.setForceSerialization();
        }

        // Not when a subclass, barcodes or flow mode may read pairSort. pairSort is made when the rows refuse an end.
        final boolean pairsAsRows = !useBarcodes && calcHelper == this && getClass() == MarkDuplicates.class;
        final Supplier<SortingCollection<ReadEndsForMarkDuplicates>> newPairSort = () ->
                SortingCollection.newInstance(ReadEndsForMarkDuplicates.class,
                        pairCodec,
                        new ReadEndsMDComparator(useBarcodes),
                        maxInMemory,
                        TMP_DIR);
        this.pairedEndRows = pairsAsRows ? new PairedReadEndRows(maxInMemory) : null;
        this.pairSort = pairsAsRows ? null : newPairSort.get();
        this.fragSort = SortingCollection.newInstance(ReadEndsForMarkDuplicates.class,
                fragCodec,
                new ReadEndsMDComparator(useBarcodes),
                maxInMemory,
                TMP_DIR);

        // The decoding threads build each read's ends as buildReadEnds does, but for its library and index in the file.
        // Not with a subclass, which may override buildReadEnds, nor with barcodes or flow mode.
        final boolean endsWhileDecoding = !useBarcodes && calcHelper == this && getClass() == MarkDuplicates.class &&
                opticalDuplicateFinder != null && opticalDuplicateFinder.getClass() == OpticalDuplicateFinder.class;
        final String readNameRegex = endsWhileDecoding ? opticalDuplicateFinder.getReadNameRegex() : null;
        final boolean collectPgIds = PROGRAM_RECORD_ID != null;
        // Not htsjdk's async BGZF read, even for the header: in htsjdk 5.0.0 it hangs on a corrupt block.
        final SamHeaderAndIterator headerAndIterator = openInputs(true, false, pass1DecodeThreads(), null,
                !endsWhileDecoding ? null : () -> {
                    final ReadNameParser parser = new ReadNameParser(readNameRegex);
                    final ReadGroupSlots readGroupSlots = new ReadGroupSlots();
                    return rec -> decodeRead(rec, parser, collectPgIds, readGroupSlots);
                });
        final SAMFileHeader.SortOrder assumedSortOrder = headerAndIterator.header.getSortOrder();
        final SAMFileHeader header = headerAndIterator.header;
        this.pairedEndKeys = useBarcodes || flowBasedArguments.FLOW_MODE ? null :
                new PairedReadEndKeys(header.getSequenceDictionary().size(),
                        (int) Math.min(Runtime.getRuntime().maxMemory() * SORTING_COLLECTION_SIZE_RATIO / Long.BYTES,
                                Integer.MAX_VALUE - 32),
                        this.fragSort::add);
        final ReadEndsForMarkDuplicatesMap tmp = assumedSortOrder == SAMFileHeader.SortOrder.queryname ?
                new MemoryBasedReadEndsForMarkDuplicatesMap() :  new DiskBasedReadEndsForMarkDuplicatesMap(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, diskCodec);
        long index = 0;
        final ProgressLoggerInterface progress = new SingleThreadProgressLogger(log, (int) 1e6, "Read");
        // Records are decoded and validated on other threads while this one builds read ends.
        final CloseableIterator<SAMRecord> iterator = headerAndIterator.iterator instanceof ParallelBAMRecordIterator ?
                headerAndIterator.iterator :
                new ReadAheadIterator<>(headerAndIterator.iterator, READ_AHEAD_BATCH, "MarkDuplicates read-ahead");
        final ParallelBAMRecordIterator decodingIterator = endsWhileDecoding &&
                iterator instanceof ParallelBAMRecordIterator ? (ParallelBAMRecordIterator) iterator : null;
        // Set once the optical duplicate finder has failed to parse a read name: it warns only the first time.
        boolean readNameUnparsed = false;

        if (null == this.libraryIdGenerator) {
            this.libraryIdGenerator = new LibraryIdGenerator(header, flowBasedArguments.FLOW_MODE);
        }
        // Library IDs by DecodedRead.readGroupSlot, each set by getLibraryId on its read group's first read.
        final int[] libraryIdsBySlot = new int[header.getReadGroups().size()];
        Arrays.fill(libraryIdsBySlot, NO_LIBRARY_ID);

        // The pair pass beside pass 1: for one coordinate-sorted input whose pairs go to the rows, without duplicate
        // set members, with no optical duplicate finder but OpticalDuplicateFinder, which keeps no state between sets.
        PairPassDuringRead pairPass = null;
        // Pass 2 beside them, writing a new BAM file next to OUTPUT, on an input read by ParallelBAMRecordIterator, with
        // a heap of 1 GiB or more: both passes' buffers together can exhaust a small heap that pass 1 alone fits in.
        SecondPassDuringRead secondPass = null;
        if (pairsAsRows && INPUT.size() == 1 && assumedSortOrder == SAMFileHeader.SortOrder.coordinate &&
                !TAG_DUPLICATE_SET_MEMBERS && Runtime.getRuntime().availableProcessors() >= 4 &&
                (READ_NAME_REGEX == null || opticalDuplicateFinder != null &&
                        opticalDuplicateFinder.getClass() == OpticalDuplicateFinder.class)) {
            if (iterator instanceof ParallelBAMRecordIterator && decodeThreads() > 0 &&
                    Runtime.getRuntime().maxMemory() >= SECOND_PASS_MIN_HEAP && ParallelBAMFileWriter.canWrite(OUTPUT) &&
                    Files.notExists(OUTPUT.toPath(), LinkOption.NOFOLLOW_LINKS)) {
                secondPass = new SecondPassDuringRead(header.getSequenceDictionary().size(),
                        this.REMOVE_SEQUENCING_DUPLICATES || this.TAGGING_POLICY != DuplicateTaggingPolicy.DontTag,
                        Math.max(1 << 16, maxInMemory / 8));
                this.secondPassDuringRead = secondPass;
            }
            pairPass = new PairPassDuringRead(maxInMemory, secondPass);
            this.pairPassDuringRead = pairPass;
        }

        String duplicateQueryName = null;
        long duplicateIndex = NO_SUCH_INDEX;
        boolean stoppedAtUnplacedReads = false;
        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            // Null when the record was read on this thread or a lookup on the decoding thread threw.
            final DecodedRead decoded = decodingIterator == null ? null : (DecodedRead) decodingIterator.annotation();

            // This doesn't have anything to do with building sorted ReadEnd lists, but it can be done in the same pass
            // over the input
            if (PROGRAM_RECORD_ID != null) {
                // Gather all PG IDs seen in merged input files in first pass.  These are gathered for two reasons:
                // - to know how many different PG records to create to represent this program invocation.
                // - to know what PG IDs are already used to avoid collisions when creating new ones.
                // Note that if there are one or more records that do not have a PG tag, then a null value
                // will be stored in this set.
                pgIdsSeen.add(decoded != null ? decoded.pg : rec.getStringAttribute(SAMTag.PG.name()));
                // Pass 2 gave its output the header these IDs made when it started.
                if (secondPass != null && secondPass.pgIdsSeenChanged(pgIdsSeen.size())) {
                    secondPass.abandon();
                }
            }

            // If working in query-sorted, need to keep index of first record with any given query-name.
            if (assumedSortOrder == SAMFileHeader.SortOrder.queryname && !rec.getReadName().equals(duplicateQueryName)) {
                duplicateQueryName = rec.getReadName();
                duplicateIndex = index;
            }

            if (rec.getReadUnmappedFlag()) {
                if (rec.getReferenceIndex() == -1 && assumedSortOrder == SAMFileHeader.SortOrder.coordinate) {
                    // When we hit the unmapped reads with no coordinate, no reason to continue (only in coordinate sort).
                    stoppedAtUnplacedReads = true;
                    break;
                }
                // If this read is unmapped but sorted with the mapped reads, just skip it.

            } else if (!rec.isSecondaryOrSupplementary()) {
                final long indexForRead = assumedSortOrder == SAMFileHeader.SortOrder.queryname ? duplicateIndex : index;
                final ReadEndsForMarkDuplicates fragmentEnd;
                if (decoded == null) {
                    fragmentEnd = calcHelper.buildReadEnds(header, indexForRead, rec, useBarcodes);
                } else {
                    fragmentEnd = decoded.ends;
                    fragmentEnd.read1IndexInFile = indexForRead;
                    final int slot = decoded.readGroupSlot;
                    if (slot < 0) {
                        fragmentEnd.libraryId = libraryIdGenerator.getLibraryId(rec);
                    } else {
                        if (libraryIdsBySlot[slot] == NO_LIBRARY_ID) {
                            libraryIdsBySlot[slot] = libraryIdGenerator.getLibraryId(rec);
                        }
                        fragmentEnd.libraryId = (short) libraryIdsBySlot[slot];
                    }
                    // No location: parse the name here too, so that the finder prints its warning, in file order.
                    if (fragmentEnd.readGroup == -1 && !readNameUnparsed) {
                        if (this.opticalDuplicateFinder.addLocationInformation(rec.getReadName(), fragmentEnd)) {
                            fragmentEnd.readGroup = readGroupIndex(header, rec);
                        } else {
                            readNameUnparsed = true;
                        }
                    }
                }
                final boolean keptAsKey = this.pairedEndKeys != null && fragmentEnd.isPaired() &&
                        this.pairedEndKeys.add(fragmentEnd);
                if (!keptAsKey) {
                    this.fragSort.add(fragmentEnd);
                }

                final boolean withMappedMate = rec.getReadPairedFlag() && !rec.getMateUnmappedFlag();
                if (secondPass != null && pairPass != null && secondPass.isActive()) {
                    if (withMappedMate) {
                        pairPass.addPairedKey(fragmentEnd);
                    } else {
                        pairPass.addFragment(fragmentEnd);
                    }
                }

                if (withMappedMate) {
                    if (pairPass != null && !pairPass.atOrAboveFrontier(fragmentEnd)) {
                        pairPass.abandon();
                        pairPass = null;
                    }
                    final String key;
                    if (decoded != null && decoded.pairKey != null) {
                        key = decoded.pairKey;
                    } else {
                        final StringBuilder keyBuilder = new StringBuilder();
                        keyBuilder.append(rec.getReadGroup().getReadGroupId());
                        keyBuilder.append(rec.getReadName());
                        key = keyBuilder.toString();
                    }
                    ReadEndsForMarkDuplicates pairedEnds = tmp.remove(rec.getReferenceIndex(), key);

                    // See if we've already seen the first end or not
                    if (pairedEnds == null) {
                        // pairedEnds is modified when the mate comes along, so it is a copy when fragSort holds fragmentEnd.
                        pairedEnds = keptAsKey ? fragmentEnd : fragmentEnd.clone();
                        tmp.put(pairedEnds.read2ReferenceIndex, key, pairedEnds);
                        if (secondPass != null && secondPass.isActive()) {
                            secondPass.firstMate(indexForRead, rec);
                        }
                    } else {
                        if (secondPass != null && secondPass.isActive()) {
                            secondPass.completed(pairedEnds.read1IndexInFile);
                        }
                        final int matesRefIndex = fragmentEnd.read1ReferenceIndex;
                        final int matesCoordinate = fragmentEnd.read1Coordinate;

                        // Set orientationForOpticalDuplicates, which always goes by the first then the second end for the strands.  NB: must do this
                        // before updating the orientation later.
                        if (rec.getFirstOfPairFlag()) {
                            pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(rec.getReadNegativeStrandFlag(), pairedEnds.orientation == ReadEnds.R);
                            if (useBarcodes) {
                                ((ReadEndsForMarkDuplicatesWithBarcodes) pairedEnds).readOneBarcode = getReadOneBarcodeValue(rec);
                            }
                        } else {
                            pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R, rec.getReadNegativeStrandFlag());
                            if (useBarcodes) {
                                ((ReadEndsForMarkDuplicatesWithBarcodes) pairedEnds).readTwoBarcode = getReadTwoBarcodeValue(rec);
                            }
                        }

                        // If the other read is actually later, simply add the other read's data as read2, else flip the reads
                        if (matesRefIndex > pairedEnds.read1ReferenceIndex ||
                                (matesRefIndex == pairedEnds.read1ReferenceIndex && matesCoordinate >= pairedEnds.read1Coordinate)) {
                            pairedEnds.read2ReferenceIndex = matesRefIndex;
                            pairedEnds.read2Coordinate = matesCoordinate;
                            pairedEnds.read2IndexInFile = indexForRead;
                            pairedEnds.orientation = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R,
                                    rec.getReadNegativeStrandFlag());

                            // if the two read ends are in the same position, pointing in opposite directions,
                            // the orientation is undefined and the procedure above
                            // will depend on the order of the reads in the file.
                            // To avoid this, we set it explicitly (to FR):
                            if (pairedEnds.read2ReferenceIndex == pairedEnds.read1ReferenceIndex &&
                                    pairedEnds.read2Coordinate == pairedEnds.read1Coordinate &&
                                    pairedEnds.orientation == ReadEnds.RF) {
                                pairedEnds.orientation = ReadEnds.FR;
                            }
                        } else {
                            pairedEnds.read2ReferenceIndex = pairedEnds.read1ReferenceIndex;
                            pairedEnds.read2Coordinate = pairedEnds.read1Coordinate;
                            pairedEnds.read2IndexInFile = pairedEnds.read1IndexInFile;
                            pairedEnds.read1ReferenceIndex = matesRefIndex;
                            pairedEnds.read1Coordinate = matesCoordinate;
                            pairedEnds.read1IndexInFile = indexForRead;
                            pairedEnds.orientation = ReadEnds.getOrientationByte(rec.getReadNegativeStrandFlag(),
                                    pairedEnds.orientation == ReadEnds.R);
                        }

                        // A read's ends from a decoding thread hold its score, which getReadDuplicateScore would return.
                        pairedEnds.score += decoded != null ? fragmentEnd.score :
                                calcHelper.getReadDuplicateScore(rec, pairedEnds);
                        if (this.pairedEndRows != null && !this.pairedEndRows.add(pairedEnds)) {
                            log.debug("Pair rows refused a pair after ", this.pairedEndRows.size(),
                                    " rows; completed pairs go to pairSort.");
                            this.pairSort = newPairSort.get();
                            this.pairedEndRows.drainTo(this.pairSort::add);
                            this.pairedEndRows = null;
                            if (pairPass != null) {
                                pairPass.abandon();
                                pairPass = null;
                            }
                        }
                        if (this.pairedEndRows == null) {
                            this.pairSort.add(pairedEnds);
                        } else if (pairPass != null) {
                            pairPass.add(pairedEnds);
                        }
                    }
                }
            }

            // Print out some stats every 1m reads
            ++index;
            if (pairPass != null && pairPass.readyToSend(index)) {
                if (!pairPass.send(rec, index)) {
                    pairPass = null;
                } else if (secondPass != null && !secondPass.isStarted() && secondPass.isActive()) {
                    startSecondPass(secondPass, index);
                }
            }
            if (progress.record(rec)) {
                log.info("Tracking " + tmp.size() + " as yet unmatched pairs. " + tmp.sizeInRam() + " records in RAM.");
            }
        }

        log.info("Read " + index + " records. " + tmp.size() + " pairs never matched.");
        this.numRecordsInPass1 = index;
        this.numRecordsValidatedInPass1 = stoppedAtUnplacedReads ? index + 1 : index;
        if (secondPass != null) {
            secondPass.passOneDone();
            // Unstarted, it would only now begin; pass 1 ended before its first batch.
            if (!secondPass.isStarted() || pairPass == null) {
                secondPass.abandon();
            }
            secondPass.validated.finish(this.numRecordsValidatedInPass1);
        }
        if (pairPass != null) {
            pairPass.doneAdding();
        }
        iterator.close();

        // Tell these collections to free up memory if possible.
        if (this.pairedEndKeys != null) {
            this.pairedEndKeys.doneAdding();
            if (this.pairedEndKeys.keysHandedOver() > 0) {
                log.debug("Paired-end keys handed to the fragment sort: ", this.pairedEndKeys.keysHandedOver());
            }
        }
        if (this.pairSort != null) {
            this.pairSort.doneAdding();
        }
        this.fragSort.doneAdding();
    }

    /**
     * Calculates score for the duplicate read
     * @param rec - read
     * @param pairedEnds - location of the read ends
     * @return - read score calculated according to the DUPLICATE_SCORING_STRATEGY:
     *         SUM_OF_BASE_QUALITIES, (default)
     *         TOTAL_MAPPED_REFERENCE_LENGTH,
     *         RANDOM
     *
     */
    public short getReadDuplicateScore(final SAMRecord rec, final ReadEndsForMarkDuplicates pairedEnds) {
        return DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);
    }

    /**
     * Builds a read ends object that represents a single read.
     */
    public ReadEndsForMarkDuplicates buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec, final boolean useBarcodes) {
        final ReadEndsForMarkDuplicates ends;

        if (useBarcodes) {
            ends = new ReadEndsForMarkDuplicatesWithBarcodes();
        } else {
            ends = new ReadEndsForMarkDuplicates();
        }
        ends.read1ReferenceIndex = rec.getReferenceIndex();
        ends.read1Coordinate = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();
        ends.orientation = rec.getReadNegativeStrandFlag() ? ReadEnds.R : ReadEnds.F;
        ends.read1IndexInFile = index;
        ends.score = DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);

        // Doing this lets the ends object know that it's part of a pair
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            ends.read2ReferenceIndex = rec.getMateReferenceIndex();
        }

        // Fill in the library ID
        ends.libraryId = libraryIdGenerator.getLibraryId(rec);

        // Fill in the location information for optical duplicates
        if (this.opticalDuplicateFinder.addLocationInformation(rec.getReadName(), ends)) {
            ends.readGroup = readGroupIndex(header, rec);
        }

        if (useBarcodes) {
            final ReadEndsForMarkDuplicatesWithBarcodes endsWithBarcode = (ReadEndsForMarkDuplicatesWithBarcodes) ends;
            final String topStrandNormalizedUmi = UmiUtil.getTopStrandNormalizedUmi(rec, BARCODE_TAG, DUPLEX_UMI);
            endsWithBarcode.barcode = Objects.hash(topStrandNormalizedUmi);

            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                endsWithBarcode.readOneBarcode = getReadOneBarcodeValue(rec);
            } else {
                endsWithBarcode.readTwoBarcode = getReadTwoBarcodeValue(rec);
            }
        }

        return ends;
    }

    /** The RG number of a read: its read group's place in the header's list, the list's length when absent. */
    private static short readGroupIndex(final SAMFileHeader header, final SAMRecord rec) {
        short readGroupIndex = 0;
        final String rg = (String) rec.getAttribute(ReservedTagConstants.READ_GROUP_ID);
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();

        if (rg != null && readGroups != null) {
            for (final SAMReadGroupRecord readGroup : readGroups) {
                if (readGroup.getReadGroupId().equals(rg)) {
                    break;
                } else {
                    readGroupIndex++;
                }
            }
        }
        return readGroupIndex;
    }

    /**
     * htsjdk 5.0.0's ProgressLogger for one thread: the same lines at the same records and the same return values.
     * htsjdk's record(SAMRecord) is synchronized and decodes every record's read name; this one reads the name only
     * for a line that prints it.
     */
    private static final class SingleThreadProgressLogger implements ProgressLoggerInterface {
        // AbstractProgressLogger prints the last read name once more than this many positions went backwards.
        private static final long PRINT_READ_NAME_THRESHOLD = 1000;

        // Takes each line's parts as Log.info does.
        private final Consumer<Object[]> sink;
        private final int n;
        private final String verb;
        private final NumberFormat fmt = new DecimalFormat("#,###");
        private final NumberFormat timeFmt = new DecimalFormat("00");
        private long startTime;
        private long processed;
        // The count at which the next line is logged: the next multiple of n.
        private long nextLine;
        private long lastStartTime;
        private String lastChrom;
        private int lastPos;
        private long countNonIncreasing;

        SingleThreadProgressLogger(final Log log, final int n, final String verb) {
            this(log::info, n, verb);
        }

        /** With each line's parts, as Log.info takes them, given to sink rather than logged. */
        SingleThreadProgressLogger(final Consumer<Object[]> sink, final int n, final String verb) {
            if (n <= 0) {
                throw new IllegalArgumentException("n must be positive: " + n);
            }
            this.sink = sink;
            this.n = n;
            this.verb = verb;
            reset();
        }

        @Override
        public boolean record(final String chrom, final int pos) {
            return record(chrom, pos, null);
        }

        @Override
        public boolean record(final SAMRecord rec) {
            // A BAM record with l_read_name 0 has a name that fails to decode; htsjdk decodes it here and throws.
            if (rec instanceof BAMRecord && rec.getReadNameLength() < 0) {
                rec.getReadName();
            }
            final String chrom = rec.getReferenceName();
            return SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(chrom) ? record(null, 0, rec) :
                    record(chrom, rec.getAlignmentStart(), rec);
        }

        @Override
        public boolean record(final SAMRecord... recs) {
            boolean logged = false;
            for (final SAMRecord rec : recs) {
                logged = record(rec) || logged;
            }
            return logged;
        }

        @Override
        public void reset() {
            startTime = System.currentTimeMillis();
            processed = 0;
            nextLine = n;
            lastStartTime = -1;
            lastChrom = null;
            lastPos = 0;
            countNonIncreasing = 0;
        }

        private boolean record(final String chrom, final int pos, final SAMRecord rec) {
            if (chrom != null && chrom.equals(lastChrom) && pos < lastPos) {
                countNonIncreasing++;
            } else {
                lastChrom = chrom;
            }
            lastPos = pos;
            if (lastStartTime == -1) {
                lastStartTime = System.currentTimeMillis();
            }
            if (++processed != nextLine) {
                return false;
            }
            nextLine += n;
            logLine(rec == null || countNonIncreasing <= PRINT_READ_NAME_THRESHOLD ? null : rec.getReadName());
            return true;
        }

        private void logLine(final String readName) {
            final long now = System.currentTimeMillis();
            final long lastPeriodSeconds = (now - lastStartTime) / 1000;
            lastStartTime = now;
            final long seconds = (now - startTime) / 1000;
            final long minutes = seconds / 60;
            final String elapsed = timeFmt.format(minutes / 60) + ":" + timeFmt.format(minutes % 60) + ":" +
                    timeFmt.format(seconds % 60);
            final String readInfo = lastChrom == null ? "*/*" : lastChrom + ":" + fmt.format(lastPos);
            final String readNameInfo = readName == null ? "" : ".  Last read name: " + readName;
            sink.accept(new Object[]{verb, " ", pad(fmt.format(processed), 13), " records.  Elapsed time: ", elapsed,
                    "s.  Time for last ", fmt.format(n), ": ", pad(fmt.format(lastPeriodSeconds), 4),
                    "s.  Last read position: ", readInfo, readNameInfo});
        }

        private static String pad(final String in, final int length) {
            final StringBuilder out = new StringBuilder(Math.max(length, in.length()));
            while (out.length() < length - in.length()) {
                out.append(' ');
            }
            return out.append(in).toString();
        }
    }

    /** What buildSortedReadEndLists looks up for one record, looked up on a decoding thread. */
    private static final class DecodedRead {
        // A record that is not a mapped primary read and has no PG tag (or pgIdsSeen is not collected).
        static final DecodedRead NOTHING = new DecodedRead(null, null, null, -1);

        // The PG tag; null when it is absent or pgIdsSeen is not collected.
        final String pg;
        // A mapped primary read's ends, as decodeReadEnds builds them; null for other records.
        final ReadEndsForMarkDuplicates ends;
        // The pair map's key, read group ID then read name, hashed; null when the read has no mapped mate or no
        // read group in the header.
        final String pairKey;
        // For a mapped primary read, its read group's first place in the header's read-group list (ReadGroupSlots);
        // -1 for other records and when the RG tag is absent or names no read group in the header.
        final int readGroupSlot;

        DecodedRead(final String pg, final ReadEndsForMarkDuplicates ends, final String pairKey, final int readGroupSlot) {
            this.pg = pg;
            this.ends = ends;
            this.pairKey = pairKey;
            this.readGroupSlot = readGroupSlot;
        }
    }

    /**
     * One decoding thread's map from read-group ID to the ID's first place in the records' header's read-group list,
     * made from the first record's header: every record of a ParallelBAMRecordIterator carries the same header.
     */
    private static final class ReadGroupSlots {
        private Map<String, Integer> slots = null;

        /** The ID's first place in the header's read-group list; -1 when the header has no read group with it. */
        int slot(final SAMFileHeader header, final String readGroupId) {
            if (slots == null) {
                slots = new HashMap<>();
                final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
                for (int i = 0; i < readGroups.size(); i++) {
                    final String id = readGroups.get(i).getReadGroupId();
                    if (header.getReadGroup(id) != null) {
                        slots.putIfAbsent(id, i);
                    }
                }
            }
            final Integer slot = slots.get(readGroupId);
            return slot == null ? -1 : slot;
        }
    }

    /**
     * The lookups buildSortedReadEndLists makes for one record, on a decoding thread: the PG tag when collectPgIds,
     * and for a mapped primary read its ends (decodeReadEnds), the pair map's key, with the key's hash computed, and
     * its read group's slot, all from one lookup of the RG tag.
     * A lookup that throws makes the decoding thread drop the whole result, and the main thread then repeats every
     * lookup itself, so a failure surfaces as it does there.
     */
    private DecodedRead decodeRead(final SAMRecord rec, final ReadNameParser parser, final boolean collectPgIds,
                                   final ReadGroupSlots readGroupSlots) {
        final String pg = collectPgIds ? rec.getStringAttribute(SAMTag.PG.name()) : null;
        if (pg != null) {
            pg.hashCode();
        }
        if (rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) {
            return pg == null ? DecodedRead.NOTHING : new DecodedRead(pg, null, null, -1);
        }
        final String readGroupId = (String) rec.getAttribute(ReservedTagConstants.READ_GROUP_ID);
        final int slot = readGroupId == null ? -1 : readGroupSlots.slot(rec.getHeader(), readGroupId);
        final ReadEndsForMarkDuplicates ends = decodeReadEnds(rec, parser, slot);
        String pairKey = null;
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            if (slot >= 0) {
                pairKey = readGroupId + rec.getReadName();
                pairKey.hashCode();
            } else {
                final SAMReadGroupRecord readGroup = rec.getReadGroup();
                if (readGroup != null) {
                    pairKey = readGroup.getReadGroupId() + rec.getReadName();
                    pairKey.hashCode();
                }
            }
        }
        return new DecodedRead(pg, ends, pairKey, slot);
    }

    /**
     * buildReadEnds for a mapped primary read on a decoding thread, but for libraryId and read1IndexInFile, which
     * depend on the reads before it. parser is the thread's own, without a log, so it never warns; readGroup stays -1
     * when it finds no location. The records' header is the one buildSortedReadEndLists reads. readGroupSlot is
     * DecodedRead's; at -1 the read group's place is looked up as buildReadEnds does.
     */
    private ReadEndsForMarkDuplicates decodeReadEnds(final SAMRecord rec, final ReadNameParser parser,
                                                     final int readGroupSlot) {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        ends.read1ReferenceIndex = rec.getReferenceIndex();
        ends.read1Coordinate = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();
        ends.orientation = rec.getReadNegativeStrandFlag() ? ReadEnds.R : ReadEnds.F;
        ends.score = DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            ends.read2ReferenceIndex = rec.getMateReferenceIndex();
        }
        if (parser.addLocationInformation(rec.getReadName(), ends)) {
            ends.readGroup = readGroupSlot >= 0 ? (short) readGroupSlot : readGroupIndex(rec.getHeader(), rec);
        }
        return ends;
    }

    /**
     * Goes through the accumulated ReadEndsForMarkDuplicates objects and determines which of them are
     * to be marked as duplicates.
     */
    protected void sortIndicesForDuplicates(final boolean indexOpticalDuplicates){
        final int maxInMemory = duplicateIndicesInMemory(indexOpticalDuplicates);
        if (indexOpticalDuplicates) {
            this.opticalDuplicateIndexes = new SortingLongCollection(maxInMemory, TMP_DIR.toArray(new File[TMP_DIR.size()]));
        }
        log.info("Will retain up to " + maxInMemory + " duplicate indices before spilling to disk.");
        this.duplicateIndexes = new SortingLongCollection(maxInMemory, TMP_DIR.toArray(new File[TMP_DIR.size()]));
        if (TAG_DUPLICATE_SET_MEMBERS) {
            final RepresentativeReadIndexerCodec representativeIndexCodec = new RepresentativeReadIndexerCodec();
            this.representativeReadIndicesForDuplicates = SortingCollection.newInstance(RepresentativeReadIndexer.class,
                    representativeIndexCodec,
                    Comparator.comparing(read -> read.readIndexInFile),
                    maxInMemory,
                    TMP_DIR);
        }
    }
    /** The duplicate indices sortIndicesForDuplicates' collections hold before spilling, as it logs. */
    private int duplicateIndicesInMemory(final boolean indexOpticalDuplicates) {
        final int entryOverhead;
        if (TAG_DUPLICATE_SET_MEMBERS) {
            // Memory requirements for RepresentativeReadIndexer:
            // three int entries + overhead: (3 * 4) + 4 = 16 bytes
            entryOverhead = 16;
        } else {
            entryOverhead = SortingLongCollection.SIZEOF;
        }
        // Keep this number from getting too large even if there is a huge heap.
        int maxInMemory = (int) Math.min((Runtime.getRuntime().maxMemory() * 0.25) / entryOverhead, (double) (Integer.MAX_VALUE - 32));
        // If we're also tracking optical duplicates, reduce maxInMemory, since we'll need two sorting collections
        if (indexOpticalDuplicates) {
            maxInMemory /= ((entryOverhead + SortingLongCollection.SIZEOF) / entryOverhead);
        }
        // The collections allocate and zero maxInMemory entries up front, so cap it at what can be added.
        return (int) Math.max(1, Math.min(maxInMemory, this.numRecordsInPass1));
    }

    public void generateDuplicateIndexes(final boolean useBarcodes, final boolean indexOpticalDuplicates) {
        // Before sortIndicesForDuplicates allocates its large arrays: the rows are not needed once the pair pass has
        // handled every pair, and G1 cannot move the keys' runs out of the way of a large allocation.
        final PairPassDuringRead pairPass = this.pairPassDuringRead;
        this.pairPassDuringRead = null;
        final boolean pairPassFinished = pairPass != null && pairPass.finish();
        if (pairPassFinished && this.pairedEndRows != null) {
            this.pairedEndRows = null;
            log.debug("Pair rows released before the duplicate indices are allocated.");
        }
        if (this.pairedEndKeys != null) {
            this.pairedEndKeys.splitRuns();
            log.debug("Paired-end keys split before the duplicate indices are allocated.");
        }
        sortIndicesForDuplicates(indexOpticalDuplicates);

        ReadEndsForMarkDuplicates firstOfNextChunk = null;
        final List<ReadEndsForMarkDuplicates> nextChunk = new ArrayList<>(200);

        // First just do the pairs
        log.info("Traversing read pair information and detecting duplicates.");
        if (pairPassFinished) {
            pairPass.addResults();
        } else if (this.pairedEndRows != null) {
            // The chunks of the loop below, in the same order; a chunk of one pair is only counted.
            final PairedReadEndRows rows = this.pairedEndRows;
            this.pairedEndRows = null;
            rows.sort();
            final int numberOfRows = rows.size();
            for (int start = 0, end; start < numberOfRows; start = end) {
                end = start + 1;
                while (end < numberOfRows && !rows.startsDuplicateSet(end)) {
                    end++;
                }
                if (end - start == 1) {
                    addSingletonToCount(libraryIdGenerator);
                } else {
                    nextChunk.clear();
                    for (int position = start; position < end; position++) {
                        nextChunk.add(rows.getAt(position));
                    }
                    handleChunk(nextChunk);
                }
            }
        } else {
            for (final ReadEndsForMarkDuplicates next : this.pairSort) {
                if (firstOfNextChunk != null && areComparableForDuplicates(firstOfNextChunk, next, true, useBarcodes)) {
                    nextChunk.add(next);
                } else {
                    handleChunk(nextChunk);
                    nextChunk.clear();
                    nextChunk.add(next);
                    firstOfNextChunk = next;
                }
            }
            handleChunk(nextChunk);

            this.pairSort.cleanup();
            this.pairSort = null;
        }

        // Now deal with the fragments
        log.info("Traversing fragment information and detecting duplicates.");
        boolean containsPairs = false;
        boolean containsFrags = false;

        firstOfNextChunk = null;
        nextChunk.clear();

        // A chunk's paired ends are those in fragSort plus those whose key pairedEndKeys holds.
        for (final ReadEndsForMarkDuplicates next : this.fragSort) {
            if (firstOfNextChunk != null && areComparableForDuplicates(firstOfNextChunk, next, false, useBarcodes)) {
                nextChunk.add(next);
                containsPairs = containsPairs || next.isPaired();
                containsFrags = containsFrags || !next.isPaired();
            } else {
                if ((nextChunk.size() > 1 || containsPairs) && containsFrags) {
                    markDuplicateFragments(nextChunk, containsPairs);
                }
                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                containsPairs = next.isPaired() || this.pairedEndKeys != null && this.pairedEndKeys.contains(next);
                containsFrags = !next.isPaired();
            }
        }
        markDuplicateFragments(nextChunk, containsPairs);
        this.fragSort.cleanup();
        this.fragSort = null;
        this.pairedEndKeys = null;

        log.info("Sorting list of duplicate records.");
        this.duplicateIndexes.doneAddingStartIteration();
        if (this.opticalDuplicateIndexes != null) {
            this.opticalDuplicateIndexes.doneAddingStartIteration();
        }
        if (TAG_DUPLICATE_SET_MEMBERS) {
            this.representativeReadIndicesForDuplicates.doneAdding();
        }
    }

    protected void handleChunk(List<ReadEndsForMarkDuplicates> nextChunk) {
        if (nextChunk.size() > 1) {
            markDuplicatePairs(nextChunk);
            if (TAG_DUPLICATE_SET_MEMBERS) {
                addRepresentativeReadIndex(nextChunk);
            }
        } else if (nextChunk.size() == 1) {
            addSingletonToCount(libraryIdGenerator);
        }
    }

    protected boolean areComparableForDuplicates(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs, final boolean compareRead2, final boolean useBarcodes) {
        boolean areComparable = lhs.libraryId == rhs.libraryId;

        if (useBarcodes && areComparable) { // areComparable is useful here to avoid the casts below
            final ReadEndsForMarkDuplicatesWithBarcodes lhsWithBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) lhs;
            final ReadEndsForMarkDuplicatesWithBarcodes rhsWithBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) rhs;
            areComparable = lhsWithBarcodes.barcode == rhsWithBarcodes.barcode &&
                    lhsWithBarcodes.readOneBarcode == rhsWithBarcodes.readOneBarcode &&
                    lhsWithBarcodes.readTwoBarcode == rhsWithBarcodes.readTwoBarcode;
        }

        if (areComparable) {
            areComparable = lhs.read1ReferenceIndex == rhs.read1ReferenceIndex &&
                    lhs.read1Coordinate == rhs.read1Coordinate &&
                    lhs.orientation == rhs.orientation;
        }

        if (areComparable && compareRead2) {
            areComparable = lhs.read2ReferenceIndex == rhs.read2ReferenceIndex &&
                    lhs.read2Coordinate == rhs.read2Coordinate;
        }

        return areComparable;
    }

    private void addIndexAsDuplicate(final long bamIndex) {
        this.duplicateIndexes.add(bamIndex);
        ++this.numDuplicateIndices;
    }

    private void addRepresentativeReadOfDuplicateSet(final long representativeReadIndexInFile, final int setSize, final long read1IndexInFile) {
        final RepresentativeReadIndexer rri = new RepresentativeReadIndexer();
        rri.representativeReadIndexInFile = (int) representativeReadIndexInFile;
        rri.setSize = setSize;
        rri.readIndexInFile = (int) read1IndexInFile;
        this.representativeReadIndicesForDuplicates.add(rri);
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and identify the representative read based on
     * quality score. For all members of the duplicate set, add the read1 index-in-file of the representative
     * read to the records of the first and second in a pair. This value becomes is used for
     * the 'DI' tag.
     */
    private void addRepresentativeReadIndex(final List<ReadEndsForMarkDuplicates> list) {
        short maxScore = 0;
        ReadEndsForMarkDuplicates best = null;

        /** All read ends should have orientation FF, FR, RF, or RR **/
        for (final ReadEndsForMarkDuplicates end : list) {
            if (end.score > maxScore || best == null) {
                maxScore = end.score;
                best = end;
            }
        }

        // for read name (for representative read name), add the last of the pair that was examined
        for (final ReadEndsForMarkDuplicates end : list) {
            addRepresentativeReadOfDuplicateSet(best.read1IndexInFile, list.size(), end.read1IndexInFile);
            if (end.read1IndexInFile != end.read2IndexInFile) {
                addRepresentativeReadOfDuplicateSet(best.read1IndexInFile, list.size(), end.read2IndexInFile);
            }
        }
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and removes from it all objects that should
     * not be marked as duplicates.  This assumes that the list contains objects representing pairs.
     */
    private void markDuplicatePairs(final List<ReadEndsForMarkDuplicates> list) {
        markDuplicatePairs(list, this.libraryIdGenerator, this::addIndexAsDuplicate,
                this.opticalDuplicateIndexes == null ? null : this.opticalDuplicateIndexes::add);
    }

    /**
     * markDuplicatePairs with the counts going to counts' histograms, each duplicate's indices to duplicates and each
     * optical duplicate's to opticalDuplicates, unless it is null.
     */
    private void markDuplicatePairs(final List<ReadEndsForMarkDuplicates> list, final LibraryIdGenerator counts,
                                    final LongConsumer duplicates, final LongConsumer opticalDuplicates) {
        short maxScore = 0;
        ReadEndsForMarkDuplicates best = null;

        /** All read ends should have orientation FF, FR, RF, or RR **/
        for (final ReadEndsForMarkDuplicates end : list) {
            if (end.score > maxScore || best == null) {
                maxScore = end.score;
                best = end;
            }
        }

        if (this.READ_NAME_REGEX != null) {
            AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(list, best, opticalDuplicateFinder, counts);
        }

        for (final ReadEndsForMarkDuplicates end : list) {
            if (end != best) {
                duplicates.accept(end.read1IndexInFile);

                // in query-sorted case, these will be the same.
                // TODO: also in coordinate sorted, when one read is unmapped
                if (end.read2IndexInFile != end.read1IndexInFile) {
                    duplicates.accept(end.read2IndexInFile);
                }

                if (end.isOpticalDuplicate && opticalDuplicates != null) {
                    opticalDuplicates.accept(end.read1IndexInFile);
                    // We expect end.read2IndexInFile==read1IndexInFile when we are in queryname sorted files, as the read-pairs
                    // will be sorted together and nextIndexIfNeeded() will only pull one index from opticalDuplicateIndexes.
                    // This means that in queryname sorted order we will only pull from the sorting collection once,
                    // where as we would pull twice for coordinate sorted files.
                    if (end.read2IndexInFile != end.read1IndexInFile) {
                        opticalDuplicates.accept(end.read2IndexInFile);
                    }
                }
            }
        }
    }

    /**
     * Method for deciding when to pull from the SortingLongCollection for the next read based on sorting order.
     * - If file is queryname sorted then we expect one index per pair of reads, so we only want to iterate when we
     * are no longer reading from that read-pair.
     * - If file is coordinate-sorted we want to base our iteration entirely on the indexes of both reads in the pair
     * <p>
     * This logic is applied to both Optical and Library duplicates
     *
     * @param sortOrder          Sort order for the underlying bam file
     * @param recordInFileIndex  Index of the current sam record rec
     * @param nextDuplicateIndex Index of next expected duplicate (optical or otherwise) in the file
     * @param lastQueryName      Name of the last read seen (for keeping queryname sorted groups together)
     * @param rec                Current record to compare against
     * @param duplicateIndexes   DuplicateIndexes collection to iterate over
     * @return the duplicate after iteration
     */
    private long nextIndexIfNeeded(final SAMFileHeader.SortOrder sortOrder, final long recordInFileIndex, final long nextDuplicateIndex, final String lastQueryName, final SAMRecord rec, final SortingLongCollection duplicateIndexes) {
        // Manage the flagging of optical/sequencing duplicates
        // Possibly figure out the next opticalDuplicate index (if going by coordinate, if going by query name, only do this
        // if the query name has changed)
        final boolean needNextDuplicateIndex = recordInFileIndex > nextDuplicateIndex &&
                (sortOrder == SAMFileHeader.SortOrder.coordinate || !rec.getReadName().equals(lastQueryName));

        if (needNextDuplicateIndex) {
            return (duplicateIndexes.hasNext() ? duplicateIndexes.next() : NO_SUCH_INDEX);
        }
        return nextDuplicateIndex;
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and removes from it all objects that should
     * not be marked as duplicates.  This will set the duplicate index for only list items are fragments.
     *
     * @param containsPairs true if the list also contains objects containing pairs, false otherwise.
     */
    protected void markDuplicateFragments(final List<ReadEndsForMarkDuplicates> list, final boolean containsPairs) {
        if (containsPairs) {
            for (final ReadEndsForMarkDuplicates end : list) {
                if (!end.isPaired()) {
                    addIndexAsDuplicate(end.read1IndexInFile);
                }
            }
        } else {
            short maxScore = 0;
            ReadEndsForMarkDuplicates best = null;
            for (final ReadEndsForMarkDuplicates end : list) {
                if (end.score > maxScore || best == null) {
                    maxScore = end.score;
                    best = end;
                }
            }

            for (final ReadEndsForMarkDuplicates end : list) {
                if (end != best) {
                    addIndexAsDuplicate(end.read1IndexInFile);
                }
            }
        }
    }

    /** A reference index and coordinate as a long in the same order, the coordinate clamped to the int range. */
    private static long position(final int referenceIndex, final long coordinate) {
        final long clamped = Math.max(Integer.MIN_VALUE, Math.min(Integer.MAX_VALUE, coordinate));
        return ((long) referenceIndex << 32) + (clamped - Integer.MIN_VALUE);
    }

    /** Completed pairs handed to the pair pass during pass 1, with the frontier main gives them. */
    private static final class PairBatch {
        static final PairBatch STOP = new PairBatch(0, null);

        final ReadEndsForMarkDuplicates[] ends;
        int size = 0;
        long frontier;
        boolean last;
        // The lowest and highest read2 position() of the pairs not yet handled, and their lowest record index, kept
        // by the pair pass.
        long lowest;
        long highest;
        long firstIndex;

        // With pass 2 beside pass 1, null otherwise: the unpaired ends of mapped primary reads, and the keys
        // (SecondPassDuringRead.key) of the ends of reads with a mapped mate, with the keys' lowest and highest
        // position().
        final ReadEndsForMarkDuplicates[] fragments;
        int numFragments = 0;
        long[] keys;
        int numKeys = 0;
        long keysLowest = Long.MAX_VALUE;
        long keysHighest = Long.MIN_VALUE;
        // Every record below this index has been read, and each of them whose duplicate flag is not yet decided has
        // its end in this batch or an earlier one.
        long undecidedFrom;

        PairBatch(final int capacity, final long[] keys) {
            this.ends = new ReadEndsForMarkDuplicates[capacity];
            this.fragments = keys == null ? null : new ReadEndsForMarkDuplicates[capacity];
            this.keys = keys;
        }
    }

    /** Unpaired ends waiting on the pair pass for pass 2, with their lowest and highest position() and record index. */
    private static final class FragmentBatch {
        final ReadEndsForMarkDuplicates[] ends;
        int size;
        long lowest;
        long highest;
        long firstIndex;

        FragmentBatch(final ReadEndsForMarkDuplicates[] ends, final int size) {
            this.ends = ends;
            this.size = size;
            measure();
        }

        void measure() {
            long low = Long.MAX_VALUE;
            long high = Long.MIN_VALUE;
            long first = Long.MAX_VALUE;
            for (int i = 0; i < this.size; i++) {
                final ReadEndsForMarkDuplicates end = this.ends[i];
                final long at = position(end.read1ReferenceIndex, end.read1Coordinate);
                low = Math.min(low, at);
                high = Math.max(high, at);
                first = Math.min(first, end.read1IndexInFile);
            }
            this.lowest = low;
            this.highest = high;
            this.firstIndex = first;
        }
    }

    /** Keys of paired read ends, kept by the pair pass while an unpaired end with the same key may still come. */
    private static final class KeyBatch {
        final long[] keys;
        final int size;
        final long lowest;
        final long highest;
        private boolean sorted = false;

        KeyBatch(final long[] keys, final int size, final long lowest, final long highest) {
            this.keys = keys;
            this.size = size;
            this.lowest = lowest;
            this.highest = highest;
        }

        /** Whether key, of an end at this position(), is here. */
        boolean contains(final long key, final long at) {
            if (at < this.lowest || at > this.highest) {
                return false;
            }
            if (!this.sorted) {
                Arrays.sort(this.keys, 0, this.size);
                this.sorted = true;
            }
            return Arrays.binarySearch(this.keys, 0, this.size, key) >= 0;
        }
    }

    /** Bits by record index in chunks made as bits are set; one thread sets them. */
    private static final class IndexBits {
        private static final int CHUNK_BITS = 20;
        private static final int WORDS = 1 << (CHUNK_BITS - 6);

        private long[][] chunks = new long[64][];

        void set(final long index) {
            final int c = (int) (index >>> CHUNK_BITS);
            if (c >= this.chunks.length) {
                this.chunks = Arrays.copyOf(this.chunks, Math.max(c + 1, 2 * this.chunks.length));
            }
            long[] chunk = this.chunks[c];
            if (chunk == null) {
                chunk = new long[WORDS];
                this.chunks[c] = chunk;
            }
            chunk[(int) (index >>> 6) & (WORDS - 1)] |= 1L << index;
        }

        boolean get(final long index) {
            final long[][] all = this.chunks;
            final int c = (int) (index >>> CHUNK_BITS);
            if (c >= all.length) {
                return false;
            }
            final long[] chunk = all[c];
            return chunk != null && (chunk[(int) (index >>> 6) & (WORDS - 1)] & (1L << index)) != 0;
        }
    }

    /**
     * The pair pass on a thread of its own during pass 1. Main hands it the completed pairs in batches, each with a
     * frontier: a position that every later read with a mapped mate must reach with its 5' end, which main checks.
     * Both 5' ends of a pair lie at or below its read2 end, so once a duplicate set's read2 end lies below a frontier
     * every pair of the set has been handed over. The thread then sorts such pairs with ReadEndsMDComparator and
     * handles each set as handleChunk does, into its own duplicate indices, singleton count and histograms, which
     * generateDuplicateIndexes adds to its own. A read below the frontier, a pair the rows refuse, too many pairs
     * waiting or a failure on the thread abandons the result, and the pair pass runs on the rows after pass 1.
     *
     * With pass 2 beside pass 1 the thread also decides its duplicate flags: the pairs' as above, and those of the
     * unpaired ends, which main hands over too and checks against the frontier, each chunk of them once its 5' position
     * lies below a frontier, as the fragment pass would, with a paired end with the same key (which main hands over
     * as well) making it a chunk with pairs. After each batch it tells pass 2 below which record every flag is decided.
     */
    private final class PairPassDuringRead implements Runnable {
        private static final int BATCH = 4096;
        // How far below a read's alignment start the frontier it gives lies: room for a later read's leading clips,
        // and small beside a contig, so that a deep short contig such as chrM is handled while it is read.
        private static final int MARGIN = 1000;
        // With pass 2 beside pass 1, a batch also goes once it holds BATCH unpaired ends or KEYS keys, or RECORDS
        // records were read since the last one went.
        private static final int KEYS = 4 * BATCH;
        private static final int RECORDS = 4 * BATCH;

        private final ArrayBlockingQueue<PairBatch> queue;
        private final long maxWaitingPairs;
        private final Thread thread;
        private volatile boolean abandoned = false;
        // Pass 2 beside pass 1; null without it.
        private final SecondPassDuringRead secondPass;
        // Key arrays the thread is done with, for main's next batches.
        private final ConcurrentLinkedQueue<long[]> freeKeys = new ConcurrentLinkedQueue<>();

        // Main's: the batch being filled, the last frontier handed over and the records read when it was.
        private PairBatch batch;
        private long frontier = Long.MIN_VALUE;
        private long recordsAtLastSend = 0;

        // The thread's, read by main after the join.
        private boolean handledAll = false;
        private final LibraryIdGenerator counts = new LibraryIdGenerator(new SAMFileHeader());
        private long singletons = 0;
        private long[] duplicates = new long[1024];
        private int numDuplicates = 0;
        private long[] opticalDuplicates = new long[16];
        private int numOpticalDuplicates = 0;

        // The thread's: batches holding pairs not yet handled, how many, and the lowest read2 position among them.
        private final List<PairBatch> waiting = new ArrayList<>();
        private long waitingPairs = 0;
        private long lowestWaiting = Long.MAX_VALUE;
        private ReadEndsForMarkDuplicates[] ready = new ReadEndsForMarkDuplicates[BATCH];
        private final List<ReadEndsForMarkDuplicates> set = new ArrayList<>(200);
        private final Comparator<ReadEndsForMarkDuplicates> comparator = new ReadEndsMDComparator(false);
        private final LongConsumer duplicateSink = this::addDuplicate;
        private final LongConsumer opticalDuplicateSink = this::addOpticalDuplicate;

        // The thread's, for pass 2: unpaired ends not yet decided, how many, the lowest position among them, and the
        // keys of paired ends that an unpaired end not yet handed over may still share.
        private final List<FragmentBatch> waitingFragments = new ArrayList<>();
        private long waitingFragmentEnds = 0;
        private long lowestWaitingFragment = Long.MAX_VALUE;
        private final List<KeyBatch> keyBatches = new ArrayList<>();
        private long heldKeys = 0;
        private ReadEndsForMarkDuplicates[] readyFragments = new ReadEndsForMarkDuplicates[BATCH];

        /**
         * @param capacity   the rows' capacity; up to an eighth of it may wait on the thread and in its queue each
         * @param secondPass pass 2 beside pass 1, whose duplicate flags the thread decides; null for none
         */
        PairPassDuringRead(final int capacity, final SecondPassDuringRead secondPass) {
            this.maxWaitingPairs = Math.max(4L * BATCH, capacity / 8);
            this.queue = new ArrayBlockingQueue<>((int) Math.min(4096, maxWaitingPairs / BATCH));
            this.secondPass = secondPass;
            this.batch = newBatch();
            this.thread = new Thread(this, "MarkDuplicates pair pass");
            this.thread.setDaemon(true);
            this.thread.start();
        }

        /** A batch to fill, with room for pass 2's ends and keys while pass 2 goes on. */
        private PairBatch newBatch() {
            if (this.secondPass == null || !this.secondPass.isActive()) {
                return new PairBatch(BATCH, null);
            }
            final long[] keys = this.freeKeys.poll();
            return new PairBatch(BATCH, keys != null ? keys : new long[KEYS]);
        }

        /** Whether a read's 5' end, from its ends, is at or above the last frontier handed over. */
        boolean atOrAboveFrontier(final ReadEndsForMarkDuplicates end) {
            return position(end.read1ReferenceIndex, end.read1Coordinate) >= frontier;
        }

        /** Adds a completed pair to the batch. */
        void add(final ReadEndsForMarkDuplicates pair) {
            this.batch.ends[this.batch.size++] = pair;
        }

        /** For pass 2: the end of a mapped primary read without a mapped mate; one below the frontier abandons it. */
        void addFragment(final ReadEndsForMarkDuplicates end) {
            final PairBatch filling = this.batch;
            if (filling.fragments == null || !atOrAboveFrontier(end)) {
                this.secondPass.abandon();
                return;
            }
            filling.fragments[filling.numFragments++] = end;
        }

        /** For pass 2: the end of a read with a mapped mate; one whose key cannot be packed abandons it. */
        void addPairedKey(final ReadEndsForMarkDuplicates end) {
            final PairBatch filling = this.batch;
            // A read whose mate is mapped on no reference is not paired for the fragment pass, which may then mark it.
            final long key = filling.keys == null || !end.isPaired() ? -1 : this.secondPass.key(end);
            if (key < 0) {
                this.secondPass.abandon();
                return;
            }
            filling.keys[filling.numKeys++] = key;
            final long at = position(end.read1ReferenceIndex, end.read1Coordinate);
            filling.keysLowest = Math.min(filling.keysLowest, at);
            filling.keysHighest = Math.max(filling.keysHighest, at);
        }

        /**
         * Whether the batch goes once recordsRead records are read: when it holds BATCH pairs, or with pass 2 beside
         * pass 1 when an array is full or RECORDS records were read since the last batch went.
         */
        boolean readyToSend(final long recordsRead) {
            final PairBatch filling = this.batch;
            return filling.size >= BATCH || filling.keys != null && (filling.numFragments >= BATCH ||
                    filling.numKeys >= filling.keys.length || recordsRead - this.recordsAtLastSend >= RECORDS);
        }

        /**
         * Hands over the batch with a frontier MARGIN below the alignment start of rec, the last record read, unless
         * the last one was higher, and waits while the queue is full.
         *
         * @param recordsRead the records read, rec the last of them
         * @return false when the pass is abandoned; nothing more is handed over then
         */
        boolean send(final SAMRecord rec, final long recordsRead) {
            this.frontier = Math.max(this.frontier,
                    position(rec.getReferenceIndex(), (long) rec.getAlignmentStart() - MARGIN));
            return send(false, recordsRead);
        }

        /** Hands over the last pairs: pass 1 has read every read. */
        void doneAdding() {
            this.frontier = Long.MAX_VALUE;
            send(true, Long.MAX_VALUE);
        }

        private boolean send(final boolean last, final long recordsRead) {
            this.batch.frontier = this.frontier;
            this.batch.last = last;
            if (this.batch.keys != null && !last) {
                this.batch.undecidedFrom = this.secondPass.firstUndecided(recordsRead, this.frontier);
                this.secondPass.validated.advance(recordsRead);
            } else if (this.batch.keys != null) {
                // Main has given pass 2's reader the final count of records pass 1 validated.
                this.batch.undecidedFrom = Long.MAX_VALUE;
            }
            this.recordsAtLastSend = recordsRead;
            try {
                while (!this.queue.offer(this.batch, 1, TimeUnit.MILLISECONDS)) {
                    if (this.abandoned) {
                        break;
                    }
                }
            } catch (final InterruptedException e) {
                Thread.currentThread().interrupt();
                abandon();
            }
            if (this.abandoned) {
                abandon();
                return false;
            }
            this.batch = newBatch();
            return true;
        }

        /** Stops the thread and drops its result, and abandons pass 2 beside it; called on main or the thread. */
        void abandon() {
            this.abandoned = true;
            this.queue.clear();
            this.queue.offer(PairBatch.STOP);
            if (this.secondPass != null) {
                this.secondPass.abandon();
            }
        }

        /**
         * Whether the thread handled every pair, so that its results stand, waiting for it unless it was abandoned:
         * an abandoned thread may still be in a large duplicate set, and touches nothing but its own results.
         */
        boolean finish() {
            if (this.abandoned) {
                return false;
            }
            boolean interrupted = false;
            while (true) {
                try {
                    this.thread.join();
                    break;
                } catch (final InterruptedException e) {
                    interrupted = true;
                }
            }
            if (interrupted) {
                Thread.currentThread().interrupt();
            }
            return this.handledAll && !this.abandoned;
        }

        /** Adds the results to the duplicate indices and the histograms, as the pair pass would have. */
        void addResults() {
            for (int i = 0; i < this.numDuplicates; i++) {
                addIndexAsDuplicate(this.duplicates[i]);
            }
            if (opticalDuplicateIndexes != null) {
                for (int i = 0; i < this.numOpticalDuplicates; i++) {
                    opticalDuplicateIndexes.add(this.opticalDuplicates[i]);
                }
            }
            addHistograms();
        }

        /** Adds the results to the histograms, as the pair pass would have. */
        void addHistograms() {
            libraryIdGenerator.getDuplicateCountHist().addHistogram(this.counts.getDuplicateCountHist());
            libraryIdGenerator.getNonOpticalDuplicateCountHist().addHistogram(this.counts.getNonOpticalDuplicateCountHist());
            libraryIdGenerator.getOpticalDuplicateCountHist().addHistogram(this.counts.getOpticalDuplicateCountHist());
            libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap().addHistogram(this.counts.getOpticalDuplicatesByLibraryIdMap());
            // Each set of one pair, as addSingletonToCount counts it.
            if (this.singletons > 0) {
                libraryIdGenerator.getDuplicateCountHist().increment(1.0, this.singletons);
                libraryIdGenerator.getNonOpticalDuplicateCountHist().increment(1.0, this.singletons);
            }
        }

        /** The optical duplicate clusters the thread found; read after finish(). */
        long numberOfOpticalDuplicateClusters() {
            return (long) this.counts.getOpticalDuplicatesByLibraryIdMap().getSumOfValues();
        }

        @Override
        public void run() {
            try {
                while (!this.abandoned) {
                    final PairBatch received = this.queue.take();
                    if (received == PairBatch.STOP || this.abandoned) {
                        break;
                    }
                    receive(received);
                    if (received.last && !this.abandoned) {
                        this.handledAll = true;
                        break;
                    }
                }
            } catch (final InterruptedException | RuntimeException | Error e) {
                abandon();
            }
            this.waiting.clear();
            this.ready = null;
            this.waitingFragments.clear();
            this.keyBatches.clear();
            this.readyFragments = null;
        }

        private void receive(final PairBatch received) {
            long lowest = Long.MAX_VALUE;
            long highest = Long.MIN_VALUE;
            long firstIndex = Long.MAX_VALUE;
            for (int i = 0; i < received.size; i++) {
                final ReadEndsForMarkDuplicates end = received.ends[i];
                final long read2 = position(end.read2ReferenceIndex, end.read2Coordinate);
                lowest = Math.min(lowest, read2);
                highest = Math.max(highest, read2);
                firstIndex = Math.min(firstIndex, Math.min(end.read1IndexInFile, end.read2IndexInFile));
            }
            received.lowest = lowest;
            received.highest = highest;
            received.firstIndex = firstIndex;
            if (received.size > 0) {
                this.waiting.add(received);
                this.waitingPairs += received.size;
                this.lowestWaiting = Math.min(this.lowestWaiting, lowest);
            }
            if (this.waitingPairs > this.maxWaitingPairs) {
                abandon();
            } else if (received.last || received.frontier > this.lowestWaiting) {
                handleBelow(received.frontier, received.last);
            }
            if (received.keys != null) {
                if (!this.abandoned && this.secondPass.isActive()) {
                    decideForSecondPass(received);
                } else {
                    this.waitingFragments.clear();
                    this.keyBatches.clear();
                }
                received.keys = null;
            }
        }

        /** Handles the sets of the waiting pairs whose read2 end lies below the frontier, or of all of them. */
        private void handleBelow(final long frontier, final boolean all) {
            if (this.ready.length < this.waitingPairs) {
                this.ready = new ReadEndsForMarkDuplicates[(int) Math.max(this.waitingPairs, 2L * this.ready.length)];
            }
            final ReadEndsForMarkDuplicates[] ends = this.ready;
            int n = 0;
            int kept = 0;
            long lowest = Long.MAX_VALUE;
            for (int b = 0; b < this.waiting.size(); b++) {
                final PairBatch waitingBatch = this.waiting.get(b);
                if (all || waitingBatch.highest < frontier) {
                    System.arraycopy(waitingBatch.ends, 0, ends, n, waitingBatch.size);
                    n += waitingBatch.size;
                    continue;
                }
                if (waitingBatch.lowest < frontier) {
                    int keep = 0;
                    long low = Long.MAX_VALUE;
                    long high = Long.MIN_VALUE;
                    long first = Long.MAX_VALUE;
                    for (int i = 0; i < waitingBatch.size; i++) {
                        final ReadEndsForMarkDuplicates end = waitingBatch.ends[i];
                        final long read2 = position(end.read2ReferenceIndex, end.read2Coordinate);
                        if (read2 < frontier) {
                            ends[n++] = end;
                        } else {
                            waitingBatch.ends[keep++] = end;
                            low = Math.min(low, read2);
                            high = Math.max(high, read2);
                            first = Math.min(first, Math.min(end.read1IndexInFile, end.read2IndexInFile));
                        }
                    }
                    Arrays.fill(waitingBatch.ends, keep, waitingBatch.size, null);
                    waitingBatch.size = keep;
                    waitingBatch.lowest = low;
                    waitingBatch.highest = high;
                    waitingBatch.firstIndex = first;
                }
                this.waiting.set(kept++, waitingBatch);
                lowest = Math.min(lowest, waitingBatch.lowest);
            }
            this.waiting.subList(kept, this.waiting.size()).clear();
            this.waitingPairs -= n;
            this.lowestWaiting = lowest;

            Arrays.sort(ends, 0, n, this.comparator);
            for (int start = 0, end; start < n; start = end) {
                end = start + 1;
                while (end < n && areComparableForDuplicates(ends[start], ends[end], true, false)) {
                    end++;
                }
                if (end - start == 1) {
                    this.singletons++;
                } else {
                    this.set.clear();
                    for (int i = start; i < end; i++) {
                        this.set.add(ends[i]);
                    }
                    markDuplicatePairs(this.set, this.counts, this.duplicateSink, this.opticalDuplicateSink);
                }
            }
            this.set.clear();
            Arrays.fill(ends, 0, n, null);
        }

        /**
         * For pass 2, after the batch's pairs: keeps its unpaired ends and keys, decides the chunks of unpaired ends
         * below its frontier, drops the keys below it, and tells pass 2 below which record every flag is decided.
         */
        private void decideForSecondPass(final PairBatch received) {
            if (received.numFragments > 0) {
                final FragmentBatch fragments = new FragmentBatch(received.fragments, received.numFragments);
                this.waitingFragments.add(fragments);
                this.waitingFragmentEnds += fragments.size;
                this.lowestWaitingFragment = Math.min(this.lowestWaitingFragment, fragments.lowest);
            }
            if (received.numKeys > 0) {
                this.keyBatches.add(new KeyBatch(received.keys, received.numKeys, received.keysLowest,
                        received.keysHighest));
                this.heldKeys += received.numKeys;
            } else {
                this.freeKeys.offer(received.keys);
            }
            if (this.waitingFragmentEnds > this.maxWaitingPairs || this.heldKeys > 4 * this.maxWaitingPairs) {
                this.secondPass.abandon();
                this.waitingFragments.clear();
                this.keyBatches.clear();
                return;
            }
            if (received.last || received.frontier > this.lowestWaitingFragment) {
                handleFragmentsBelow(received.frontier, received.last);
            }
            // Every later read's 5' end is at or above the frontier, so keys below it match no end still to come.
            int kept = 0;
            for (int b = 0; b < this.keyBatches.size(); b++) {
                final KeyBatch keys = this.keyBatches.get(b);
                if (received.last || keys.highest < received.frontier) {
                    this.heldKeys -= keys.size;
                    this.freeKeys.offer(keys.keys);
                } else {
                    this.keyBatches.set(kept++, keys);
                }
            }
            this.keyBatches.subList(kept, this.keyBatches.size()).clear();

            long undecidedFrom = received.undecidedFrom;
            for (final PairBatch waitingBatch : this.waiting) {
                undecidedFrom = Math.min(undecidedFrom, waitingBatch.firstIndex);
            }
            for (final FragmentBatch waitingBatch : this.waitingFragments) {
                undecidedFrom = Math.min(undecidedFrom, waitingBatch.firstIndex);
            }
            this.secondPass.decided(undecidedFrom);
        }

        /**
         * Decides the chunks of the waiting unpaired ends whose 5' position lies below the frontier, or of all of
         * them, as the fragment pass does: every chunk with a paired end marks all of its unpaired ends, and any other
         * chunk all but its best.
         */
        private void handleFragmentsBelow(final long frontier, final boolean all) {
            if (this.readyFragments.length < this.waitingFragmentEnds) {
                this.readyFragments = new ReadEndsForMarkDuplicates[(int) Math.max(this.waitingFragmentEnds,
                        2L * this.readyFragments.length)];
            }
            final ReadEndsForMarkDuplicates[] ends = this.readyFragments;
            int n = 0;
            int kept = 0;
            long lowest = Long.MAX_VALUE;
            for (int b = 0; b < this.waitingFragments.size(); b++) {
                final FragmentBatch waitingBatch = this.waitingFragments.get(b);
                if (all || waitingBatch.highest < frontier) {
                    System.arraycopy(waitingBatch.ends, 0, ends, n, waitingBatch.size);
                    n += waitingBatch.size;
                    continue;
                }
                if (waitingBatch.lowest < frontier) {
                    int keep = 0;
                    for (int i = 0; i < waitingBatch.size; i++) {
                        final ReadEndsForMarkDuplicates end = waitingBatch.ends[i];
                        if (position(end.read1ReferenceIndex, end.read1Coordinate) < frontier) {
                            ends[n++] = end;
                        } else {
                            waitingBatch.ends[keep++] = end;
                        }
                    }
                    Arrays.fill(waitingBatch.ends, keep, waitingBatch.size, null);
                    waitingBatch.size = keep;
                    waitingBatch.measure();
                }
                this.waitingFragments.set(kept++, waitingBatch);
                lowest = Math.min(lowest, waitingBatch.lowest);
            }
            this.waitingFragments.subList(kept, this.waitingFragments.size()).clear();
            this.waitingFragmentEnds -= n;
            this.lowestWaitingFragment = lowest;

            Arrays.sort(ends, 0, n, this.comparator);
            for (int start = 0, end; start < n; start = end) {
                end = start + 1;
                while (end < n && areComparableForDuplicates(ends[start], ends[end], false, false)) {
                    end++;
                }
                if (heldKey(ends[start])) {
                    for (int i = start; i < end; i++) {
                        fragmentDuplicate(ends[i].read1IndexInFile);
                    }
                } else if (end - start > 1) {
                    short maxScore = 0;
                    ReadEndsForMarkDuplicates best = null;
                    for (int i = start; i < end; i++) {
                        if (ends[i].score > maxScore || best == null) {
                            maxScore = ends[i].score;
                            best = ends[i];
                        }
                    }
                    for (int i = start; i < end; i++) {
                        if (ends[i] != best) {
                            fragmentDuplicate(ends[i].read1IndexInFile);
                        }
                    }
                }
            }
            Arrays.fill(ends, 0, n, null);
        }

        /** Whether a paired end handed over has the key of this unpaired end. */
        private boolean heldKey(final ReadEndsForMarkDuplicates end) {
            final long key = this.secondPass.key(end);
            if (key < 0) {
                // Main abandons pass 2 on a paired end whose key cannot be packed, so none has this one's.
                return false;
            }
            final long at = position(end.read1ReferenceIndex, end.read1Coordinate);
            for (final KeyBatch keys : this.keyBatches) {
                if (keys.contains(key, at)) {
                    return true;
                }
            }
            return false;
        }

        private void fragmentDuplicate(final long index) {
            this.secondPass.duplicates.set(index);
            this.secondPass.fragmentDuplicates++;
        }

        private void addDuplicate(final long index) {
            if (this.numDuplicates == this.duplicates.length) {
                this.duplicates = Arrays.copyOf(this.duplicates, 2 * this.duplicates.length);
            }
            this.duplicates[this.numDuplicates++] = index;
            if (this.secondPass != null) {
                this.secondPass.duplicates.set(index);
            }
        }

        private void addOpticalDuplicate(final long index) {
            if (this.numOpticalDuplicates == this.opticalDuplicates.length) {
                this.opticalDuplicates = Arrays.copyOf(this.opticalDuplicates, 2 * this.opticalDuplicates.length);
            }
            this.opticalDuplicates[this.numOpticalDuplicates++] = index;
            if (this.secondPass != null && this.secondPass.opticalDuplicates != null) {
                this.secondPass.opticalDuplicates.set(index);
            }
        }
    }

    /**
     * Pass 2 on a thread of its own while pass 1 reads: it writes every record as doWork's pass 2 does, into a new
     * file beside OUTPUT that replaces OUTPUT at the end, reading the input again through a ParallelBAMRecordIterator
     * that leaves unvalidated only the records pass 1 validates, and it takes record i only once every duplicate flag
     * below i+1 is decided. Main tells the pair pass, for each batch, below which record every record is read and every
     * mapped primary read without a decided flag is in that batch or an earlier one: a read whose mate is not yet read
     * is not decided until the mate comes, or until the frontier passes the mate's alignment start, after which the
     * read keeps no mate and is not a duplicate. The pair pass decides the rest and raises undecidedFrom.
     *
     * Abandoned, the pass stops and its output is deleted, and doWork finishes as without it: when the pair pass is
     * abandoned, an unpaired end lies below the frontier, a key cannot be packed, a read whose mate was taken to be
     * missing finds it, too many ends or keys wait, PG IDs appear that the output header did not have, or pass 1 or
     * pass 2 fails.
     */
    private final class SecondPassDuringRead implements Runnable {
        // BGZF blocks pass 2's writer keeps in flight while pass 1 reads: deeper than its default, so that the deflaters
        // have work while pass 2 waits for decisions.
        private static final int BLOCKS_IN_FLIGHT_DURING_PASS_ONE = 96;

        // The records pass 1 has validated, for pass 2's reader.
        final ParallelBAMRecordIterator.ValidatedRecords validated = new ParallelBAMRecordIterator.ValidatedRecords();
        // Set by the pair pass thread before undecidedFrom passes them: the duplicates' and optical duplicates'
        // records (opticalDuplicates null when no optical indices are needed), and the unpaired ends marked.
        final IndexBits duplicates = new IndexBits();
        final IndexBits opticalDuplicates;
        long fragmentDuplicates = 0;
        private volatile long undecidedFrom = 0;
        private volatile boolean abandoned = false;
        private volatile boolean consumerWaiting = false;
        private volatile Thread thread = null;

        private final int referenceBits;
        private final long maxPendingFirstMates;
        // Main's: the reads whose mate is not yet read, in record order, with their mate's position(), and which of
        // them have since found their mate or been taken to have none.
        private final IndexBits completedFirstMates = new IndexBits();
        private final IndexBits orphans = new IndexBits();
        private boolean anyOrphans = false;
        private long[] pendingIndex = new long[1024];
        private long[] pendingMate = new long[1024];
        private int pendingHead = 0;
        private int pendingSize = 0;
        private int pgIdsSeenSize = -1;

        // Set by start, before the thread starts.
        private SamHeaderAndIterator input = null;
        private SAMFileHeader outputHeader;
        private Map<String, String> chainedPgIds;
        private Path temp = null;
        private LibraryIdGenerator libraryMetrics;
        // The thread's results: whether it created temp, and whether it wrote every record and closed the output.
        private volatile boolean createdTemp = false;
        // The thread's writer once made, and whether pass 1 has read every record.
        private volatile SAMFileWriter writer = null;
        private volatile boolean passOneDone = false;
        private volatile boolean written = false;
        // The thread's failure, and whether it came after main released the log.
        private volatile Throwable failure = null;
        private volatile boolean failedAfterRelease = false;
        // The thread's log lines held back until main has logged the lines before pass 2's; guarded by heldLines.
        private final List<Object[]> heldLines = new ArrayList<>();
        private boolean logReleased = false;
        private volatile boolean released = false;

        SecondPassDuringRead(final int numberOfReferences, final boolean indexOpticalDuplicates,
                             final long maxPendingFirstMates) {
            this.referenceBits = numberOfReferences <= 1 ? 0 : 32 - Integer.numberOfLeadingZeros(numberOfReferences - 1);
            this.opticalDuplicates = indexOpticalDuplicates ? new IndexBits() : null;
            this.maxPendingFirstMates = maxPendingFirstMates;
        }

        boolean isActive() {
            return !this.abandoned;
        }

        boolean isAbandoned() {
            return this.abandoned;
        }

        boolean isStarted() {
            return this.thread != null;
        }

        /** Whether pgIdsSeen has grown since the output header was made; main's. */
        boolean pgIdsSeenChanged(final int size) {
            return this.pgIdsSeenSize >= 0 && size != this.pgIdsSeenSize;
        }

        /**
         * The library, reference index, 5' coordinate and orientation of an end, packed into a non-negative long, one
         * per distinct four; -1 when they cannot be packed.
         */
        long key(final ReadEndsForMarkDuplicates end) {
            final int libraryBits = 63 - 33 - this.referenceBits;
            if (libraryBits < 0 || end.libraryId < 0 || end.libraryId >= (1L << libraryBits) ||
                    end.read1ReferenceIndex < 0 || (end.read1ReferenceIndex >>> this.referenceBits) != 0 ||
                    (end.orientation != ReadEnds.F && end.orientation != ReadEnds.R)) {
                return -1;
            }
            final long reference = ((long) end.libraryId << this.referenceBits) | end.read1ReferenceIndex;
            return (reference << 33) | (((long) end.read1Coordinate - Integer.MIN_VALUE) << 1) | end.orientation;
        }

        /** Main: the read with this index has a mapped mate not yet read. */
        void firstMate(final long index, final SAMRecord rec) {
            if (this.pendingSize == this.pendingIndex.length) {
                if (this.pendingSize >= this.maxPendingFirstMates) {
                    log.debug("Reads waiting for their mates reached ", this.pendingSize, "; pass 2 beside pass 1 stops.");
                    abandon();
                    return;
                }
                final int n = this.pendingIndex.length;
                final long[] indices = new long[2 * n];
                final long[] mates = new long[2 * n];
                for (int i = 0; i < n; i++) {
                    indices[i] = this.pendingIndex[(this.pendingHead + i) & (n - 1)];
                    mates[i] = this.pendingMate[(this.pendingHead + i) & (n - 1)];
                }
                this.pendingIndex = indices;
                this.pendingMate = mates;
                this.pendingHead = 0;
            }
            final int slot = (this.pendingHead + this.pendingSize) & (this.pendingIndex.length - 1);
            this.pendingIndex[slot] = index;
            this.pendingMate[slot] = position(rec.getMateReferenceIndex(), rec.getMateAlignmentStart());
            this.pendingSize++;
        }

        /** Main: the read with this index, whose mate was not yet read, has found it. */
        void completed(final long firstMateIndex) {
            if (this.anyOrphans && this.orphans.get(firstMateIndex)) {
                // It was taken to have no mate, and so not to be a duplicate.
                abandon();
            } else {
                this.completedFirstMates.set(firstMateIndex);
            }
        }

        /**
         * Main, when a batch goes: the first record that may lack a decided flag, below which every record is read.
         * A read still waiting for a mate whose alignment start lies below the frontier is taken to have none.
         */
        long firstUndecided(final long recordsRead, final long frontier) {
            while (this.pendingSize > 0) {
                final int head = this.pendingHead;
                final long index = this.pendingIndex[head];
                if (!this.completedFirstMates.get(index)) {
                    if (this.pendingMate[head] >= frontier) {
                        return Math.min(recordsRead, index);
                    }
                    this.orphans.set(index);
                    this.anyOrphans = true;
                }
                this.pendingHead = (head + 1) & (this.pendingIndex.length - 1);
                this.pendingSize--;
            }
            return recordsRead;
        }

        /** Main: starts the thread, which writes temp. pgIdsSeenSize is pgIdsSeen's size when chainedPgIds was made. */
        void start(final SamHeaderAndIterator input, final SAMFileHeader outputHeader,
                   final Map<String, String> chainedPgIds, final int pgIdsSeenSize, final Path temp) {
            this.input = input;
            this.outputHeader = outputHeader;
            this.chainedPgIds = chainedPgIds;
            this.pgIdsSeenSize = pgIdsSeenSize;
            this.temp = temp;
            this.libraryMetrics = new LibraryIdGenerator(input.header, flowBasedArguments.FLOW_MODE);
            final Thread pass = new Thread(this, "MarkDuplicates pass 2");
            pass.setDaemon(true);
            pass.start();
            // Until it is set, the thread waits for decisions in 1 ms steps rather than being woken.
            this.thread = pass;
        }

        /** Main: pass 1 has read every record; pass 2's writer goes back to its default number of blocks in flight. */
        void passOneDone() {
            this.passOneDone = true;
            final SAMFileWriter made = this.writer;
            if (made != null) {
                ParallelBAMFileWriter.setBlocksInFlightLimit(made, 0);
            }
        }

        /** The pair pass thread: every duplicate flag below this record is decided. */
        void decided(final long undecidedFrom) {
            if (undecidedFrom < this.undecidedFrom) {
                abandon();
                return;
            }
            this.undecidedFrom = undecidedFrom;
            if (this.consumerWaiting) {
                final Thread pass = this.thread;
                if (pass != null) {
                    LockSupport.unpark(pass);
                }
            }
        }

        /** Stops pass 2; its output is then never used. Any thread. */
        void abandon() {
            this.abandoned = true;
            this.validated.abort();
            final Thread pass = this.thread;
            if (pass != null) {
                LockSupport.unpark(pass);
            }
        }

        /** The input's sort order, which pass 2 logs. */
        SAMFileHeader.SortOrder sortOrder() {
            return this.input.header.getSortOrder();
        }

        /** A line for the log, held back until main releases it. */
        private void logLine(final Object[] parts) {
            synchronized (this.heldLines) {
                if (!this.logReleased) {
                    this.heldLines.add(parts);
                    return;
                }
            }
            log.info(parts);
        }

        /**
         * Main, once it has logged the lines before pass 2's: logs the lines the thread held back, and lets it read on
         * past the records pass 1 decided while reading.
         */
        void releaseLog() {
            synchronized (this.heldLines) {
                for (final Object[] parts : this.heldLines) {
                    log.info(parts);
                }
                this.heldLines.clear();
                this.logReleased = true;
            }
            this.released = true;
            final Thread pass = this.thread;
            if (pass != null) {
                LockSupport.unpark(pass);
            }
        }

        /** The thread: waits for releaseLog; false when abandoned first. */
        private boolean awaitReleased() {
            while (!this.released) {
                if (this.abandoned) {
                    return false;
                }
                LockSupport.parkNanos(this, 1_000_000L);
            }
            return !this.abandoned;
        }

        /**
         * The thread: the first record whose flag is not decided once it is above index; -1 when abandoned. Records
         * decided only once pass 1 has ended wait for releaseLog too, and so do the unplaced reads after them, whose
         * validation may print warnings.
         */
        private long awaitDecided(final long index) {
            while (true) {
                final long decided = this.undecidedFrom;
                if (decided == Long.MAX_VALUE) {
                    return awaitReleased() ? decided : -1;
                }
                if (decided > index) {
                    return decided;
                }
                if (this.abandoned) {
                    return -1;
                }
                this.consumerWaiting = true;
                if (this.undecidedFrom <= index && !this.abandoned) {
                    LockSupport.parkNanos(this, 1_000_000L);
                }
                this.consumerWaiting = false;
            }
        }

        @Override
        public void run() {
            boolean wrote = false;
            final CloseableIterator<SAMRecord> iterator = this.input.iterator;
            try {
                final SAMFileHeader header = this.input.header;
                final SAMFileWriter writer = ParallelBAMFileWriter.makeAsyncWriter(this.outputHeader, true,
                        this.temp.toFile(), true);
                this.createdTemp = true;
                this.temp.toFile().deleteOnExit();
                this.writer = writer;
                if (!this.passOneDone) {
                    ParallelBAMFileWriter.setBlocksInFlightLimit(writer, BLOCKS_IN_FLIGHT_DURING_PASS_ONE);
                    if (this.passOneDone) {
                        ParallelBAMFileWriter.setBlocksInFlightLimit(writer, 0);
                    }
                }
                try (SAMFileWriter out = writer) {
                    final ProgressLoggerInterface progress = new SingleThreadProgressLogger(this::logLine, (int) 1e7,
                            "Written");
                    long recordInFileIndex = 0;
                    long decided = 0;
                    while (iterator.hasNext() && !this.abandoned) {
                        final SAMRecord rec = iterator.next();
                        if (recordInFileIndex >= decided) {
                            decided = awaitDecided(recordInFileIndex);
                            if (decided <= recordInFileIndex) {
                                break;
                            }
                        }

                        // As doWork's pass 2 does for a coordinate-sorted input.
                        final DuplicationMetrics metrics = AbstractMarkDuplicatesCommandLineProgram.addReadToLibraryMetrics(
                                rec, header, this.libraryMetrics, flowBasedArguments.FLOW_MODE);
                        final boolean isDuplicate = this.duplicates.get(recordInFileIndex);
                        if (isDuplicate) {
                            rec.setDuplicateReadFlag(true);
                            metrics.addDuplicateReadToMetrics(rec);
                        } else {
                            rec.setDuplicateReadFlag(false);
                        }
                        final boolean isOpticalDuplicate = this.opticalDuplicates != null &&
                                this.opticalDuplicates.get(recordInFileIndex);
                        if (CLEAR_DT) {
                            rec.setAttribute(DUPLICATE_TYPE_TAG, null);
                        }
                        if (TAGGING_POLICY != DuplicateTaggingPolicy.DontTag && rec.getDuplicateReadFlag()) {
                            if (isOpticalDuplicate) {
                                rec.setAttribute(DUPLICATE_TYPE_TAG, DuplicateType.SEQUENCING.code());
                            } else if (TAGGING_POLICY == DuplicateTaggingPolicy.All) {
                                rec.setAttribute(DUPLICATE_TYPE_TAG, DuplicateType.LIBRARY.code());
                            }
                        }
                        recordInFileIndex++;
                        if (REMOVE_DUPLICATES && rec.getDuplicateReadFlag()) {
                            continue;
                        }
                        if (REMOVE_SEQUENCING_DUPLICATES && isOpticalDuplicate) {
                            continue;
                        }
                        if (PROGRAM_RECORD_ID != null && pgTagArgumentCollection.ADD_PG_TAG_TO_READS) {
                            rec.setAttribute(SAMTag.PG.name(), this.chainedPgIds.get(rec.getStringAttribute(SAMTag.PG.name())));
                        }
                        out.addAlignment(rec);
                        progress.record(rec);
                    }
                    if (awaitReleased()) {
                        log.info("Writing complete. Closing input iterator.");
                        iterator.close();
                        log.info("Duplicate Index cleanup.");
                        log.info("Getting Memory Stats.");
                        reportMemoryStats("Before output close");
                        wrote = true;
                    }
                }
            } catch (final Throwable t) {
                wrote = false;
                this.failedAfterRelease = this.released;
                this.failure = t;
            } finally {
                try {
                    iterator.close();
                } catch (final RuntimeException e) {
                    wrote = false;
                }
                this.written = wrote && !this.abandoned;
            }
        }

        /** Main: whether the thread wrote every record and closed the output, waiting for it. */
        boolean awaitWritten() {
            join();
            return this.written && !this.abandoned;
        }

        /**
         * Main, after awaitWritten: the RuntimeException the thread failed with after releaseLog, unabandoned, which
         * pass 2 after pass 1 would meet at the same record; null for none or any other failure.
         */
        RuntimeException failureAfterRelease() {
            return this.failure instanceof RuntimeException && this.failedAfterRelease && !this.abandoned &&
                    this.createdTemp ? (RuntimeException) this.failure : null;
        }

        private void join() {
            final Thread pass = this.thread;
            if (pass == null) {
                return;
            }
            boolean interrupted = false;
            while (true) {
                try {
                    pass.join();
                    break;
                } catch (final InterruptedException e) {
                    interrupted = true;
                }
            }
            if (interrupted) {
                Thread.currentThread().interrupt();
            }
        }

        /** Main: moves the output written into place; false, leaving it, when that fails. */
        boolean moveTo(final File output) {
            final Path target = output.toPath();
            try {
                try {
                    Files.move(this.temp, target, StandardCopyOption.ATOMIC_MOVE);
                } catch (final AtomicMoveNotSupportedException e) {
                    Files.move(this.temp, target, StandardCopyOption.REPLACE_EXISTING);
                }
                this.temp = null;
                return true;
            } catch (final IOException | RuntimeException e) {
                return false;
            }
        }

        /** Main: puts the library metrics pass 2 counted in place of target's, which pass 2 has not touched. */
        void addMetricsTo(final LibraryIdGenerator target) {
            for (final Map.Entry<String, DuplicationMetrics> entry : this.libraryMetrics.getMetricsByLibraryMap().entrySet()) {
                target.addMetricsByLibrary(entry.getKey(), entry.getValue());
            }
        }

        /** Main: stops pass 2, waits for it and deletes its output. */
        void discard() {
            abandon();
            if (this.thread != null) {
                join();
            } else if (this.input != null) {
                CloserUtil.close(this.input.iterator);
            }
            if (this.temp != null && this.createdTemp) {
                try {
                    Files.deleteIfExists(this.temp);
                } catch (final IOException | RuntimeException e) {
                    // Left behind; deleteOnExit may still remove it.
                }
            }
            this.temp = null;
        }
    }

    /**
     * Starts pass 2 beside pass 1 once recordsRead records are read, writing a new file beside OUTPUT; abandons it if
     * that fails, which pass 2 after pass 1 then meets in its turn.
     */
    private void startSecondPass(final SecondPassDuringRead secondPass, final long recordsRead) {
        secondPass.validated.advance(recordsRead);
        Path temp = null;
        SamHeaderAndIterator input = null;
        try {
            temp = siblingName(OUTPUT);
            input = openInputs(false, false, decodeThreads(), RawBlockBAMRecord.FACTORY, null, secondPass.validated);
            if (input == null || !(input.iterator instanceof ParallelBAMRecordIterator)) {
                throw new IllegalStateException("not a ParallelBAMRecordIterator");
            }
            final SAMFileHeader outputHeader = input.header.clone();
            COMMENT.forEach(outputHeader::addComment);
            final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);
            secondPass.start(input, outputHeader, chainedPgIds, pgIdsSeen.size(), temp);
            log.debug("Pass 2 started beside pass 1 after ", recordsRead, " records.");
        } catch (final IOException | RuntimeException e) {
            secondPass.abandon();
            if (input != null && !secondPass.isStarted()) {
                CloserUtil.close(input.iterator);
            }
        }
    }

    /** A hidden name beside output that no file has; pass 2's writer creates the file, failing if one appeared since. */
    private static Path siblingName(final File output) throws IOException {
        final Path target = output.toPath().toAbsolutePath();
        final Path directory = target.getParent();
        for (int attempt = 0; ; attempt++) {
            final Path candidate = directory.resolve("." + target.getFileName() + "." +
                    Long.toHexString(ThreadLocalRandom.current().nextLong()) + ".tmp");
            if (Files.notExists(candidate, LinkOption.NOFOLLOW_LINKS)) {
                return candidate;
            }
            if (attempt >= 16) {
                throw new FileAlreadyExistsException(candidate.toString());
            }
        }
    }

    /**
     * Comparator for ReadEndsForMarkDuplicates that orders by read1 position then pair orientation then read2 position.
     */
    static class ReadEndsMDComparator implements Comparator<ReadEndsForMarkDuplicates> {

        final boolean useBarcodes;

        public ReadEndsMDComparator(final boolean useBarcodes) {
            this.useBarcodes = useBarcodes;
        }

        public int compare(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs) {
            int compareDifference = lhs.libraryId - rhs.libraryId;
            if (useBarcodes) {
                final ReadEndsForMarkDuplicatesWithBarcodes lhsWithBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) lhs;
                final ReadEndsForMarkDuplicatesWithBarcodes rhsWithBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) rhs;
                if (compareDifference == 0) {
                    compareDifference = Integer.compare(lhsWithBarcodes.barcode, rhsWithBarcodes.barcode);
                }
                if (compareDifference == 0) {
                    compareDifference = Integer.compare(lhsWithBarcodes.readOneBarcode, rhsWithBarcodes.readOneBarcode);
                }
                if (compareDifference == 0) {
                    compareDifference = Integer.compare(lhsWithBarcodes.readTwoBarcode, rhsWithBarcodes.readTwoBarcode);
                }
            }
            if (compareDifference == 0) {
                compareDifference = lhs.read1ReferenceIndex - rhs.read1ReferenceIndex;
            }
            if (compareDifference == 0) {
                compareDifference = lhs.read1Coordinate - rhs.read1Coordinate;
            }
            if (compareDifference == 0) {
                compareDifference = lhs.orientation - rhs.orientation;
            }
            if (compareDifference == 0) {
                compareDifference = lhs.read2ReferenceIndex - rhs.read2ReferenceIndex;
            }
            if (compareDifference == 0) {
                compareDifference = lhs.read2Coordinate - rhs.read2Coordinate;
            }

            if (compareDifference == 0) {
                compareDifference = lhs.getTile() - rhs.getTile();
            }

            if (compareDifference == 0) {
                compareDifference = lhs.getX() - rhs.getX();
            }

            if (compareDifference == 0) {
                compareDifference = lhs.getY() - rhs.getY();
            }

            // The following is arbitrary and is only included for completeness.
            // Other implementations may chose to forgo this tiebreak if they do not have
            // access to the index-in-file of the records (e.g. SPARK implmentations)

            if (compareDifference == 0) {
                compareDifference = Long.compare(lhs.read1IndexInFile, rhs.read1IndexInFile);
            }
            if (compareDifference == 0) {
                compareDifference = Long.compare(lhs.read2IndexInFile, rhs.read2IndexInFile);
            }

            return compareDifference;
        }
    }
}
