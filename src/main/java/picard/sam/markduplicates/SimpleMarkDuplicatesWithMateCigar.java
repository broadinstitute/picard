/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordDuplicateComparator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.*;
import picard.PicardException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.Testing;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;
import picard.sam.markduplicates.util.LibraryIdGenerator;
import picard.sam.markduplicates.util.ReadEnds;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * This is a simple tool to mark duplicates using the DuplicateSetIterator, DuplicateSet, and SAMRecordDuplicateComparator.
 * 
 * Users should continue to use MarkDuplicates in general.  The main motivation of this tool was the fact that 
 * MarkDuplicates has many, many, many useful test cases, but few unit tests for validating individual duplicate sets. To
 * test the DuplicateSetIterator, DuplicateSet, and SAMRecordDuplicateComparator, the most expedient method was to write
 * this tool and make sure it behaves similarly to MarkDuplicates.  Not the best, I know, but good enough.  NH 06/25/2015.
 *  
 * 
 * See MarkDuplicates for more details.
 *
 * @author nhomer
 */
@CommandLineProgramProperties(
        summary = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. " +
                "All records are then written to the output file with the duplicate records flagged.",
        oneLineSummary = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules.",
        programGroup = Testing.class
)
public class SimpleMarkDuplicatesWithMateCigar extends MarkDuplicates {
    private final Log log = Log.getInstance(MarkDuplicatesWithMateCigar.class);

    /** Stock main method. */
    public static void main(final String[] args) {
        new MarkDuplicatesWithMateCigar().instanceMainWithExit(args);
    }

    private class ReadEndsForSimpleMarkDuplicatesWithMateCigar extends ReadEnds {
    }
    
    private static boolean isPairedAndBothMapped(final SAMRecord record) {
        return record.getReadPairedFlag() &&
                !record.getReadUnmappedFlag() &&
                !record.getMateUnmappedFlag();
        
    }

    /**
     * Main work method.
     */
    protected int doWork() {
        IOUtil.assertInputsAreValid(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(METRICS_FILE);

        // Open the inputs
        final SamHeaderAndIterator headerAndIterator = openInputs();
        final SAMFileHeader header = headerAndIterator.header;

        // Create the output header
        final SAMFileHeader outputHeader = header.clone();

        if (outputHeader.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException("This program requires inputs in coordinate SortOrder");
        }

        COMMENT.forEach(outputHeader::addComment);
        
        // Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
        final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);

        // Open the output
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader,
                false,
                OUTPUT);

        final SAMRecordDuplicateComparator comparator = new SAMRecordDuplicateComparator(Collections.singletonList(headerAndIterator.header));
        comparator.setScoringStrategy(this.DUPLICATE_SCORING_STRATEGY);

        final CloseableIterator<DuplicateSet> iterator = getDuplicateSetIterator(headerAndIterator, comparator);

        // progress logger!
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");

        int numDuplicates = 0;
        
        libraryIdGenerator = new LibraryIdGenerator(headerAndIterator.header);
        
        for (final DuplicateSet duplicateSet : new IterableAdapter<>(iterator)) {
            final SAMRecord representative = duplicateSet.getRepresentative();
            final boolean doOpticalDuplicateTracking = (this.READ_NAME_REGEX != null) &&
                    isPairedAndBothMapped(representative) &&
                    representative.getFirstOfPairFlag();
            final Set<String> duplicateReadEndsSeen = new HashSet<>();
            
            final List<ReadEnds> duplicateReadEnds = new ArrayList<>();
            for (final SAMRecord record : duplicateSet.getRecords()) {

                // get the metrics for the library of this read (creating a new one if needed)
                final String library = LibraryIdGenerator.getLibraryName(header, record);
                DuplicationMetrics metrics = libraryIdGenerator.getMetricsByLibrary(library);
                if (metrics == null) {
                    metrics = new DuplicationMetrics();
                    metrics.LIBRARY = library;
                    libraryIdGenerator.addMetricsByLibrary(library, metrics);
                }

                if (record.isSecondaryOrSupplementary()) {
                    ++metrics.SECONDARY_OR_SUPPLEMENTARY_RDS;
                } else {

                    // First bring the simple metrics up to date
                    if (record.getReadUnmappedFlag()) {
                        ++metrics.UNMAPPED_READS;
                    } else if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                        ++metrics.UNPAIRED_READS_EXAMINED;
                    } else {
                        ++metrics.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
                    }

                    if (record.getDuplicateReadFlag()) {
                        // Update the duplication metrics
                        if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                            ++metrics.UNPAIRED_READ_DUPLICATES;
                        } else {
                            ++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
                        }
                        numDuplicates++;
                    }
                    
                    // To track optical duplicates, store a set of locations for mapped pairs, first end only.  We care about orientation relative
                    // to the first end of the pair for optical duplicate tracking, which is more stringent than PCR duplicate tracking.
                    if (doOpticalDuplicateTracking &&
                            isPairedAndBothMapped(record) &&
                            !duplicateReadEndsSeen.contains(record.getReadName())) {
                        
                        final ReadEndsForSimpleMarkDuplicatesWithMateCigar readEnd = new ReadEndsForSimpleMarkDuplicatesWithMateCigar();
                        // set orientation for optical duplicates
                        if (record.getFirstOfPairFlag()) {
                            readEnd.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(record.getReadNegativeStrandFlag(), record.getMateNegativeStrandFlag());
                        } else {
                            readEnd.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(record.getMateNegativeStrandFlag(), record.getReadNegativeStrandFlag());
                        }
                        if (opticalDuplicateFinder.addLocationInformation(record.getReadName(), readEnd)) {
                            if (null != record.getReadGroup()) {
                                final short index = libraryIdGenerator.getLibraryId(record);
                                readEnd.setLibraryId(index);
                            }
                        }
                        duplicateReadEnds.add(readEnd);
                        duplicateReadEndsSeen.add(record.getReadName());
                    }
                }

                if (!this.REMOVE_DUPLICATES || !record.getDuplicateReadFlag()) {
                    if (PROGRAM_RECORD_ID != null) {
                        record.setAttribute(SAMTag.PG.name(), chainedPgIds.get(record.getStringAttribute(SAMTag.PG.name())));
                    }
                    out.addAlignment(record);
                    progress.record(record);
                }
            }

            // Track the optical duplicates
            if (this.READ_NAME_REGEX != null && 1 < duplicateReadEnds.size()) {
                AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(duplicateReadEnds, duplicateReadEnds.get(0), opticalDuplicateFinder, libraryIdGenerator);
            }
        }

        // remember to close the inputs
        iterator.close();

        out.close();

        if (this.READ_NAME_REGEX == null) {
            log.warn("Skipped optical duplicate cluster discovery; library size estimation may be inaccurate!");
        } else {
            log.info("Found " + (this.libraryIdGenerator.getNumberOfOpticalDuplicateClusters()) + " optical duplicate clusters.");
        }

        // Log info
        log.info("Processed " + progress.getCount() + " records");
        log.info("Marking " + numDuplicates + " records as duplicates.");

        // Write out the metrics
        finalizeAndWriteMetrics(libraryIdGenerator);

        return 0;
    }

    protected CloseableIterator<DuplicateSet> getDuplicateSetIterator(final SamHeaderAndIterator headerAndIterator, final SAMRecordDuplicateComparator comparator) {
        return new DuplicateSetIterator(headerAndIterator.iterator,
                    headerAndIterator.header,
                    false,
                    comparator);
    }
}
