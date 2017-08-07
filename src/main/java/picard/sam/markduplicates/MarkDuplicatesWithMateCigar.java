/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.*;
import picard.cmdline.programgroups.SamOrBam;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;


import java.util.*;

/**
 * An even better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments.
 * <p/>
 * This tool differs with MarkDuplicates as it may break ties differently.  Furthermore,
 * as it is a one-pass algorithm, it cannot know the program records contained in the file
 * that should be chained in advance.  Therefore it will only be able to examine the header
 * to attempt to infer those program group records that have no associated previous program
 * group record. If a read is encountered without a program record, or not one as previously
 * defined, it will not be updated.
 * <p/>
 * This tool will also not work with alignments that have large gaps or skips, such as those
 * from RNA-seq data.  This is due to the need to buffer small genomic windows to ensure
 * integrity of the duplicate marking, while large skips (ex. skipping introns) in the
 * alignment records would force making that window very large, thus exhausting memory.
 *
 * @author Nils Homer
 */
@CommandLineProgramProperties(
        summary = MarkDuplicatesWithMateCigar.USAGE_SUMMARY + MarkDuplicatesWithMateCigar.USAGE_DETAILS,
        oneLineSummary =  MarkDuplicatesWithMateCigar.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class MarkDuplicatesWithMateCigar extends AbstractMarkDuplicatesCommandLineProgram {
    static final String USAGE_SUMMARY = "Identifies duplicate reads, accounting for mate CIGAR.  ";
    static final String USAGE_DETAILS = "This tool locates and tags duplicate reads (both PCR and optical) in a BAM or SAM file, where " +
            "duplicate reads are defined as originating from the same original fragment of DNA, taking into account the CIGAR string of " +
            "read mates. <br /><br />" +
            "" +
            "It is intended as an improvement upon the original MarkDuplicates algorithm, from which it differs in several ways, including" +
            "differences in how it breaks ties. It may be the most effective duplicate marking program available, as it handles all cases " +
            "including clipped and gapped alignments and locates duplicate molecules using mate cigar information. However, please note " +
            "that it is not yet used in the Broad's production pipeline, so use it at your own risk. <br /><br />" +
            "" +
            "Note also that this tool will not work with alignments that have large gaps or deletions, such as those from RNA-seq data.  " +
            "This is due to the need to buffer small genomic windows to ensure integrity of the duplicate marking, while large skips " +
            "(ex. skipping introns) in the alignment records would force making that window very large, thus exhausting memory. <br />" +
            "" +
            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar MarkDuplicatesWithMateCigar \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=mark_dups_w_mate_cig.bam \\<br />" +
            "      M=mark_dups_w_mate_cig_metrics.txt" +
            "</pre>" +
            "<hr />";

    private final Log log = Log.getInstance(MarkDuplicatesWithMateCigar.class);

    @Argument(doc = "The minimum distance to buffer records to account for clipping on the 5' end of the records. " +
            "For a given alignment, this parameter controls the width of the window to search for duplicates of that alignment. " +
            "Due to 5' read clipping, duplicates do not necessarily have the same 5' alignment coordinates, so the algorithm " +
            "needs to search around the neighborhood. For single end sequencing data, the neighborhood is only determined by " +
            "the amount of clipping (assuming no split reads), thus setting MINIMUM_DISTANCE to twice the sequencing read length " +
            "should be sufficient. For paired end sequencing, the neighborhood is also determined by the fragment insert size, " +
            "so you may want to set MINIMUM_DISTANCE to something like twice the 99.5% percentile of the fragment insert size " +
            "distribution (see CollectInsertSizeMetrics). Or you can set this number to -1 to use either a) twice the first read's " +
            "read length, or b) 100, whichever is smaller. Note that the larger the window, the greater the RAM requirements, so " +
            "you could run into performance limitations if you use a value that is unnecessarily large.", optional = true)
    public int MINIMUM_DISTANCE = -1;

    @Argument(doc = "Skip record pairs with no mate cigar and include them in the output.")
    boolean SKIP_PAIRS_WITH_NO_MATE_CIGAR = true;

    @Argument(doc = "The block size for use in the coordinate-sorted record buffer.", optional = true)
    public int BLOCK_SIZE = 100000;

    /** Warnings that will only be emitted once */
    private boolean warnedNullProgramRecords = false;
    private boolean warnedMissingProgramRecords = false;

    /** Stock main method. */
    public static void main(final String[] args) {
        new MarkDuplicatesWithMateCigar().instanceMainWithExit(args);
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

        // Since this is one-pass, unlike MarkDuplicates, we cannot only chain together program
        // group records we have seen, we have to assume all of them may be seen.  We can perhaps
        // filter out any program groups which have been referenced previously.
        setPGIdsSeen(outputHeader);
        // Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
        final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);

        // Open the output
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader,
                true,
                OUTPUT);

        // Create the mark duplicate iterator.  The duplicate marking is handled by the iterator, conveniently.
        final MarkDuplicatesWithMateCigarIterator iterator = new MarkDuplicatesWithMateCigarIterator(headerAndIterator.header,
                headerAndIterator.iterator,
                this.opticalDuplicateFinder,
                this.DUPLICATE_SCORING_STRATEGY,
                this.MINIMUM_DISTANCE,
                this.REMOVE_DUPLICATES,
                this.SKIP_PAIRS_WITH_NO_MATE_CIGAR,
                this.MAX_RECORDS_IN_RAM,
                this.BLOCK_SIZE,
                this.TMP_DIR);

        // progress logger!
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");

        // Go through the records
        for (final SAMRecord record : new IterableAdapter<SAMRecord>(iterator)) {
            if (progress.record(record)) {
                iterator.logMemoryStats(log);
            }

            // Update the program record if necessary
            updateProgramRecord(record, chainedPgIds);

            // Write the alignment
            out.addAlignment(record);
        }

        // remember to close the inputs
        iterator.close();

        out.close();

        // For convenience to reference
        final Histogram<Short> opticalDupesByLibraryId = iterator.getOpticalDupesByLibraryId();

        // Log info
        log.info("Processed " + progress.getCount() + " records");
        log.info("Found " + iterator.getNumRecordsWithNoMateCigar() + " records with no mate cigar optional tag.");
        log.info("Marking " + iterator.getNumDuplicates() + " records as duplicates.");
        log.info("Found " + ((long) opticalDupesByLibraryId.getSumOfValues()) + " optical duplicate clusters."); // cast as long due to returning a double

        // Write out the metrics
        finalizeAndWriteMetrics(iterator.getLibraryIdGenerator());

        return 0;
    }

    /**
     * Updates the program record if necessary.
     */
    private void updateProgramRecord(final SAMRecord record, final Map<String, String> chainedPgIds) {
        if (PROGRAM_RECORD_ID != null) {
            final String pgId = record.getStringAttribute(SAMTag.PG.name());
            if (null == pgId) {
                if (!warnedNullProgramRecords) {
                    warnedNullProgramRecords = true;
                    log.warn("Encountered a record with no program record, program group chaining will not occur for this read: " + record);
                } // else already warned!
            } else if (!chainedPgIds.containsKey(pgId)) {
                if (!warnedMissingProgramRecords) {
                    warnedMissingProgramRecords = true;
                    log.warn("Encountered a record with an intermediate program record, program group chaining will not occur for this read: " + record);
                } // else already warned!
            } else {
                record.setAttribute(SAMTag.PG.name(), chainedPgIds.get(pgId));
            }
        }
    }

    /**
     * Generate the list of program records seen in the SAM file, approximating this with those in the header that were not
     * themselves mentioned elsewhere.
     */
    private void setPGIdsSeen(final SAMFileHeader header) {
        final Set<String> pgIdsSeenAsPrevious = new HashSet<String>();

        // get all program record ids that are mentioned as previously seen
        for (final SAMProgramRecord samProgramRecord : header.getProgramRecords()) {
            final String previousProgramGroupID = samProgramRecord.getPreviousProgramGroupId();
            if (null != previousProgramGroupID) pgIdsSeenAsPrevious.add(previousProgramGroupID);
        }

        // ignore those that were previously seen
        for (final SAMProgramRecord samProgramRecord : header.getProgramRecords()) {
            final String pgId = samProgramRecord.getId();
            if (!pgIdsSeenAsPrevious.contains(pgId)) this.pgIdsSeen.add(pgId);
        }
    }
}