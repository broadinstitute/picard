/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package picard.sam;

import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SamRecordIntervalIteratorFactory;
import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Reads a SAM or BAM file and combines the output to one file
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = MergeSamFiles.USAGE_SUMMARY + MergeSamFiles.USAGE_DETAILS,
        oneLineSummary = MergeSamFiles.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class MergeSamFiles extends CommandLineProgram {
    private static final Log log = Log.getInstance(MergeSamFiles.class);

    static final String USAGE_SUMMARY = "Merges multiple SAM and/or BAM files into a single file.  ";
    static final String USAGE_DETAILS = "This tool is used for combining SAM and/or BAM files from different runs or read groups, similarly " +
            "to the \"merge\" function of Samtools (http://www.htslib.org/doc/samtools.html).  " +
            "<br /><br />Note that to prevent errors in downstream processing, it is critical to identify/label read groups appropriately. " +
            "If different samples contain identical read group IDs, this tool will avoid collisions by modifying the read group IDs to be " +
            "unique. For more information about read groups, see the " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=6472'>GATK Dictionary entry.</a> <br /><br />" +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar MergeSamFiles \\<br />" +
            "      I=input_1.bam \\<br />" +
            "      I=input_2.bam \\<br />" +
            "      O=merged_files.bam" +
            "</pre>" +
            "<hr />"
           ;
    @Argument(shortName = "I", doc = "SAM or BAM input file", minElements = 1)
    public List<File> INPUT = new ArrayList<File>();

    @Argument(shortName = "O", doc = "SAM or BAM file to write merged result to")
    public File OUTPUT;

    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, doc = "Sort order of output file", optional = true)
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    @Argument(doc = "If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise.",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = false;

    @Argument(shortName = "MSD", doc = "Merge the sequence dictionaries", optional = true)
    public boolean MERGE_SEQUENCE_DICTIONARIES = false;

    @Argument(doc = "Option to create a background thread to encode, " +
            "compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases " +
            "runtime by ~20% when writing out a compressed BAM file.")
    public boolean USE_THREADING = false;

    @Argument(doc = "Comment(s) to include in the merged output file's header.", optional = true, shortName = "CO")
    public List<String> COMMENT = new ArrayList<String>();

    @Argument(shortName = "RGN", doc = "An interval list file that contains the locations of the positions to merge. "+
            "Assume bam are sorted and indexed. "+
            "The resulting file will contain alignments that may overlap with genomic regions outside the requested region. "+
            "Unmapped reads are discarded.",
            optional = true)
    public File INTERVALS = null;

    private static final int PROGRESS_INTERVAL = 1000000;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new MergeSamFiles().instanceMain(argv));
    }

    /** Combines multiple SAM/BAM files into one. */
    @Override
    protected int doWork() {
        boolean matchedSortOrders = true;
        
        // read interval list if it is defined
        final List<Interval> intervalList = (INTERVALS == null ? null : IntervalList.fromFile(INTERVALS).uniqued().getIntervals() );
        // map reader->iterator used if INTERVALS is defined
        final Map<SamReader, CloseableIterator<SAMRecord> > samReaderToIterator = new HashMap<SamReader, CloseableIterator<SAMRecord> >(INPUT.size());
        
        // Open the files for reading and writing
        final List<SamReader> readers = new ArrayList<SamReader>();
        final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
        {
            SAMSequenceDictionary dict = null; // Used to try and reduce redundant SDs in memory

            for (final File inFile : INPUT) {
                IOUtil.assertFileIsReadable(inFile);
                final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(inFile);
                 if ( INTERVALS != null ) {
                     if( ! in.hasIndex() ) throw new PicardException("Merging with interval but Bam file is not indexed "+ inFile);
                     final CloseableIterator<SAMRecord> samIterator = new SamRecordIntervalIteratorFactory().makeSamRecordIntervalIterator(in, intervalList, true);
                     samReaderToIterator.put(in, samIterator);
                 }

                readers.add(in);
                headers.add(in.getFileHeader());

                // A slightly hackish attempt to keep memory consumption down when merging multiple files with
                // large sequence dictionaries (10,000s of sequences). If the dictionaries are identical, then
                // replace the duplicate copies with a single dictionary to reduce the memory footprint. 
                if (dict == null) {
                    dict = in.getFileHeader().getSequenceDictionary();
                } else if (dict.equals(in.getFileHeader().getSequenceDictionary())) {
                    in.getFileHeader().setSequenceDictionary(dict);
                }

                matchedSortOrders = matchedSortOrders && in.getFileHeader().getSortOrder() == SORT_ORDER;
            }
        }

        // If all the input sort orders match the output sort order then just merge them and
        // write on the fly, otherwise setup to merge and sort before writing out the final file
        IOUtil.assertFileIsWritable(OUTPUT);
        final boolean presorted;
        final SAMFileHeader.SortOrder headerMergerSortOrder;
        final boolean mergingSamRecordIteratorAssumeSorted;

        if (matchedSortOrders || SORT_ORDER == SAMFileHeader.SortOrder.unsorted || ASSUME_SORTED || INTERVALS != null ) {
            log.info("Input files are in same order as output so sorting to temp directory is not needed.");
            headerMergerSortOrder = SORT_ORDER;
            mergingSamRecordIteratorAssumeSorted = ASSUME_SORTED;
            presorted = true;
        } else {
            log.info("Sorting input files using temp directory " + TMP_DIR);
            headerMergerSortOrder = SAMFileHeader.SortOrder.unsorted;
            mergingSamRecordIteratorAssumeSorted = false;
            presorted = false;
        }
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(headerMergerSortOrder, headers, MERGE_SEQUENCE_DICTIONARIES);
        final MergingSamRecordIterator iterator;
        // no interval defined, get an iterator for the whole bam
        if( intervalList == null) {
            iterator = new MergingSamRecordIterator(headerMerger, readers, mergingSamRecordIteratorAssumeSorted);
        }
        else {
            // show warning related to https://github.com/broadinstitute/picard/pull/314/files
            log.info("Warning: merged bams from different interval lists may contain the same read in both files");
            iterator = new MergingSamRecordIterator(headerMerger, samReaderToIterator, true);
        }
        final SAMFileHeader header = headerMerger.getMergedHeader();
        for (final String comment : COMMENT) {
            header.addComment(comment);
        }
        header.setSortOrder(SORT_ORDER);
        final SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
        if (USE_THREADING) {
            samFileWriterFactory.setUseAsyncIo(true);
        }
        final SAMFileWriter out = samFileWriterFactory.makeSAMOrBAMWriter(header, presorted, OUTPUT);

        // Lastly loop through and write out the records
        final ProgressLogger progress = new ProgressLogger(log, PROGRESS_INTERVAL);
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            out.addAlignment(record);
            progress.record(record);
        }

        log.info("Finished reading inputs.");
        for(final CloseableIterator<SAMRecord> iter : samReaderToIterator.values())  CloserUtil.close(iter);
        CloserUtil.close(readers);
        out.close();
        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (CREATE_INDEX && SORT_ORDER != SAMFileHeader.SortOrder.coordinate) {
            return new String[]{"Can't CREATE_INDEX unless SORT_ORDER is coordinate"};
        }
        return null;
    }

}
