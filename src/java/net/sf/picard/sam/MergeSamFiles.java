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
package net.sf.picard.sam;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;

/**
 * Reads a SAM or BAM file and combines the output to one file
 *
 * @author Dave Tefft
 */
public class MergeSamFiles extends CommandLineProgram {
    private static final Log log = Log.getInstance(MergeSamFiles.class);

    // Usage and parameters
    @Usage(programVersion="1.0")
    public String USAGE = "Merges multiple SAM/BAM files into one file.\n";

    @Option(shortName="I", doc="SAM or BAM input file", minElements=1)
    public List<File> INPUT = new ArrayList<File>();

    @Option(shortName="O", doc="SAM or BAM file to write merged result to")
    public File OUTPUT;

    @Option(shortName=StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, doc="Sort order of output file", optional=true)
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    @Option(doc="If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise.",
    shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = false;

    @Option(shortName="MSD", doc="Merge the seqeunce dictionaries", optional=true)
    public boolean MERGE_SEQUENCE_DICTIONARIES = false;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new MergeSamFiles().instanceMain(argv));
    }

    /** Combines multiple SAM/BAM files into one. */
    @Override
	protected int doWork() {
        boolean matchedSortOrders = true;

        // Open the files for reading and writing
        final List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        {
            SAMSequenceDictionary dict = null; // Used to try and reduce redundant SDs in memory

            for (final File inFile : INPUT) {
                IoUtil.assertFileIsReadable(inFile);
                final SAMFileReader in = new SAMFileReader(inFile);
                readers.add(in);

                // A slightly hackish attempt to keep memory consumption down when merging multiple files with
                // large sequence dictionaries (10,000s of sequences). If the dictionaries are identical, then
                // replace the duplicate copies with a single dictionary to reduce the memory footprint. 
                if (dict == null) {
                    dict = in.getFileHeader().getSequenceDictionary();
                }
                else if (dict.equals(in.getFileHeader().getSequenceDictionary())) {
                    in.getFileHeader().setSequenceDictionary(dict);
                }

                matchedSortOrders = matchedSortOrders && in.getFileHeader().getSortOrder() == SORT_ORDER;
            }
        }

        // If all the input sort orders match the output sort order then just merge them and
        // write on the fly, otherwise setup to merge and sort before writing out the final file
        IoUtil.assertFileIsWritable(OUTPUT);
        final MergingSamRecordIterator iterator;
        final SAMFileWriter out;

        if (matchedSortOrders || SORT_ORDER == SAMFileHeader.SortOrder.unsorted || ASSUME_SORTED) {
            log.info("Input files are in same order as output so sorting to temp directory is not needed.");
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SORT_ORDER, MERGE_SEQUENCE_DICTIONARIES);
            iterator = new MergingSamRecordIterator(headerMerger, ASSUME_SORTED);
            out = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerMerger.getMergedHeader(), true, OUTPUT);
        }
        else {
            log.info("Sorting input files using temp directory " + TMP_DIR);
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers,
                                                                             SAMFileHeader.SortOrder.unsorted,
                                                                             MERGE_SEQUENCE_DICTIONARIES);
            iterator = new MergingSamRecordIterator(headerMerger, false);
            final SAMFileHeader header = headerMerger.getMergedHeader();
            header.setSortOrder(SORT_ORDER);
            out = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);
        }

        // Lastly loop through and write out the records
        for (long numRecords = 1; iterator.hasNext(); ++numRecords) {
            final SAMRecord record = iterator.next();
            out.addAlignment(record);
            if (numRecords % 10000000 == 0) {
                log.info(numRecords + " records read.");
            }
        }
        log.info("Finished reading inputs.");
        out.close();
        return 0;
    }

}