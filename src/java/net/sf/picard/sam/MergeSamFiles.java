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
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.PicardException;
import net.sf.samtools.*;

/**
 * Reads a SAM or BAM file and combines the output to one file
 *
 * @author Tim Fennell
 */
public class MergeSamFiles extends CommandLineProgram {
    private static final Log log = Log.getInstance(MergeSamFiles.class);

    // Usage and parameters
    @Usage
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

    @Option(doc="Option to enable a simple two-thread producer consumer version of the merge algorithm that " +
            "uses one thread to read and merge the records from the input files and another thread to encode, " +
            "compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases " +
            "runtime by ~20% when writing out a compressed BAM file.")
    public boolean USE_THREADING = false;

    @Option(doc="Comment(s) to include in the merged output file's header.", optional=true, shortName="CO")
    public List<String> COMMENT = new ArrayList<String>();

    private static final int PROGRESS_INTERVAL = 1000000;

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
        final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
        {
            SAMSequenceDictionary dict = null; // Used to try and reduce redundant SDs in memory

            for (final File inFile : INPUT) {
                IoUtil.assertFileIsReadable(inFile);
                final SAMFileReader in = new SAMFileReader(inFile);
                readers.add(in);
                headers.add(in.getFileHeader());

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
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SORT_ORDER, headers, MERGE_SEQUENCE_DICTIONARIES);
            iterator = new MergingSamRecordIterator(headerMerger, readers, ASSUME_SORTED);
            out = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerMerger.getMergedHeader(), true, OUTPUT);
        }
        else {
            log.info("Sorting input files using temp directory " + TMP_DIR);
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.unsorted,
                                                                             headers,
                                                                             MERGE_SEQUENCE_DICTIONARIES);
            iterator = new MergingSamRecordIterator(headerMerger, readers, false);
            final SAMFileHeader header = headerMerger.getMergedHeader();
            header.setSortOrder(SORT_ORDER);
            for (String comment : COMMENT) {
                header.addComment(comment);
            }
            out = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);
        }

        // Lastly loop through and write out the records
        if (USE_THREADING) {
            final BlockingQueue<SAMRecord> queue = new ArrayBlockingQueue<SAMRecord>(10000);
            final AtomicBoolean producerSuccceeded = new AtomicBoolean(false);
            final AtomicBoolean consumerSuccceeded = new AtomicBoolean(false);
            Runnable producer = new Runnable() {
                public void run() {
                    try {
                        while (iterator.hasNext()) {
                            queue.put(iterator.next());
                        }
                        producerSuccceeded.set(true);
                    }
                    catch (InterruptedException ie) {
                        throw new PicardException("Interrupted reading SAMRecord to merge.", ie);
                    }
                }
            };

            Runnable consumer = new Runnable() {
                public void run() {
                    try {
                        long i = 0;
                        SAMRecord rec = null;

                        while ((rec = queue.poll(15, TimeUnit.SECONDS)) != null) {
                            out.addAlignment(rec);
                            if (++i % PROGRESS_INTERVAL == 0) log.info(i + " records processed.");
                        }
                        consumerSuccceeded.set(true);
                    }
                    catch (InterruptedException ie) {
                        throw new PicardException("Interrupted writing SAMRecord to output file.", ie);
                    }
                }
            };

            Thread producerThread = new Thread(producer);
            Thread consumerThread = new Thread(consumer);
            producerThread.start();
            consumerThread.start();

            try {
                consumerThread.join();
                producerThread.join();
            }
            catch (InterruptedException ie) {
                throw new PicardException("Interrupted while waiting for threads to finished writing.", ie);
            }
            if (!producerSuccceeded.get()) {
                throw new PicardException("Error reading or merging inputs.");
            }
            if (!consumerSuccceeded.get()) {
                throw new PicardException("Error writing output");
            }
        }
        else {
            for (long numRecords = 1; iterator.hasNext(); ++numRecords) {
                final SAMRecord record = iterator.next();
                out.addAlignment(record);
                if (numRecords % PROGRESS_INTERVAL == 0) {
                    log.info(numRecords + " records read.");
                }
            }

        }

        log.info("Finished reading inputs.");
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
