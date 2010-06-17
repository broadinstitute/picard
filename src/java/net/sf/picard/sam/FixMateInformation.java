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

import net.sf.picard.cmdline.Usage;

import java.io.File;
import java.io.IOException;

import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.RuntimeIOException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.util.Log;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.util.SortingCollection;

/**
 * Class to fix mate pair information for all reads in a SAM file.  Will run in fairly limited
 * memory unless there are lots of mate pairs that are far apart from each other in the file.
 *
 * @author Tim Fennell
 */
public class FixMateInformation extends CommandLineProgram {
    @Usage public final String USAGE = "Ensure that all mate-pair information is in sync between each read " +
            " and it's mate pair.  If no OUTPUT file is supplied then the output is written to a temporary file " +
            " and then copied over the INPUT file.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input file to fix.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional=true,
            doc="The output file to write to. If no output file is supplied, the input file is overwritten.")
    public File OUTPUT;

    @Option(shortName=StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, optional=true,
    doc="Optional sort order if the OUTPUT file should be sorted differently than the INPUT file.")
    public SortOrder SORT_ORDER;

    private static final Log log = Log.getInstance(FixMateInformation.class);

    public static void main(final String[] args) {
        new FixMateInformation().instanceMainWithExit(args);
    }

    protected int doWork() {
        // Open up the input
        INPUT = INPUT.getAbsoluteFile();
        if (OUTPUT != null) OUTPUT = OUTPUT.getAbsoluteFile();        
        IoUtil.assertFileIsReadable(INPUT);
        final SAMFileReader in = new SAMFileReader(INPUT);
        final SAMFileHeader header = in.getFileHeader();

        // Deal with the various sorting complications
        final boolean inputIsQueryNameSorted = in.getFileHeader().getSortOrder() == SortOrder.queryname;
        final SortOrder outputSortOrder = SORT_ORDER == null ? in.getFileHeader().getSortOrder() : SORT_ORDER;
        final boolean outputIsQueryNameSorted = outputSortOrder == SortOrder.queryname;
        header.setSortOrder(outputSortOrder);

        // Decide where to write the fixed file - into the specified output file
        // or into a temporary file that will overwrite the INPUT file eventually
        final boolean differentOutputSpecified = OUTPUT != null;

        if (differentOutputSpecified) {
            IoUtil.assertFileIsWritable(OUTPUT);
        }
        else {
            final File dir = INPUT.getParentFile();
            try {
                IoUtil.assertDirectoryIsWritable(dir);
                IoUtil.assertFileIsWritable(INPUT);
                OUTPUT = File.createTempFile(INPUT.getName() + ".being_fixed.", ".bam", dir);
            }
            catch (IOException ioe) {
                throw new RuntimeIOException("Could not create tmp file in " + dir.getAbsolutePath());
            }
        }

        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(header,
                                                                                header.getSortOrder() == SortOrder.queryname,
                                                                                OUTPUT);

        // Get the input records sorted by query name
        final PeekableIterator<SAMRecord> iterator;
        if (inputIsQueryNameSorted) {
            iterator = new PeekableIterator<SAMRecord>(in.iterator());
        }
        else {
            log.info("Sorting input into queryname order.");
            final SortingCollection<SAMRecord> sorter = SortingCollection.newInstance(SAMRecord.class,
                                                                                      new BAMRecordCodec(in.getFileHeader()),
                                                                                      new SAMRecordQueryNameComparator(),
                                                                                      MAX_RECORDS_IN_RAM,
                                                                                      TMP_DIR);

            for (final SAMRecord rec : in) {
                sorter.add(rec);
            }

            iterator = new PeekableIterator<SAMRecord>(sorter.iterator());
            in.close();

        }

        // Now go through the input file and fix up the records
        int count = 0;
        while (iterator.hasNext()) {
            final SAMRecord rec1 = iterator.next();
            final SAMRecord rec2 = iterator.hasNext() ? iterator.peek() : null;

            if (rec2 != null && rec1.getReadName().equals(rec2.getReadName())) {
                iterator.next(); // consume the peeked record
                SamPairUtil.setMateInfo(rec1, rec2, header);
                out.addAlignment(rec1);
                out.addAlignment(rec2);
                count += 2;
            }
            else {
                out.addAlignment(rec1);
                ++count;
            }

            if (count % 1000000 == 0) {
                log.info("Processed " + count + " records.");
            }
        }

        if (header.getSortOrder() == SortOrder.queryname) {
            log.info("Closing output file.");
        }
        else {
            log.info("Finished processing reads; re-sorting output file.");
        }

        out.close();

        // Lastly if we're fixing in place, swap the files
        if (!differentOutputSpecified) {
            log.info("Replacing input file with fixed file.");

            final File old = new File(INPUT.getParentFile(), INPUT.getName() + ".old");
            if (!old.exists() && INPUT.renameTo(old)) {
                if (OUTPUT.renameTo(INPUT)) {

                    if (!old.delete()) {
                        log.warn("Could not delete old file: " + old.getAbsolutePath());
                        return 1;
                    }
                }
                else {
                    log.error("Could not move new file to " + INPUT.getAbsolutePath());
                    log.error("Input file preserved as: " + old.getAbsolutePath());
                    log.error("New file preserved as: " + OUTPUT.getAbsolutePath());
                    return 1;
                }
            }
            else {
                log.error("Could not move input file out of the way: " + INPUT.getAbsolutePath());

                if (!OUTPUT.delete()) {
                    log.error("Could not delete temporary file: " + OUTPUT.getAbsolutePath());
                }

                return 1;
            }

        }

        return 0;
    }
}
