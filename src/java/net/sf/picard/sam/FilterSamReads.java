/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

/**
 * $Id$
 */
package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.filter.AlignedFilter;
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.ReadNameFilter;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 * From a SAM or BAM file, produce a new SAM or BAM by filtering aligned reads or a list of read
 * names provided in a file (one readname per line)
 */
public class FilterSamReads extends CommandLineProgram {

    private static final Log log = Log.getInstance(FilterSamReads.class);

    @Usage
    public String USAGE =
        "Produces a new SAM or BAM file by including or excluding aligned reads " +
            "or a list of specific reads names from a SAM or BAM file.\n";

    @Option(doc = "The SAM or BAM file that will be read excluded",
        optional = false,
        shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc = "Include aligned reads in the OUTPUT SAM or BAM\n" +
        "Note that *both* first and second of paired reads must be aligned to be included in the OUTPUT SAM or BAM",
        mutex = { "EXCLUDE_ALIGNED", "EXCLUDE_READS", "INCLUDE_READS" }, optional = false,
        shortName = "IA")
    public boolean INCLUDE_ALIGNED = false;

    @Option(doc = "Exclude aligned reads from the OUTPUT SAM or BAM\n" +
        "Note that *both* first and second of pair must be aligned to be excluded from the OUTPUT SAM or BAM",
        mutex = { "INCLUDE_ALIGNED", "EXCLUDE_READS", "INCLUDE_READS" }, optional = false,
        shortName = "EA")
    public boolean EXCLUDE_ALIGNED = false;

    @Option(
        doc = "File containing read names that will be included in the OUTPUT SAM or BAM",
        mutex = { "EXCLUDE_READS", "INCLUDE_ALIGNED", "EXCLUDE_ALIGNED" }, optional = false,
        shortName = "IR")
    public File INCLUDE_READS = null;

    @Option(
        doc = "File containing read names that will be excluded from the OUTPUT SAM or BAM",
        mutex = { "INCLUDE_READS", "INCLUDE_ALIGNED", "EXCLUDE_ALIGNED" }, optional = false,
        shortName = "ER")
    public File EXCLUDE_READS = null;

    @Option(
        doc = "SortOrder of the OUTPUT SAM or BAM file, otherwise use the SortOrder of the INPUT file.",
        optional = true, shortName = "SO")
    public SAMFileHeader.SortOrder SORT_ORDER;

    @Option(
        doc = "Create .reads files (for debugging purposes)",
        optional = true)
    public boolean WRITE_READS_FILE = true;

    @Option(doc = "SAM or BAM file to write read excluded results to",
        optional = false, shortName = "O")
    public File OUTPUT;

    private void filterReads(final FilteringIterator filteringIterator) {
        // get OUTPUT header from INPUT and owerwrite it if necessary
        final SAMFileReader inputReader = new SAMFileReader(INPUT);
        final SAMFileHeader.SortOrder inputSortOrder = inputReader.getFileHeader().getSortOrder();
        final SAMFileHeader outPutHeader = inputReader.getFileHeader();
        if (SORT_ORDER != null) {
            outPutHeader.setSortOrder(SORT_ORDER);
        }
        final boolean presorted = inputSortOrder.equals(outPutHeader.getSortOrder());
        log.info("Filtering [presorted=" + presorted + "] " + INPUT.getName() + " -> OUTPUT=" +
            OUTPUT.getName() + " [sortorder=" + outPutHeader.getSortOrder().name() + "]");

        // create OUTPUT file
        final SAMFileWriter outputWriter =
            new SAMFileWriterFactory().makeSAMOrBAMWriter(outPutHeader, presorted, OUTPUT);

        int count = 0;
        while (filteringIterator.hasNext()) {
            outputWriter.addAlignment(filteringIterator.next());
            count++;
            if (count != 0 && (count % 1000000 == 0)) {
                log.info(new DecimalFormat("#,###").format(count) + " SAMRecords written to " +
                    OUTPUT.getName());
            }
        }

        filteringIterator.close();
        outputWriter.close();
        inputReader.close();
        log.info(new DecimalFormat("#,###").format(count) + " SAMRecords written to " +
            OUTPUT.getName());
    }

    private void writeReadsFile() {
        try {
            final SAMFileReader reader = new SAMFileReader(OUTPUT);
            final File readsFile = new File(OUTPUT.getParentFile(), IoUtil.basename(OUTPUT) + ".reads");
            IoUtil.assertFileIsWritable(readsFile);
            final BufferedWriter bw = IoUtil.openFileForBufferedWriting(readsFile, false);

            for (final SAMRecord rec : reader) {
                bw.write(rec.toString() + "\n");
            }

            bw.close();
            reader.close();
            IoUtil.assertFileIsReadable(readsFile);
        } catch (IOException e) {
            log.error(e);
        }
    }

    @Override
    protected int doWork() {

        if (INPUT.equals(OUTPUT)) {
            throw new PicardException("INPUT file and OUTPUT file must differ!");
        }

        try {
            IoUtil.assertFileIsReadable(INPUT);
            IoUtil.assertFileIsWritable(OUTPUT);

            if (INCLUDE_READS != null) {
                filterReads(new FilteringIterator(new SAMFileReader(INPUT).iterator(),
                    new ReadNameFilter(INCLUDE_READS, true)));
            } else if (EXCLUDE_READS != null) {
                filterReads(new FilteringIterator(new SAMFileReader(INPUT).iterator(),
                    new ReadNameFilter(EXCLUDE_READS, false)));
            } else if (INCLUDE_ALIGNED) {
                filterReads(new FilteringIterator(new SAMFileReader(INPUT).iterator(),
                    new AlignedFilter(true), true));
            } else if (EXCLUDE_ALIGNED) {
                filterReads(new FilteringIterator(new SAMFileReader(INPUT).iterator(),
                    new AlignedFilter(false), true));
            }

            IoUtil.assertFileIsReadable(OUTPUT);
            if (WRITE_READS_FILE) writeReadsFile();
            return 0;

        } catch (Exception e) {
            if (OUTPUT.exists() && !OUTPUT.delete()) {
                log.warn("Failed to delete " + OUTPUT.getAbsolutePath());
            }

            log.error(e, "Failed to filter " + INPUT.getName());
            return 1;
        }
    }

    /**
     * Stock main method.
     *
     * @param args main arguements
     */
    public static void main(final String[] args) {
        System.exit(new FilterSamReads().instanceMain(args));
    }

}