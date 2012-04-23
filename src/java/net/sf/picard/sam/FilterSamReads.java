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
 * <p/>
 * $Id$
 */
public class FilterSamReads extends CommandLineProgram {

    private static final Log log = Log.getInstance(FilterSamReads.class);

    private static enum Filter {
        includeAligned("OUTPUT SAM/BAM will contain aligned reads only. INPUT SAM/BAM must be in queryname SortOrder. (Note that *both* first and second of paired reads must be aligned to be included in the OUTPUT SAM or BAM)"),
        excludeAligned("OUTPUT SAM/BAM will contain un-mapped reads only. INPUT SAM/BAM must be in queryname SortOrder. (Note that *both* first and second of pair must be aligned to be excluded from the OUTPUT SAM or BAM)"),
        includeReadList("OUTPUT SAM/BAM will contain reads that are supplied in the READ_LIST_FILE file"),
        excludeReadList("OUTPUT bam will contain reads that are *not* supplied in the READ_LIST_FILE file");

        private final String description;

        Filter(final String description) {
            this.description = description;
        }

        @Override
        public String toString() {
            return this.name() + " [" + description + "]";
        }
    }

    @Usage
    public String USAGE =
            "Produces a new SAM or BAM file by including or excluding aligned reads " +
                    "or a list of reads names supplied in the READ_LIST_FILE from the INPUT SAM or BAM file.\n";

    @Option(doc = "The SAM or BAM file that will be filtered.",
        optional = false,
        shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc = "Filter.", optional = false)
    public Filter FILTER = null;

    @Option(doc = "Read List File containing reads that will be included or excluded from the OUTPUT SAM or BAM file.",
            optional = true,
            shortName = "RLF")
    public File READ_LIST_FILE;

    @Option(
        doc = "SortOrder of the OUTPUT SAM or BAM file, otherwise use the SortOrder of the INPUT file.",
        optional = true, shortName = "SO")
    public SAMFileHeader.SortOrder SORT_ORDER;

    @Option(
        doc = "Create .reads files (for debugging purposes)",
        optional = true)
    public boolean WRITE_READS_FILES = true;

    @Option(doc = "SAM or BAM file to write read excluded results to",
        optional = false, shortName = "O")
    public File OUTPUT;

    private void filterReads(final FilteringIterator filteringIterator) {

        // get OUTPUT header from INPUT and owerwrite it if necessary
        final SAMFileReader inputReader = new SAMFileReader(INPUT);
        final SAMFileHeader.SortOrder inputSortOrder = inputReader.getFileHeader().getSortOrder();
        final SAMFileHeader outputHeader = inputReader.getFileHeader();
        if (SORT_ORDER != null) {
            outputHeader.setSortOrder(SORT_ORDER);
        }
        final boolean presorted = inputSortOrder.equals(outputHeader.getSortOrder());
        log.info("Filtering [presorted=" + presorted + "] " + INPUT.getName() + " -> OUTPUT=" +
            OUTPUT.getName() + " [sortorder=" + outputHeader.getSortOrder().name() + "]");

        // create OUTPUT file
        final SAMFileWriter outputWriter =
            new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader, presorted, OUTPUT);

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

    /**
     * Write out a file of read names for debugging purposes.
     *
     * @param samOrBamFile The SAM or BAM file for which we are going to write out a file of its
     *                     containing read names
     */
    private void writeReadsFile(final File samOrBamFile) throws IOException {
        final SAMFileReader reader = new SAMFileReader(samOrBamFile);
        final File readsFile =
            new File(OUTPUT.getParentFile(), IoUtil.basename(samOrBamFile) + ".reads");
        IoUtil.assertFileIsWritable(readsFile);
        final BufferedWriter bw = IoUtil.openFileForBufferedWriting(readsFile, false);

        for (final SAMRecord rec : reader) {
            bw.write(rec.toString() + "\n");
        }

        bw.close();
        reader.close();
        IoUtil.assertFileIsReadable(readsFile);
    }

    @Override
    protected int doWork() {

        try {
            IoUtil.assertFileIsReadable(INPUT);
            IoUtil.assertFileIsWritable(OUTPUT);
            if (WRITE_READS_FILES) writeReadsFile(INPUT);

            switch (FILTER) {
                case includeAligned:
                    filterReads(new FilteringIterator(new SAMFileReader(INPUT).iterator(),
                    new AlignedFilter(true), true));
                    break;
                case excludeAligned:
                    filterReads(new FilteringIterator(new SAMFileReader(INPUT).iterator(),
                    new AlignedFilter(false), true));
                    break;
                case includeReadList:
                    filterReads(new FilteringIterator(new SAMFileReader(INPUT).iterator(),
                    new ReadNameFilter(READ_LIST_FILE, true)));
                    break;
                case excludeReadList:
                    filterReads(new FilteringIterator(new SAMFileReader(INPUT).iterator(),
                    new ReadNameFilter(READ_LIST_FILE, false)));
                    break;
                default:
                    throw new UnsupportedOperationException(FILTER.name() + " has not been implemented!");
            }

            IoUtil.assertFileIsReadable(OUTPUT);
            if (WRITE_READS_FILES) writeReadsFile(OUTPUT);
            return 0;

        } catch (Exception e) {
            if (OUTPUT.exists() && !OUTPUT.delete()) {
                log.warn("Failed to delete " + OUTPUT.getAbsolutePath());
            }

            log.error(e, "Failed to filter " + INPUT.getName());
            return 1;
        }
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (INPUT.equals(OUTPUT)) {
            return new String[]{"INPUT file and OUTPUT file must differ!"};
        }

        if ((FILTER.equals(Filter.includeReadList) ||
                FILTER.equals(Filter.excludeReadList)) &&
                READ_LIST_FILE == null) {
            return new String[]{"A READ_LIST_FILE must be specified when using the " + FILTER.name() + " option"};

        }

        return super.customCommandLineValidation();
    }

    /**
     * Stock main method.
     *
     * @param args main arguments
     */
    public static void main(final String[] args) {
        System.exit(new FilterSamReads().instanceMain(args));
    }

}