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
package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.FilteringIterator;
import htsjdk.samtools.filter.ReadNameFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

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
@CommandLineProgramProperties(
        usage = "Produces a new SAM or BAM file by including or excluding aligned reads " +
                "or a list of reads names supplied in the READ_LIST_FILE from the INPUT SAM or BAM file.\n",
        usageShort = "Creates a new SAM or BAM file by including or excluding aligned reads",
        programGroup = SamOrBam.class
)
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
        final SamReader inputReader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileHeader.SortOrder inputSortOrder = inputReader.getFileHeader().getSortOrder();
        final SAMFileHeader outputHeader = inputReader.getFileHeader();
        if (SORT_ORDER != null) {
            outputHeader.setSortOrder(SORT_ORDER);
        }
        final boolean presorted = inputSortOrder.equals(outputHeader.getSortOrder());
        log.info("Filtering [presorted=" + presorted + "] " + INPUT.getName() + " -> OUTPUT=" +
                OUTPUT.getName() + " [sortorder=" + outputHeader.getSortOrder().name() + "]");

        // create OUTPUT file
        final SAMFileWriter outputWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outputHeader, presorted, OUTPUT);

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Written");

        while (filteringIterator.hasNext()) {
            final SAMRecord rec = filteringIterator.next();
            outputWriter.addAlignment(rec);
            progress.record(rec);
        }

        filteringIterator.close();
        outputWriter.close();
        CloserUtil.close(inputReader);
        log.info(new DecimalFormat("#,###").format(progress.getCount()) + " SAMRecords written to " + OUTPUT.getName());
    }

    /**
     * Write out a file of read names for debugging purposes.
     *
     * @param samOrBamFile The SAM or BAM file for which we are going to write out a file of its
     *                     containing read names
     */
    private void writeReadsFile(final File samOrBamFile) throws IOException {
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(samOrBamFile);
        final File readsFile =
                new File(OUTPUT.getParentFile(), IOUtil.basename(samOrBamFile) + ".reads");
        IOUtil.assertFileIsWritable(readsFile);
        final BufferedWriter bw = IOUtil.openFileForBufferedWriting(readsFile, false);

        for (final SAMRecord rec : reader) {
            bw.write(rec.toString() + "\n");
        }

        bw.close();
        reader.close();
        IOUtil.assertFileIsReadable(readsFile);
    }

    @Override
    protected int doWork() {

        try {
            IOUtil.assertFileIsReadable(INPUT);
            IOUtil.assertFileIsWritable(OUTPUT);
            if (WRITE_READS_FILES) writeReadsFile(INPUT);

            switch (FILTER) {
                case includeAligned:
                    filterReads(new FilteringIterator(SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT).iterator(),
                            new AlignedFilter(true), true));
                    break;
                case excludeAligned:
                    filterReads(new FilteringIterator(SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT).iterator(),
                            new AlignedFilter(false), true));
                    break;
                case includeReadList:
                    filterReads(new FilteringIterator(SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT).iterator(),
                            new ReadNameFilter(READ_LIST_FILE, true)));
                    break;
                case excludeReadList:
                    filterReads(new FilteringIterator(SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT).iterator(),
                            new ReadNameFilter(READ_LIST_FILE, false)));
                    break;
                default:
                    throw new UnsupportedOperationException(FILTER.name() + " has not been implemented!");
            }

            IOUtil.assertFileIsReadable(OUTPUT);
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
