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
 * 
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
import htsjdk.samtools.filter.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * From a SAM or BAM file, produce a new SAM or BAM by filtering aligned reads or a list of read
 * names provided in a file (one readname per line)
 * <p/>
 * $Id$
 */
@CommandLineProgramProperties(
        summary =  FilterSamReads.USAGE_SUMMARY + FilterSamReads.USAGE_DETAILS,
        oneLineSummary = FilterSamReads.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class FilterSamReads extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Subset read data from a SAM or BAM file";
    static final String USAGE_DETAILS = "This tool takes a SAM or BAM file and subsets it to a new file that either excludes or " +
            "only includes either aligned or unaligned reads (set using FILTER), or specific reads based on a list of reads names " +
            "supplied in the READ_LIST_FILE.  " +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar FilterSamReads \\<br /> " +
            "      I=input.bam \\ <br /> " +
            "      O=output.bam \\<br /> " +
            "      READ_LIST_FILE=read_names.txt" +
            "      FILTER=filter_value" +
            "</pre> " +
            "For information on the SAM format, please see: http://samtools.sourceforge.net" +
            "<hr />";
    private static final Log log = Log.getInstance(FilterSamReads.class);
    
    protected /* <- used in test */ enum Filter {
        includeAligned("OUTPUT SAM/BAM will contain aligned reads only. INPUT SAM/BAM must be in queryname SortOrder. (Note that *both* first and second of paired reads must be aligned to be included in the OUTPUT SAM or BAM)"),
        excludeAligned("OUTPUT SAM/BAM will contain un-mapped reads only. INPUT SAM/BAM must be in queryname SortOrder. (Note that *both* first and second of pair must be aligned to be excluded from the OUTPUT SAM or BAM)"),
        includeReadList("OUTPUT SAM/BAM will contain reads that are supplied in the READ_LIST_FILE file"),
        excludeReadList("OUTPUT bam will contain reads that are *not* supplied in the READ_LIST_FILE file"),
    	includeJavascript("OUTPUT bam will contain reads that hava been accepted by the JAVASCRIPT_FILE script."),
        includePairedIntervals("OUTPUT SAM/BAM will contain any reads (and their mate) that overlap with an interval. INPUT SAM/BAM and INTERVAL_LIST must be in coordinate SortOrder. Only aligned reads will be output.");
        private final String description;

        Filter(final String description) {
            this.description = description;
        }

        @Override
        public String toString() {
            return this.name() + " [" + description + "]";
        }
    }

    @Argument(doc = "The SAM or BAM file that will be filtered.",
            optional = false,
            shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Filter.", optional = false)
    public Filter FILTER = null;

    @Argument(doc = "Read List File containing reads that will be included or excluded from the OUTPUT SAM or BAM file.",
            optional = true,
            shortName = "RLF")
    public File READ_LIST_FILE;

    @Argument(doc = "Interval List File containing intervals that will be included or excluded from the OUTPUT SAM or BAM file.",
            optional = true,
            shortName = "IL")
    public File INTERVAL_LIST;

    @Argument(
            doc = "SortOrder of the OUTPUT SAM or BAM file, otherwise use the SortOrder of the INPUT file.",
            optional = true, shortName = "SO")
    public SAMFileHeader.SortOrder SORT_ORDER;

    @Argument(
            doc = "Create .reads files (for debugging purposes)",
            optional = true)
    public boolean WRITE_READS_FILES = true;

    @Argument(doc = "SAM or BAM file to write read excluded results to",
            optional = false, shortName = "O")
    public File OUTPUT;
    
	@Argument(shortName = "JS",
			doc = "Filters a SAM or BAM file with a javascript expression using the java javascript-engine. "
	        + " The script puts the following variables in the script context: "
	        + " 'record' a SamRecord ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html ) and "
	        + " 'header' a SAMFileHeader ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html )."
	        + " Last value of the script should be a boolean to tell wether we should accept or reject the record.",
	        optional = true)
	public File JAVASCRIPT_FILE = null;

    
    private void filterReads(final FilteringSamIterator filteringIterator) {

        // get OUTPUT header from INPUT and overwrite it if necessary
        final SAMFileHeader fileHeader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(INPUT);
        final SAMFileHeader.SortOrder inputSortOrder = fileHeader.getSortOrder();
        if (SORT_ORDER != null) {
            fileHeader.setSortOrder(SORT_ORDER);
        }

        if (FILTER == Filter.includePairedIntervals && fileHeader.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new UnsupportedOperationException("Input must be coordinate sorted to use includePairedIntervals");
        }

        final boolean presorted = inputSortOrder.equals(fileHeader.getSortOrder());
        log.info("Filtering [presorted=" + presorted + "] " + INPUT.getName() + " -> OUTPUT=" +
                OUTPUT.getName() + " [sortorder=" + fileHeader.getSortOrder().name() + "]");

        // create OUTPUT file
        final SAMFileWriter outputWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(fileHeader, presorted, OUTPUT);

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Written");

        while (filteringIterator.hasNext()) {
            final SAMRecord rec = filteringIterator.next();
            outputWriter.addAlignment(rec);
            progress.record(rec);
        }

        filteringIterator.close();
        outputWriter.close();
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

    private List<Interval> getIntervalList (final File intervalFile) throws IOException {
        IOUtil.assertFileIsReadable(intervalFile);
        return IntervalList.fromFile(intervalFile).getIntervals();
    }

    @Override
    protected int doWork() {

        try {
            IOUtil.assertFileIsReadable(INPUT);
            IOUtil.assertFileIsWritable(OUTPUT);
            if (WRITE_READS_FILES) writeReadsFile(INPUT);

            List<Interval> intervalList = new ArrayList<>();

            if (INTERVAL_LIST != null) {
                intervalList = getIntervalList(INTERVAL_LIST);
            }

            final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
            final FilteringSamIterator filteringIterator;
            
            switch (FILTER) {
                case includeAligned:
                	filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new AlignedFilter(true), true);
                    break;
                case excludeAligned:
                	filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new AlignedFilter(false), true);
                    break;
                case includeReadList:
                	filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new ReadNameFilter(READ_LIST_FILE, true));
                    break;
                case excludeReadList:
                	filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new ReadNameFilter(READ_LIST_FILE, false));
                    break;
                case includeJavascript:
                	filteringIterator = new FilteringSamIterator(samReader.iterator(),
                			new JavascriptSamRecordFilter(
                			        JAVASCRIPT_FILE,
                					samReader.getFileHeader()));
                    break;
                case includePairedIntervals:
                    filteringIterator = new FilteringSamIterator(samReader.iterator(),
                            new IntervalKeepPairFilter(intervalList), false);
                    break;
                default:
                    throw new UnsupportedOperationException(FILTER.name() + " has not been implemented!");
            }
            
            filterReads(filteringIterator);

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

        if (FILTER.equals(Filter.includePairedIntervals) && INTERVAL_LIST == null) {
            return new String[]{"A INTERVAL_LIST must be specified when using the " + FILTER.name() + " option"};
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
