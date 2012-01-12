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
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeSet;

/**
 * From a SAM or BAM file, produce a new SAM or BAM by filtering aligned reads or a list of read
 * names provided in a file (one readname per line)
 */
public class FilterSamReads extends CommandLineProgram {

    private static final Log log = Log.getInstance(FilterSamReads.class);

    @Usage
    public String USAGE =
        "Produces a new SAM or BAM file by including or excluding aligned reads " +
            "or a list of specific reads names from a SAM or BAM file.\n" +
            "Note that paired reads are catagorized as aligned when both reads are aligned";

    @Option(doc = "The SAM or BAM file that will be read excluded",
        optional = false,
        shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc = "Include aligned reads in the OUTPUT SAM or BAM",
        mutex = { "EXCLUDE_ALIGNED", "EXCLUDE_READS", "INCLUDE_READS" }, optional = false,
        shortName = "IA")
    public boolean INCLUDE_ALIGNED = false;

    @Option(doc = "Exclude aligned reads from the OUTPUT SAM or BAM",
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

    @Option(doc = "SortOrder of the OUTPUT SAM or BAM file, " +
        "otherwise use the SortOrder of the INPUT file.",
        optional = true, shortName = "SO")
    public SAMFileHeader.SortOrder SORT_ORDER;

    @Option(doc = "SAM or BAM file to write read excluded results to",
        optional = false, shortName = "O")
    public File OUTPUT;

    /**
     * If necessary produce a new bam file sorted by queryname
     *
     * @param fileToSort File requiring sorting
     *
     * @return the queryname sorted file
     */
    private File sortByQueryName(final File fileToSort) {
        File outputFile = fileToSort;

        if (outputFile != null) {
            // sort INPUT file by queryname
            final SAMFileReader reader = new SAMFileReader(fileToSort);
            final SAMFileHeader header = reader.getFileHeader();

            if (header.getSortOrder() == null ||
                !header.getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
                // if not sorted by queryname then sort it by queryname
                log.info("Sorting " + fileToSort.getName() + " by queryname");

                outputFile = new File(OUTPUT.getParentFile(),
                    IoUtil.basename(fileToSort) + ".queryname.sorted.bam");
                IoUtil.assertFileIsWritable(outputFile);

                header.setSortOrder(SAMFileHeader.SortOrder.queryname);
                final SAMFileWriter writer =
                    new SAMFileWriterFactory().makeBAMWriter(header, false, outputFile);

                int count = 0;
                for (final SAMRecord rec : reader) {
                    writer.addAlignment(rec);

                    if (++count % 1000000 == 0) {
                        log.info(
                            new DecimalFormat("#,###").format(count) + " SAMRecords written to " +
                                outputFile.getName());
                    }
                }

                reader.close();
                writer.close();
                log.info(new DecimalFormat("#,###").format(count) + " SAMRecords written to " +
                    outputFile.getName());
            }
        }

        return outputFile;
    }

    private SAMRecord obtainAssertedMate(final SAMRecordIterator samRecordIterator,
                                         final SAMRecord firstOfPair) {

        if (samRecordIterator.hasNext()) {

            final SAMRecord secondOfPair = samRecordIterator.next();

            // Validate paired reads arrive as first of pair, then second of pair
            if (!firstOfPair.getReadName().equals(secondOfPair.getReadName())) {
                throw new PicardException(
                    "Second read from pair not found: " + firstOfPair.getReadName() + ", " +
                        secondOfPair.getReadName());
            } else if (!firstOfPair.getFirstOfPairFlag()) {
                throw new PicardException("First record in unmapped bam is not first of pair: " +
                    firstOfPair.getReadName());
            } else if (!secondOfPair.getReadPairedFlag()) {
                throw new PicardException(
                    "Second record is not marked as paired: " + secondOfPair.getReadName());
            } else if (!secondOfPair.getSecondOfPairFlag()) {
                throw new PicardException(
                    "Second record is not second of pair: " + secondOfPair.getReadName());
            } else {
                return secondOfPair;
            }

        } else {
            throw new PicardException(
                "A second record does not exist: " + firstOfPair.getReadName());
        }
    }

    private void filterReads(final File readFilterFile, final boolean include) {

        // read RLF file into sorted Tree Set - this will remove dups
        IoUtil.assertFileIsReadable(readFilterFile);
        IoUtil.assertFileSizeNonZero(readFilterFile);
        final BufferedReader is;

        try {
            is = IoUtil.openFileForBufferedReading(readFilterFile);
        } catch (IOException e) {
            throw new PicardException(e.getMessage(), e);
        }

        final Scanner scanner = new Scanner(is);
        final Set<String> readNameFilterSet = new TreeSet<String>();

        while (scanner.hasNext()) {
            final String line = scanner.nextLine();

            if (!line.trim().isEmpty()) {
                readNameFilterSet.add(line.trim());
            }
        }

        scanner.close();

        // get OUTPUT header from INPUT and owerwrite it if necessary
        final SAMFileReader inputReader = new SAMFileReader(INPUT);
        final SAMFileHeader.SortOrder inputSortOrder = inputReader.getFileHeader().getSortOrder();
        final SAMFileHeader outPutHeader = inputReader.getFileHeader();
        if (SORT_ORDER != null) {
            outPutHeader.setSortOrder(SORT_ORDER);
        }

        final boolean presorted = inputSortOrder.equals(outPutHeader.getSortOrder());
        log.info("Filtering [presorted=" + presorted + "] " + INPUT.getName());

        // create OUTPUT file
        final SAMFileWriter outputWriter =
            new SAMFileWriterFactory().makeSAMOrBAMWriter(outPutHeader, presorted, OUTPUT);

        int count = 0;
        for (final SAMRecord rec : inputReader) {

            if (include) {
                if (readNameFilterSet.contains(rec.getReadName())) {
                    outputWriter.addAlignment(rec);
                    count++;
                }
            } else {
                if (!readNameFilterSet.contains(rec.getReadName())) {
                    outputWriter.addAlignment(rec);
                    count++;
                }
            }

            if (count != 0 && (count % 1000000 == 0)) {
                log.info(new DecimalFormat("#,###").format(count) + " SAMRecords written to " +
                    OUTPUT.getName());
            }
        }

        inputReader.close();
        outputWriter.close();

        log.info(new DecimalFormat("#,###").format(count) + " SAMRecords written to " +
            OUTPUT.getName());
    }

    /**
     * Writes reads to OUTPUT file. If read is paired and INCLUDE_ALIGNED was specified both reads
     * must be aligned to be included in the OUTPUT file. If read is paired and EXCLUDE_ALIGNED was
     * specified at least one of the reads must be unmapped to be included in the OUTPUT file.
     *
     * @param include if true include aligned reads otherwise exclude aligned reads
     */
    private void filterAlignedReads(final boolean include) {

        // get OUTPUT header from INPUT and owerwrite it if necessary
        final SAMFileReader inputReader = new SAMFileReader(INPUT);
        final SAMFileHeader outPutHeader = inputReader.getFileHeader();
        if (SORT_ORDER != null) {
            outPutHeader.setSortOrder(SORT_ORDER);
        }
        inputReader.close();

        log.info("SORT_ORDER of " + OUTPUT.getName() + " OUTPUT file will be " +
            outPutHeader.getSortOrder().name());

        // sort INPUT by queryname
        final File inputQnSortedFile = sortByQueryName(INPUT);
        final SAMFileReader sortedInputReader = new SAMFileReader(inputQnSortedFile);

        final boolean presorted =
            sortedInputReader.getFileHeader().getSortOrder().equals(outPutHeader.getSortOrder());
        log.info("Filtering [presorted=" + presorted + "] " + inputQnSortedFile.getName());

        // create OUTPUT file
        final SAMFileWriter outputWriter =
            new SAMFileWriterFactory().makeSAMOrBAMWriter(outPutHeader, presorted, OUTPUT);

        int count = 0;
        final SAMRecordIterator it = sortedInputReader.iterator();
        while (it.hasNext()) {
            final SAMRecord r1 = it.next();
            boolean writeToOutput = false;

            if (r1.getReadPairedFlag()) {

                final SAMRecord r2 = obtainAssertedMate(it, r1);

                if (include) {
                    // both r1 and r2 must be mapped for it to be sent to the OUTPUT file
                    if (!r1.getReadUnmappedFlag() && !r2.getReadUnmappedFlag()) {
                        writeToOutput = true;
                    }
                } else {
                    // EXCLUDE_ALIGNED - if either r1 or r2 is unmapped send it to OUTPUT file
                    if (r1.getReadUnmappedFlag() || r2.getReadUnmappedFlag()) {
                        writeToOutput = true;
                    }
                }

                if (writeToOutput) {
                    outputWriter.addAlignment(r1);
                    count++;
                    outputWriter.addAlignment(r2);
                    count++;
                } else {
                    log.debug("Skipping " + r1.toString() + " and " + r2.toString());
                }

            } else {
                if (include) {
                    if (!r1.getReadUnmappedFlag()) {
                        writeToOutput = true;
                    }
                } else {
                    // EXCLUDE_ALIGNED - if read is unmapped send it to OUTPUT file
                    if (r1.getReadUnmappedFlag()) {
                        writeToOutput = true;
                    }
                }

                if (writeToOutput) {
                    outputWriter.addAlignment(r1);
                    count++;
                } else {
                    log.info("Skipping " + r1.toString());
                }
            }

            if (count != 0 && (count % 1000000 == 0)) {
                log.info(new DecimalFormat("#,###").format(count) + " SAMRecords written to " +
                    OUTPUT.getName());
            }
        }

        sortedInputReader.close();
        outputWriter.close();

        log.info(new DecimalFormat("#,###").format(count) + " SAMRecords written to " +
            OUTPUT.getName());

        if (!inputQnSortedFile.equals(INPUT)) {
            log.debug("Removing temporary file " + inputQnSortedFile.getAbsolutePath());
            if (!inputQnSortedFile.delete()) {
                throw new PicardException(
                    "Failed to delete " + inputQnSortedFile.getAbsolutePath());
            }
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
                filterReads(INCLUDE_READS, true);
            } else if (EXCLUDE_READS != null) {
                filterReads(EXCLUDE_READS, false);
            } else {
                filterAlignedReads(INCLUDE_ALIGNED);
            }

            IoUtil.assertFileIsReadable(OUTPUT);
            return 0;

        } catch (Exception e) {
            if (!OUTPUT.delete()) {
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