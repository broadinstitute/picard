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
import net.sf.picard.util.CollectionUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeSet;

/**
 * Filters out mapped or unmapped reads from an INPUT sam/bam file and writes
 * out the remaining reads to the specified OUTPUT sam/bam file.
 */
public class FilterSamReads extends CommandLineProgram {

private static final Log log = Log.getInstance(FilterSamReads.class);

public enum ReadFilterType {INCLUDE, EXCLUDE}
public enum ReadMappingType {MAPPED, UNMAPPED}

@Usage
public String USAGE =
    "Produces a new SAM or BAM file by filtering (including or excluding) " +
        "mapped or unmapped reads from a specified SAM or BAM file";

@Option(doc = "The SAM or BAM file that will be read excluded",
    optional = false,
    shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
public File INPUT;

@Option(doc = "SAM or BAM file to write read excluded results to",
    optional = false, shortName = "O")
public File OUTPUT;

@Option(
    doc = "Determines if reads should be included or excluded from OUTPUT SAM or BAM file",
    optional = false, shortName = "RFT")
public ReadFilterType READ_FILTER_TYPE;

@Option(doc = "A file of read names that will be excluded from the SAM or BAM",
    mutex = { "READ_MAPPING_TYPE" }, optional = false, shortName = "RLF")
public File READNAME_LIST_FILE;

@Option(doc = "Exclude mapped or umapped reads from the SAM or BAM",
    mutex = { "READNAME_LIST_FILE" }, optional = false, shortName = "RMT")
public ReadMappingType READ_MAPPING_TYPE;

@Option(doc = "SortOrder of the OUTPUT SAM or BAM file",
    optional = true, shortName = "SO")
public SAMFileHeader.SortOrder SORT_ORDER;

private void filterByReadNameList() {

    // read RLF file into sorted Tree Set - this will remove dups
    IoUtil.assertFileIsReadable(READNAME_LIST_FILE);
    IoUtil.assertFileSizeNonZero(READNAME_LIST_FILE);
    final BufferedReader is;

    try {
        is = IoUtil.openFileForBufferedReading(READNAME_LIST_FILE);
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

    final SAMFileReader reader = new SAMFileReader(INPUT);
    final SAMFileHeader header = reader.getFileHeader();

    if (SORT_ORDER != null) {
        header.setSortOrder(SORT_ORDER);
    }

    log.info("SORT_ORDER of OUTPUT file will be " + header.getSortOrder().name());

    final SAMFileWriter writer =
        new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);

    int count = 0;
    for (final SAMRecord rec : reader) {

        if (READ_FILTER_TYPE.equals(ReadFilterType.INCLUDE)) {
            if (readNameFilterSet.contains(rec.getReadName())) {
                writer.addAlignment(rec);
                readNameFilterSet.remove(rec.getReadName());
            }
        } else {
            if (readNameFilterSet.contains(rec.getReadName())) {
                readNameFilterSet.remove(rec.getReadName());
            } else {
                writer.addAlignment(rec);
            }
        }

        if (++count % 1000000 == 0) {
            log.info(new DecimalFormat("#,###").format(count) +
                " SAMRecords written to " + OUTPUT.getName());
        }
    }

    reader.close();
    writer.close();

    if (!readNameFilterSet.isEmpty()) {
        throw new PicardException(readNameFilterSet.size() +
            " reads could not be found within INPUT=" + INPUT.getPath() + " [" +
            CollectionUtil.join(readNameFilterSet, ",") + "]");
    }

    IoUtil.assertFileIsReadable(OUTPUT);
}

private void validatePair(final SAMRecord firstOfPair,
                          final SAMRecord secondOfPair) {
    // Validate that paired reads arrive as first of pair followed by second of pair
    if (!firstOfPair.getReadName().equals(secondOfPair.getReadName()))
        throw new PicardException(
            "Second read from pair not found: " + firstOfPair.getReadName() +
                ", " + secondOfPair.getReadName());

    if (!firstOfPair.getFirstOfPairFlag()) throw new PicardException(
        "First record in unmapped bam is not first of pair: " +
            firstOfPair.getReadName());
    if (!secondOfPair.getReadPairedFlag()) throw new PicardException(
        "Second record is not marked as paired: " + secondOfPair.getReadName());
    if (!secondOfPair.getSecondOfPairFlag()) throw new PicardException(
        "Second record is not second of pair: " + secondOfPair.getReadName());
}

private void filterByMappingType() {
    final SAMFileReader reader = new SAMFileReader(INPUT);
    final SAMFileHeader header = reader.getFileHeader();
    if (SORT_ORDER != null) {
        header.setSortOrder(SORT_ORDER);
    }

    log.info("SORT_ORDER of OUTPUT file will be " + header.getSortOrder().name());

    final SAMFileWriter writer =
        new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);

    int count = 0;
    final CloseableIterator<SAMRecord> it = reader.iterator();
    while (it.hasNext()) {
        final SAMRecord rec = it.next();
        boolean writeToOutput = false;

        if (rec.getReadPairedFlag()) {
            final SAMRecord secondOfPair = it.next();
            validatePair(rec, secondOfPair);

            if ((READ_FILTER_TYPE.equals(ReadFilterType.INCLUDE) &&
                READ_MAPPING_TYPE.equals(ReadMappingType.MAPPED)) ||
                (READ_FILTER_TYPE.equals(ReadFilterType.EXCLUDE) &&
                    READ_MAPPING_TYPE.equals(ReadMappingType.UNMAPPED))) {

                // include mapped or exclude unmapped
                if (!rec.getReadUnmappedFlag() &&
                    !secondOfPair.getReadUnmappedFlag()) {
                    writeToOutput = true;
                }
            } else {
                // include unmapped or exclude mapped
                if (rec.getReadUnmappedFlag() &&
                    secondOfPair.getReadUnmappedFlag()) {
                    writeToOutput = true;
                }
            }

            if (writeToOutput) {
                writer.addAlignment(rec);
                writer.addAlignment(secondOfPair);
            } else {
                log.debug(
                    "Skipping " + READ_FILTER_TYPE + " " + READ_MAPPING_TYPE +
                        " " + rec.toString() + " and " +
                        secondOfPair.toString());
            }

        } else {
            if ((READ_FILTER_TYPE.equals(ReadFilterType.INCLUDE) &&
                READ_MAPPING_TYPE.equals(ReadMappingType.MAPPED)) ||
                (READ_FILTER_TYPE.equals(ReadFilterType.EXCLUDE) &&
                    READ_MAPPING_TYPE.equals(ReadMappingType.UNMAPPED))) {

                // include mapped or exclude unmapped
                if (!rec.getReadUnmappedFlag()) {
                    writeToOutput = true;
                }
            } else {
                // include unmapped or exclude mapped
                if (rec.getReadUnmappedFlag()) {
                    writeToOutput = true;
                }
            }

            if (writeToOutput) {
                writer.addAlignment(rec);
            } else {
                log.info(
                    "Skipping " + READ_FILTER_TYPE + " " + READ_MAPPING_TYPE +
                        " " + rec.toString());
            }
        }

        if (writeToOutput && ++count % 1000000 == 0) {
            log.info(new DecimalFormat("#,###").format(count) +
                " SAMRecords written to " + OUTPUT.getName());
        }
    }

    reader.close();
    writer.close();
}

@Override
protected int doWork() {

    if (INPUT.equals(OUTPUT)) {
        throw new PicardException("INPUT file and OUTPUT file must differ!");
    }

    try {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

        if (READNAME_LIST_FILE != null) {
            filterByReadNameList();
        } else {
            filterByMappingType();
        }

        IoUtil.assertFileIsReadable(OUTPUT);

        return 0;

    } catch (Exception e) {
        if (!OUTPUT.delete()) {
           log.warn("Failed to delete " + OUTPUT.getAbsolutePath());
        }

        log.error("Failed to filter " + INPUT.getName() + " " + e.getMessage(), e);
        return 1;
    }
}

/**
 * Stock main method.
 * @param args main arguements
 */
public static void main(final String[] args) {
    System.exit(new FilterSamReads().instanceMain(args));
}

}