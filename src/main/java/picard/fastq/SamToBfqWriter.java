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
package picard.fastq;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.TagFilter;
import htsjdk.samtools.filter.WholeReadClippedFilter;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.utils.ValidationUtils;
import picard.PicardException;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Class to take unmapped reads in SAM/BAM/CRAM file format and create Maq binary fastq format file(s) --
 * one or two of them, depending on whether it's a paired-end read.  This relies on the unmapped
 * file having all paired reads together in order.
 */
public class SamToBfqWriter {

    private static final int SEED_REGION_LENGTH = 28;
    private static final int MAX_SEED_REGION_NOCALL_FIXES = 2;

    private final File samFile;
    private final String outputPrefix;
    private final String namePrefix;
    private final int nameTrim;
    private final File referenceSequence;
    private final boolean pairedReads;
    private int wrote = 0;
    private int increment = 1;
    private int chunk = 0;
    private BinaryCodec codec1;
    private BinaryCodec codec2;
    private final Log log = Log.getInstance(SamToBfqWriter.class);
    private final boolean includeNonPfReads;
    private final boolean clipAdapters;
    private final Integer basesToWrite;

    /**
     * Constructor
     *
     * @param samFile           the BAM file to read from
     * @param outputPrefix      the directory and file prefix for the binary fastq files
     * @param total             the total number of records that should be written, drawn evenly
     *                          from throughout the file (null for all).
     * @param chunk             the maximum number of records that should be written to any one file
     * @param pairedReads       whether these reads are from  a paired-end run
     * @param namePrefix        The string to be stripped off the read name
     *                          before writing to the bfq file. May be null, in which case
     *                          the name will not be trimmed.
     * @param includeNonPfReads whether to include non pf-reads
     * @param clipAdapters      whether to replace adapters as marked with XT:i clipping position attribute
     */
    public SamToBfqWriter(final File samFile, final File referenceSequence, final String outputPrefix, final Integer total,
                          final Integer chunk, final boolean pairedReads, String namePrefix,
                          boolean includeNonPfReads, boolean clipAdapters, Integer basesToWrite) {

        IOUtil.assertFileIsReadable(samFile);
        this.samFile = samFile;
        this.outputPrefix = outputPrefix;
        this.pairedReads = pairedReads;
        this.referenceSequence = referenceSequence;
        if (total != null) {
            final double writeable = countWritableRecords();
            this.increment = (int) Math.floor(writeable / total.doubleValue());
            if (this.increment == 0) {
                this.increment = 1;
            }
        }
        if (chunk != null) {
            this.chunk = chunk;
        }
        this.namePrefix = namePrefix;
        this.nameTrim = namePrefix != null ? namePrefix.length() : 0;
        this.includeNonPfReads = includeNonPfReads;
        this.clipAdapters = clipAdapters;
        this.basesToWrite = basesToWrite;
    }

    /**
     * Constructor
     *
     * @param samFile           the BAM file to read from
     * @param outputPrefix      the directory and file prefix for the binary fastq files
     * @param pairedReads       whether these reads are from  a paired-end run
     * @param namePrefix        the barcode of the run (to be stripped off the read name
     *                          before writing to the bfq file)
     * @param includeNonPfReads whether to include non pf-reads
     */
    public SamToBfqWriter(final File samFile, final String outputPrefix, final boolean pairedReads,
                          String namePrefix, boolean includeNonPfReads) {
        this(samFile, null, outputPrefix, null, null, pairedReads, namePrefix, includeNonPfReads, true, null);
    }

    /**
     * Writes the binary fastq file(s) to the output directory
     */
    public void writeBfqFiles() {

        final SamReader reader = SamReaderFactory.makeDefault().open(samFile);
        final Iterator<SAMRecord> iterator = reader.iterator();

        // Filter out noise reads and reads that fail the quality filter
        final TagFilter tagFilter = new TagFilter(ReservedTagConstants.XN, 1);
        final FailsVendorReadQualityFilter qualityFilter = new FailsVendorReadQualityFilter();
        final WholeReadClippedFilter clippedFilter = new WholeReadClippedFilter();


        if (!pairedReads) {
            List<SamRecordFilter> filters = new ArrayList<>();
            filters.add(tagFilter);
            filters.add(clippedFilter);
            if (!this.includeNonPfReads) {
                filters.add(qualityFilter);
            }
            writeSingleEndBfqs(iterator, filters);
            codec1.close();
        } else {
            writePairedEndBfqs(iterator, tagFilter, qualityFilter, clippedFilter);
            codec1.close();
            codec2.close();
        }
        log.info("Wrote " + wrote + " bfq records.");
        CloserUtil.close(reader);
    }

    /**
     * Path for writing bfqs for paired-end reads
     *
     * @param iterator      the iterator witht he SAM Records to write
     * @param tagFilter     the filter for noise reads
     * @param qualityFilter the filter for PF reads
     */
    private void writePairedEndBfqs(final Iterator<SAMRecord> iterator, final TagFilter tagFilter,
                                    final FailsVendorReadQualityFilter qualityFilter,
                                    SamRecordFilter... otherFilters) {
        // Open the codecs for writing
        int fileIndex = 0;
        initializeNextBfqFiles(fileIndex++);

        int records = 0;

        RECORD_LOOP:
        while (iterator.hasNext()) {
            final SAMRecord first = iterator.next();
            if (!iterator.hasNext()) {
                throw new PicardException("Mismatched number of records in " + this.samFile.getAbsolutePath());
            }
            final SAMRecord second = iterator.next();
            if (!second.getReadName().equals(first.getReadName()) ||
                    first.getFirstOfPairFlag() == second.getFirstOfPairFlag()) {
                throw new PicardException("Unmatched read pairs in " + this.samFile.getAbsolutePath() +
                        ": " + first.getReadName() + ", " + second.getReadName() + ".");
            }

            // If *both* are noise reads, filter them out
            if (tagFilter.filterOut(first) && tagFilter.filterOut(second)) {
                continue;
            }

            // If either fails to pass filter, then exclude them as well
            if (!includeNonPfReads && (qualityFilter.filterOut(first) || qualityFilter.filterOut(second))) {
                continue;
            }

            // If either fails any of the other filters, exclude them both
            for (SamRecordFilter filter : otherFilters) {
                if (filter.filterOut(first) || filter.filterOut(second)) {
                    continue RECORD_LOOP;
                }
            }

            // Otherwise, write them out
            records++;
            if (records % increment == 0) {
                first.setReadName(first.getReadName() + "/1");
                writeFastqRecord(first.getFirstOfPairFlag() ? codec1 : codec2, first);
                second.setReadName(second.getReadName() + "/2");
                writeFastqRecord(second.getFirstOfPairFlag() ? codec1 : codec2, second);
                wrote++;
                if (wrote % 1000000 == 0) {
                    log.info(wrote + " records written.");
                }
                if (chunk > 0 && wrote % chunk == 0) {
                    initializeNextBfqFiles(fileIndex++);
                }
            }
        }
    }

    /**
     * Path for writing bfqs for single-end reads
     *
     * @param iterator the iterator with he SAM Records to write
     * @param filters  the list of filters to be applied
     */
    private void writeSingleEndBfqs(final Iterator<SAMRecord> iterator, final List<SamRecordFilter> filters) {

        // Open the codecs for writing
        int fileIndex = 0;
        initializeNextBfqFiles(fileIndex++);

        int records = 0;

        final FilteringSamIterator it = new FilteringSamIterator(iterator, new AggregateFilter(filters));
        while (it.hasNext()) {
            final SAMRecord record = it.next();
            records++;
            if (records % increment == 0) {

                record.setReadName(record.getReadName() + "/1");
                writeFastqRecord(codec1, record);
                wrote++;
                if (wrote % 1000000 == 0) {
                    log.info(wrote + " records processed.");
                }
                if (chunk > 0 && wrote % chunk == 0) {
                    initializeNextBfqFiles(fileIndex++);
                }
            }
        }
    }

    /**
     * Closes any open bfq file(s), and opens the new one(s)
     *
     * @param fileIndex the index (counter) of the files to write
     */
    private void initializeNextBfqFiles(final int fileIndex) {
        // Close the codecs if they were writing before
        if (codec1 != null) {
            codec1.close();
            if (pairedReads) {
                codec2.close();
            }
        }

        // Open new file, using the fileIndex.
        final File bfq1 = getOutputFile(this.outputPrefix, 1, fileIndex);
        codec1 = new BinaryCodec(IOUtil.openFileForWriting(bfq1));
        log.info("Now writing to file " + bfq1.getAbsolutePath());
        if (pairedReads) {
            final File bfq2 = getOutputFile(this.outputPrefix, 2, fileIndex);
            codec2 = new BinaryCodec(IOUtil.openFileForWriting(bfq2));
            log.info("Now writing to file " + bfq2.getAbsolutePath());
        }
    }

    /**
     * Writes out a SAMRecord in Maq fastq format
     *
     * @param codec the code to write to
     * @param rec   the SAMRecord to write
     */
    private void writeFastqRecord(final BinaryCodec codec, final SAMRecord rec) {

        // Trim the run barcode off the read name
        String readName = rec.getReadName();
        if (namePrefix != null && readName.startsWith(namePrefix)) {
            readName = readName.substring(nameTrim);
        }
        // Writes the length of the read name and then the name (null-terminated)
        codec.writeString(readName, true, true);

        final char[] seqs = rec.getReadString().toCharArray();
        final char[] quals = rec.getBaseQualityString().toCharArray();

        int retainedLength = seqs.length;
        if (clipAdapters) {
            // adjust to a shorter length iff clipping tag exists
            Integer trimPoint = rec.getIntegerAttribute(ReservedTagConstants.XT);
            if (trimPoint != null) {
                ValidationUtils.validateArg(rec.getReadLength() == seqs.length, () -> "length of read and seqs differ. Found " + rec.getReadLength() + " and '" + seqs.length + ".");

                retainedLength = Math.min(seqs.length, Math.max(SEED_REGION_LENGTH, trimPoint - 1));
            }
        }

        // Write the length of the sequence
        codec.writeInt(basesToWrite != null ? basesToWrite : seqs.length);

        // Calculate and write the sequence and qualities
        final byte[] seqsAndQuals = encodeSeqsAndQuals(seqs, quals, retainedLength);
        codec.writeBytes(seqsAndQuals);
    }

    private byte[] encodeSeqsAndQuals(char[] seqs, char[] quals, int retainedLength) {
        final byte[] seqsAndQuals = new byte[basesToWrite == null ? seqs.length : basesToWrite];

        int seedRegionNoCallFixes = 0;
        for (int i = 0; i < retainedLength && i < seqsAndQuals.length; i++) {
            int quality = Math.min(quals[i] - 33, 63);
            final int base;
            switch (seqs[i]) {
                case 'A':
                case 'a':
                    base = 0;
                    break;
                case 'C':
                case 'c':
                    base = 1;
                    break;
                case 'G':
                case 'g':
                    base = 2;
                    break;
                case 'T':
                case 't':
                    base = 3;
                    break;
                case 'N':
                case 'n':
                case '.':
                    base = 0;
                    if (i < SEED_REGION_LENGTH) {
                        if (seedRegionNoCallFixes < MAX_SEED_REGION_NOCALL_FIXES) {
                            quality = 1;
                            seedRegionNoCallFixes++;
                        } else {
                            quality = 0;
                        }
                    } else {
                        quality = 1;
                    }
                    break;
                default:
                    throw new PicardException("Unknown base when writing bfq file: " + seqs[i]);
            }
            seqsAndQuals[i] = encodeBaseAndQuality(base, quality);
        }
        // rewrite clipped adapter with all A's of quality 1
        for (int i = retainedLength; i < seqsAndQuals.length; i++) {
            seqsAndQuals[i] = encodeBaseAndQuality(0, 1);
        }

        return seqsAndQuals;
    }

    private byte encodeBaseAndQuality(int base, int quality) {
        return (byte) ((base << 6) | quality);
    }

    /**
     * Count the number of records in the bamFile that could potentially be written
     *
     * @return the number of records in the Bam file
     */
    private int countWritableRecords() {
        int count = 0;

        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(referenceSequence).open(samFile);
        if (!reader.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
            //this is a fix for issue PIC-274: It looks like BamToBfqWriter requires that the input BAM is queryname sorted,
            //but it doesn't check this early, nor produce an understandable error message."
            throw new PicardException("Input file (" + samFile.getAbsolutePath() + ") needs to be sorted by queryname.");
        }
        final PeekableIterator<SAMRecord> it = new PeekableIterator<SAMRecord>(reader.iterator());
        if (!this.pairedReads) {
            // Filter out noise reads and reads that fail the quality filter
            final List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
            filters.add(new TagFilter(ReservedTagConstants.XN, 1));
            if (!this.includeNonPfReads) {
                filters.add(new FailsVendorReadQualityFilter());
            }
            final FilteringSamIterator itr = new FilteringSamIterator(it, new AggregateFilter(filters));
            while (itr.hasNext()) {
                itr.next();
                count++;
            }
        } else {
            while (it.hasNext()) {
                final SAMRecord first = it.next();
                final SAMRecord second = it.next();
                // If both are noise reads, filter them out
                if (first.getAttribute(ReservedTagConstants.XN) != null &&
                        second.getAttribute(ReservedTagConstants.XN) != null) {
                    // skip it
                }
                // If either fails to pass filter, then exclude them as well
                else if (!this.includeNonPfReads && (first.getReadFailsVendorQualityCheckFlag() || second.getReadFailsVendorQualityCheckFlag())) {
                    // skip it
                }
                // Otherwise, write them out
                else {
                    count++;
                }
            }
        }
        it.close();
        CloserUtil.close(reader);
        return count;
    }

    /**
     * Constructs the name for the output file and returns the file
     *
     * @param outputPrefix the directory and file prefix for the output bfq file
     * @param read         whether this is the file for the first or second read
     * @param index        used in file name
     * @return a new File object for the bfq file.
     */
    private File getOutputFile(final String outputPrefix, final int read, final int index) {
        final File result = new File(outputPrefix + index + "." + read + ".bfq");
        IOUtil.assertFileIsWritable(result);
        return result;
    }
}
