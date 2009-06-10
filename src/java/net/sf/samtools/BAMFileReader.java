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
package net.sf.samtools;


import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.StringLineReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Internal class for reading and querying BAM files.
 */
class BAMFileReader
    extends SAMFileReader.ReaderImplementation {
    // True if reading from a File rather than an InputStream
    private boolean mIsSeekable = false;
    // For converting bytes into other primitive types
    private BinaryCodec mStream = null;
    // Underlying compressed data stream.
    private final BlockCompressedInputStream mCompressedInputStream;
    private SAMFileHeader mFileHeader = null;
    // Populated if the file is seekable and an index exists
    private BAMFileIndex mFileIndex = null;
    private long mFirstRecordPointer = 0;
    private CloseableIterator<SAMRecord> mCurrentIterator = null;
    // If true, all SAMRecords are fully decoded as they are read.
    private final boolean eagerDecode;
    // For error-checking.
    private ValidationStringency mValidationStringency = SAMFileReader.ValidationStringency.SILENT;

    /**
     * Prepare to read BAM from a stream (not seekable)
     * @param stream source of bytes.
     * @param eagerDecode if true, decode all BAM fields as reading rather than lazily.
     */
    BAMFileReader(final InputStream stream, final boolean eagerDecode)
        throws IOException {
        mIsSeekable = false;
        mCompressedInputStream = new BlockCompressedInputStream(stream);
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        this.eagerDecode = eagerDecode;
        readHeader(null);
    }

    /**
     * Prepare to read BAM from a file (seekable)
     * @param file source of bytes.
     * @param eagerDecode if true, decode all BAM fields as reading rather than lazily.
     */
    BAMFileReader(final File file, final boolean eagerDecode)
        throws IOException {
        mIsSeekable = true;
        mCompressedInputStream = new BlockCompressedInputStream(file);
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        this.eagerDecode = eagerDecode;
        readHeader(file);
        mFirstRecordPointer = mCompressedInputStream.getFilePointer();
    }

    void close() {
        if (mStream != null) {
            mStream.close();
        }
        mStream = null;
        mFileHeader = null;
        mFileIndex = null;
    }

    /**
     * @return the file index, if one exists, else null.
     */
    BAMFileIndex getFileIndex() {
        return mFileIndex;
    }

    void setFileIndex(final BAMFileIndex fileIndex) {
        mFileIndex = fileIndex;
    }

    SAMFileHeader getFileHeader() {
        return mFileHeader;
    }

    /**
     * error-checking level for subsequent SAMRecord reads.
     */
    void setValidationStringency(final SAMFileReader.ValidationStringency validationStringency) {
        this.mValidationStringency = validationStringency;
    }

    SAMFileReader.ValidationStringency getValidationStringency() {
        return this.mValidationStringency;
    }

    /**
     * Prepare to iterate through the SAMRecords in file order.
     * Only a single iterator on a BAM file can be extant at a time.  If getIterator() or a query method has been called once,
     * that iterator must be closed before getIterator() can be called again.
     * A somewhat peculiar aspect of this method is that if the file is not seekable, a second call to
     * getIterator() begins its iteration where the last one left off.  That is the best that can be
     * done in that situation.
     */
    CloseableIterator<SAMRecord> getIterator() {
        if (mStream == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        if (mIsSeekable) {
            try {
                mCompressedInputStream.seek(mFirstRecordPointer);
            } catch (IOException exc) {
                throw new RuntimeException(exc.getMessage(), exc);
            }
        }
        mCurrentIterator = new BAMFileIterator();
        return mCurrentIterator;
    }

    /**
     * Prepare to iterate through the SAMRecords that match the given interval.
     * Only a single iterator on a BAMFile can be extant at a time.  The previous one must be closed
     * before calling any of the methods that return an iterator.
     *
     * Note that an unmapped SAMRecord may still have a reference name and an alignment start for sorting
     * purposes (typically this is the coordinate of its mate), and will be found by this method if the coordinate
     * matches the specified interval.
     *
     * Note that this method is not necessarily efficient in terms of disk I/O.  The index does not have perfect
     * resolution, so some SAMRecords may be read and then discarded because they do not match the specified interval.
     *
     * @param sequence Reference sequence sought.
     * @param start Desired SAMRecords must overlap or be contained in the interval specified by start and end.
     * A value of zero implies the start of the reference sequence.
     * @param end A value of zero implies the end of the reference sequence.
     * @param contained If true, the alignments for the SAMRecords must be completely contained in the interval
     * specified by start and end.  If false, the SAMRecords need only overlap the interval.
     * @return Iterator for the matching SAMRecords
     */
    CloseableIterator<SAMRecord> query(final String sequence, final int start, final int end, final boolean contained) {
        if (mStream == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        if (!mIsSeekable) {
            throw new UnsupportedOperationException("Cannot query stream-based BAM file");
        }
        if (mFileIndex == null) {
            throw new IllegalStateException("No BAM file index is available");
        }
        mCurrentIterator = new BAMFileIndexIterator(sequence, start, end, contained);
        return mCurrentIterator;
    }

    public CloseableIterator<SAMRecord> queryUnmapped() {
        if (mStream == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        if (!mIsSeekable) {
            throw new UnsupportedOperationException("Cannot query stream-based BAM file");
        }
        if (mFileIndex == null) {
            throw new IllegalStateException("No BAM file index is available");
        }
        try {
            final long startOfLastLinearBin = mFileIndex.getStartOfLastLinearBin();
            if (startOfLastLinearBin != -1) {
                mCompressedInputStream.seek(startOfLastLinearBin);
            } else {
                // No mapped reads in file, just start at the first read in file.
                mCompressedInputStream.seek(mFirstRecordPointer);
            }
            mCurrentIterator = new BAMFileIndexUnmappedIterator();
            return mCurrentIterator;
        } catch (IOException e) {
            throw new RuntimeException("IOException seeking to unmapped reads", e);
        }
    }

    /**
     * Reads the header from the file or stream
     * @param file Note that this is used only for reporting errors.
     */
    private void readHeader(final File file)
        throws IOException {

        final byte[] buffer = new byte[4];
        mStream.readBytes(buffer);
        if (!Arrays.equals(buffer, BAMFileConstants.BAM_MAGIC)) {
            throw new IOException("Invalid BAM file header");
        }

        final int headerTextLength = mStream.readInt();
        final String textHeader = mStream.readString(headerTextLength);
        final SAMTextHeaderCodec headerCodec = new SAMTextHeaderCodec();
        headerCodec.setValidationStringency(mValidationStringency);
        mFileHeader = headerCodec.decode(new StringLineReader(textHeader),
                file);

        final int sequenceCount = mStream.readInt();
        if (mFileHeader.getSequenceDictionary().size() > 0) {
            // It is allowed to have binary sequences but no text sequences, so only validate if both are present
            if (sequenceCount != mFileHeader.getSequenceDictionary().size()) {
                throw new SAMFormatException("Number of sequences in text header (" +
                        mFileHeader.getSequenceDictionary().size() +
                        ") != number of sequences in binary header (" + sequenceCount + ") for file " + file);
            }
            for (int i = 0; i < sequenceCount; i++) {
                final SAMSequenceRecord binarySequenceRecord = readSequenceRecord(file);
                final SAMSequenceRecord sequenceRecord = mFileHeader.getSequence(i);
                if (!sequenceRecord.getSequenceName().equals(binarySequenceRecord.getSequenceName())) {
                    throw new SAMFormatException("For sequence " + i + ", text and binary have different names in file " +
                            file);
                }
                if (sequenceRecord.getSequenceLength() != binarySequenceRecord.getSequenceLength()) {
                    throw new SAMFormatException("For sequence " + i + ", text and binary have different lengths in file " +
                            file);
                }
            }
        } else {
            // If only binary sequences are present, copy them into mFileHeader
            final List<SAMSequenceRecord> sequences = new ArrayList<SAMSequenceRecord>(sequenceCount);
            for (int i = 0; i < sequenceCount; i++) {
                sequences.add(readSequenceRecord(file));
            }
            mFileHeader.setSequenceDictionary(new SAMSequenceDictionary(sequences));
        }
    }

    /**
     * Reads a single binary sequence record from the file or stream
     * @param file Note that this is used only for reporting errors.
     */
    private SAMSequenceRecord readSequenceRecord(final File file) {
        final int nameLength = mStream.readInt();
        if (nameLength <= 1) {
            throw new SAMFormatException("Invalid BAM file header: missing sequence name in file " + file);
        }
        final String sequenceName = mStream.readString(nameLength - 1);
        // Skip the null terminator
        mStream.readByte();
        final int sequenceLength = mStream.readInt();
        return new SAMSequenceRecord(sequenceName, sequenceLength);
    }

    /**
     * Iterator for non-indexed sequential iteration through all SAMRecords in file.
     * Starting point of iteration is wherever current file position is when the iterator is constructed.
     */
    private class BAMFileIterator implements CloseableIterator<SAMRecord> {
        private SAMRecord mNextRecord = null;
        private final BAMRecordCodec bamRecordCodec = new BAMRecordCodec(getFileHeader());
        private long samRecordIndex = 0; // Records at what position (counted in records) we are at in the file

        BAMFileIterator() {
            this(true);
        }

        /**
         * @param advance Trick to enable subclass to do more setup before advancing
         */
        BAMFileIterator(final boolean advance) {
            this.bamRecordCodec.setInputStream(BAMFileReader.this.mStream.getInputStream());

            if (advance) {
                advance();
            }
        }

        public void close() {
            if (this != mCurrentIterator) {
                throw new IllegalStateException("Attempt to close non-current iterator");
            }
            mCurrentIterator = null;
        }

        public boolean hasNext() {
            return (mNextRecord != null);
        }

        public SAMRecord next() {
            final SAMRecord result = mNextRecord;
            advance();
            return result;
        }

        public void remove() {
            throw new UnsupportedOperationException("Not supported: remove");
        }

        void advance() {
            try {
                mNextRecord = getNextRecord();
                if (mNextRecord != null) {
                    ++this.samRecordIndex;
                    // Because some decoding is done lazily, the record needs to remember the validation stringency.
                    mNextRecord.setValidationStringency(mValidationStringency);

                    if (mValidationStringency != ValidationStringency.SILENT) {
                        final List<SAMValidationError> validationErrors = mNextRecord.isValid();
                        SAMUtils.processValidationErrors(validationErrors,
                                this.samRecordIndex, BAMFileReader.this.getValidationStringency());
                    }
                }
                if (eagerDecode && mNextRecord != null) {
                    mNextRecord.eagerDecode();
                }
            } catch (IOException exc) {
                throw new RuntimeException(exc.getMessage(), exc);
            }
        }

        /**
         * Read the next record from the input stream.
         */
        SAMRecord getNextRecord() throws IOException {
            return bamRecordCodec.decode();
        }

        /**
         * @return The record that will be return by the next call to next()
         */
        protected SAMRecord peek() {
            return mNextRecord;
        }
    }

    private class BAMFileIndexIterator
        extends BAMFileIterator {

        private long[] mFilePointers = null;
        private int mFilePointerIndex = 0;
        private long mFilePointerLimit = -1;
        private int mReferenceIndex = -1;
        private int mRegionStart = 0;
        private int mRegionEnd = 0;
        private boolean mReturnContained = false;


        /**
         * Prepare to iterate through SAMRecords matching the target interval.
         * @param sequence Desired reference sequence.
         * @param start 1-based start of target interval, inclusive.
         * @param end 1-based end of target interval, inclusive.
         * @param contained If true, SAMRecord must be contained in the target interval.  If false, SAMRecord need
         * only overlap the target interval.
         */
        BAMFileIndexIterator(final String sequence, final int start, final int end, final boolean contained) {
            super(false);  // delay advance() until after construction
            final SAMFileHeader fileHeader = getFileHeader();
            mReferenceIndex = fileHeader.getSequenceIndex(sequence);
            if (mReferenceIndex != -1) {
                final BAMFileIndex fileIndex = getFileIndex();
                mFilePointers = fileIndex.getSearchBins(mReferenceIndex, start, end);
            }
            mRegionStart = start;
            mRegionEnd = (end <= 0) ? Integer.MAX_VALUE : end;
            mReturnContained = contained;
            advance();
        }

        SAMRecord getNextRecord()
            throws IOException {
            while (true) {
                // Advance to next file block if necessary
                while (mCompressedInputStream.getFilePointer() >= mFilePointerLimit) {
                    if (mFilePointers == null ||
                        mFilePointerIndex >= mFilePointers.length) {
                        return null;
                    }
                    final long startOffset = mFilePointers[mFilePointerIndex++];
                    final long endOffset = mFilePointers[mFilePointerIndex++];
                    mCompressedInputStream.seek(startOffset);
                    mFilePointerLimit = endOffset;
                }
                // Pull next record from stream
                final SAMRecord record = super.getNextRecord();
                if (record == null) {
                    return null;
                }
                // If beyond the end of this reference sequence, end iteration
                final int referenceIndex = record.getReferenceIndex();
                if (referenceIndex != mReferenceIndex) {
                    if (referenceIndex < 0 ||
                        referenceIndex > mReferenceIndex) {
                        mFilePointers = null;
                        return null;
                    }
                    // If before this reference sequence, continue
                    continue;
                }
                if (mRegionStart == 0 && mRegionEnd == Integer.MAX_VALUE) {
                    // Quick exit to avoid expensive alignment end calculation
                    return record;
                }
                final int alignmentStart = record.getAlignmentStart();
                final int alignmentEnd = record.getAlignmentEnd();
                if (alignmentStart > mRegionEnd) {
                    // If scanned beyond target region, end iteration
                    mFilePointers = null;
                    return null;
                }
                // Filter for overlap with region
                if (mReturnContained) {
                    if (alignmentStart >= mRegionStart && alignmentEnd <= mRegionEnd) {
                        return record;
                    }
                } else {
                    if (alignmentEnd >= mRegionStart && alignmentStart <= mRegionEnd) {
                        return record;
                    }
                }
            }
        }
    }

    private class BAMFileIndexUnmappedIterator extends BAMFileIterator  {
        private BAMFileIndexUnmappedIterator() {
            while (this.hasNext() && peek().getReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                advance();
            }
        }
    }

}
