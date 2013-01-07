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


import net.sf.samtools.util.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.seekablestream.SeekableStream;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Internal class for reading and querying BAM files.
 */
class BAMFileReader extends SAMFileReader.ReaderImplementation {
    // True if reading from a File rather than an InputStream
    private boolean mIsSeekable = false;

    // For converting bytes into other primitive types
    private BinaryCodec mStream = null;

    // Underlying compressed data stream.
    private final BlockCompressedInputStream mCompressedInputStream;
    private SAMFileHeader mFileHeader = null;

    // One of these is populated if the file is seekable and an index exists
    private File mIndexFile = null;
    private SeekableStream mIndexStream = null;

    private BAMIndex mIndex = null;
    private long mFirstRecordPointer = 0;
    private CloseableIterator<SAMRecord> mCurrentIterator = null;

    // If true, all SAMRecords are fully decoded as they are read.
    private final boolean eagerDecode;

    // For error-checking.
    private ValidationStringency mValidationStringency;

    // For creating BAMRecords
    private SAMRecordFactory samRecordFactory;

    /**
     * Use the caching index reader implementation rather than the disk-hit-per-file model.
     */
    private boolean mEnableIndexCaching = false;

    /**
     * Use the traditional memory-mapped implementation for BAM file indexes rather than regular I/O.
     */
    private boolean mEnableIndexMemoryMapping = true;

    /**
     * Add information about the origin (reader and position) to SAM records.
     */
    private SAMFileReader mFileReader = null;

    /**
     * Prepare to read BAM from a stream (not seekable)
     * @param stream source of bytes.
     * @param eagerDecode if true, decode all BAM fields as reading rather than lazily.
     * @param validationStringency Controls how to handle invalidate reads or header lines.
     */
    BAMFileReader(final InputStream stream,
                  final File indexFile,
                  final boolean eagerDecode,
                  final ValidationStringency validationStringency,
                  final SAMRecordFactory factory)
        throws IOException {
        mIndexFile = indexFile;
        mIsSeekable = false;
        mCompressedInputStream = new BlockCompressedInputStream(stream);
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        this.eagerDecode = eagerDecode;
        this.mValidationStringency = validationStringency;
        this.samRecordFactory = factory;
        readHeader(null);
    }

    /**
     * Prepare to read BAM from a file (seekable)
     * @param file source of bytes.
     * @param eagerDecode if true, decode all BAM fields as reading rather than lazily.
     * @param validationStringency Controls how to handle invalidate reads or header lines.
     */
    BAMFileReader(final File file,
                  final File indexFile,
                  final boolean eagerDecode,
                  final ValidationStringency validationStringency,
                  final SAMRecordFactory factory)
        throws IOException {
        this(new BlockCompressedInputStream(file), indexFile!=null ? indexFile : findIndexFile(file), eagerDecode, file.getAbsolutePath(), validationStringency, factory);
        if (mIndexFile != null && mIndexFile.lastModified() < file.lastModified()) {
            System.err.println("WARNING: BAM index file " + mIndexFile.getAbsolutePath() +
                    " is older than BAM " + file.getAbsolutePath());
        }
        // Provide better error message when there is an error reading.
        mStream.setInputFileName(file.getAbsolutePath());
    }

    BAMFileReader(final SeekableStream strm,
                  final File indexFile,
                  final boolean eagerDecode,
                  final ValidationStringency validationStringency,
                  final SAMRecordFactory factory)
        throws IOException {
        this(new BlockCompressedInputStream(strm), indexFile, eagerDecode, strm.getSource(), validationStringency, factory);
    }

    BAMFileReader(final SeekableStream strm,
                  final SeekableStream indexStream,
                  final boolean eagerDecode,
                  final ValidationStringency validationStringency,
                  final SAMRecordFactory factory)
        throws IOException {
        this(new BlockCompressedInputStream(strm), indexStream, eagerDecode, strm.getSource(), validationStringency, factory);
    }

    private BAMFileReader(final BlockCompressedInputStream compressedInputStream,
                          final File indexFile,
                          final boolean eagerDecode,
                          final String source,
                          final ValidationStringency validationStringency,
                          final SAMRecordFactory factory)
        throws IOException {
        mIndexFile = indexFile;
        mIsSeekable = true;
        mCompressedInputStream = compressedInputStream;
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        this.eagerDecode = eagerDecode;
        this.mValidationStringency = validationStringency;
        this.samRecordFactory = factory;
        readHeader(source);
        mFirstRecordPointer = mCompressedInputStream.getFilePointer();
    }    

    private BAMFileReader(final BlockCompressedInputStream compressedInputStream,
                          final SeekableStream indexStream,
                          final boolean eagerDecode,
                          final String source,
                          final ValidationStringency validationStringency,
                          final SAMRecordFactory factory)
        throws IOException {
        mIndexStream = indexStream;
        mIsSeekable = true;
        mCompressedInputStream = compressedInputStream;
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        this.eagerDecode = eagerDecode;
        this.mValidationStringency = validationStringency;
        this.samRecordFactory = factory;
        readHeader(source);
        mFirstRecordPointer = mCompressedInputStream.getFilePointer();
    }

    /**
     * If true, writes the source of every read into the source SAMRecords.
     * @param enabled true to write source information into each SAMRecord.
     */
    void enableFileSource(final SAMFileReader reader, final boolean enabled) {
        this.mFileReader = enabled ? reader : null;
    }

    /**
     * If true, uses the caching version of the index reader.
     * @param enabled true to write source information into each SAMRecord.
     */
    public void enableIndexCaching(final boolean enabled) {
        if(mIndex != null)
            throw new SAMException("Unable to turn on index caching; index file has already been loaded.");
        this.mEnableIndexCaching = enabled;
    }

    /**
     * If false, disable the use of memory mapping for accessing index files (default behavior is to use memory mapping).
     * This is slower but more scalable when accessing large numbers of BAM files sequentially.
     * @param enabled True to use memory mapping, false to use regular I/O.
     */
    public void enableIndexMemoryMapping(final boolean enabled) {
        if (mIndex != null) {
            throw new SAMException("Unable to change index memory mapping; index file has already been loaded.");
        }
        this.mEnableIndexMemoryMapping = enabled;
    }

    @Override void enableCrcChecking(final boolean enabled) {
        this.mCompressedInputStream.setCheckCrcs(enabled);
    }

    @Override void setSAMRecordFactory(final SAMRecordFactory factory) { this.samRecordFactory = factory; }

    /**
     * @return true if ths is a BAM file, and has an index
     */
    public boolean hasIndex() {
        return (mIndexFile != null) || (mIndexStream != null);
    }

    /**
     * Retrieves the index for the given file type.  Ensure that the index is of the specified type.
     * @return An index of the given type.
     */
    public BAMIndex getIndex() {
        if(!hasIndex())
            throw new SAMException("No index is available for this BAM file.");
        if(mIndex == null) {
            if (mIndexFile != null)
                mIndex = mEnableIndexCaching ? new CachingBAMFileIndex(mIndexFile, getFileHeader().getSequenceDictionary(), mEnableIndexMemoryMapping)
                                             : new DiskBasedBAMFileIndex(mIndexFile, getFileHeader().getSequenceDictionary(), mEnableIndexMemoryMapping);
            else
                mIndex = mEnableIndexCaching ? new CachingBAMFileIndex(mIndexStream, getFileHeader().getSequenceDictionary())
                                             : new DiskBasedBAMFileIndex(mIndexStream, getFileHeader().getSequenceDictionary());
        }
        return mIndex;
    }

    void close() {
        if (mStream != null) {
            mStream.close();
        }
        if (mIndex != null) {
            mIndex.close();
        }
        mStream = null;
        mFileHeader = null;
        mIndex = null;
    }

    SAMFileHeader getFileHeader() {
        return mFileHeader;
    }

    /**
     * Set error-checking level for subsequent SAMRecord reads.
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

    @Override
    CloseableIterator<SAMRecord> getIterator(final SAMFileSpan chunks) {
        if (mStream == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        if (!(chunks instanceof BAMFileSpan)) {
            throw new IllegalStateException("BAMFileReader cannot handle this type of file span.");
        }

        // Create an iterator over the given chunk boundaries.
        mCurrentIterator = new BAMFileIndexIterator(((BAMFileSpan)chunks).toCoordinateArray());
        return mCurrentIterator;
    }

    /**
     * Gets an unbounded pointer to the first record in the BAM file.  Because the reader doesn't necessarily know
     * when the file ends, the rightmost bound of the file pointer will not end exactly where the file ends.  However,
     * the rightmost bound is guaranteed to be after the last read in the file.
     * @return An unbounded pointer to the first record in the BAM file.
     */
    @Override
    SAMFileSpan getFilePointerSpanningReads() {
        return new BAMFileSpan(new Chunk(mFirstRecordPointer,Long.MAX_VALUE));
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
        mCurrentIterator = createIndexIterator(sequence, start, end, contained? QueryType.CONTAINED: QueryType.OVERLAPPING);
        return mCurrentIterator;
    }

    /**
     * Prepare to iterate through the SAMRecords with the given alignment start.
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
     * @param start Alignment start sought.
     * @return Iterator for the matching SAMRecords.
     */
    CloseableIterator<SAMRecord> queryAlignmentStart(final String sequence, final int start) {
        if (mStream == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        if (!mIsSeekable) {
            throw new UnsupportedOperationException("Cannot query stream-based BAM file");
        }
        mCurrentIterator = createIndexIterator(sequence, start, -1, QueryType.STARTING_AT);
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
        try {
            final long startOfLastLinearBin = getIndex().getStartOfLastLinearBin();
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
     * @param source Note that this is used only for reporting errors.
     */
    private void readHeader(final String source)
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
                source);

        final int sequenceCount = mStream.readInt();
        if (mFileHeader.getSequenceDictionary().size() > 0) {
            // It is allowed to have binary sequences but no text sequences, so only validate if both are present
            if (sequenceCount != mFileHeader.getSequenceDictionary().size()) {
                throw new SAMFormatException("Number of sequences in text header (" +
                        mFileHeader.getSequenceDictionary().size() +
                        ") != number of sequences in binary header (" + sequenceCount + ") for file " + source);
            }
            for (int i = 0; i < sequenceCount; i++) {
                final SAMSequenceRecord binarySequenceRecord = readSequenceRecord(source);
                final SAMSequenceRecord sequenceRecord = mFileHeader.getSequence(i);
                if (!sequenceRecord.getSequenceName().equals(binarySequenceRecord.getSequenceName())) {
                    throw new SAMFormatException("For sequence " + i + ", text and binary have different names in file " +
                            source);
                }
                if (sequenceRecord.getSequenceLength() != binarySequenceRecord.getSequenceLength()) {
                    throw new SAMFormatException("For sequence " + i + ", text and binary have different lengths in file " +
                            source);
                }
            }
        } else {
            // If only binary sequences are present, copy them into mFileHeader
            final List<SAMSequenceRecord> sequences = new ArrayList<SAMSequenceRecord>(sequenceCount);
            for (int i = 0; i < sequenceCount; i++) {
                sequences.add(readSequenceRecord(source));
            }
            mFileHeader.setSequenceDictionary(new SAMSequenceDictionary(sequences));
        }
    }

    /**
     * Reads a single binary sequence record from the file or stream
     * @param source Note that this is used only for reporting errors.
     */
    private SAMSequenceRecord readSequenceRecord(final String source) {
        final int nameLength = mStream.readInt();
        if (nameLength <= 1) {
            throw new SAMFormatException("Invalid BAM file header: missing sequence name in file " + source);
        }
        final String sequenceName = mStream.readString(nameLength - 1);
        // Skip the null terminator
        mStream.readByte();
        final int sequenceLength = mStream.readInt();
        return new SAMSequenceRecord(SAMSequenceRecord.truncateSequenceName(sequenceName), sequenceLength);
    }

    /**
     * Iterator for non-indexed sequential iteration through all SAMRecords in file.
     * Starting point of iteration is wherever current file position is when the iterator is constructed.
     */
    private class BAMFileIterator implements CloseableIterator<SAMRecord> {
        private SAMRecord mNextRecord = null;
        private final BAMRecordCodec bamRecordCodec;
        private long samRecordIndex = 0; // Records at what position (counted in records) we are at in the file
        private boolean isClosed = false;

        BAMFileIterator() {
            this(true);
        }

        /**
         * @param advance Trick to enable subclass to do more setup before advancing
         */
        BAMFileIterator(final boolean advance) {
            this.bamRecordCodec = new BAMRecordCodec(getFileHeader(), samRecordFactory);
            this.bamRecordCodec.setInputStream(BAMFileReader.this.mStream.getInputStream(),
                    BAMFileReader.this.mStream.getInputFileName());

            if (advance) {
                advance();
            }
        }

        public void close() {
            if (!isClosed) {
                if (mCurrentIterator != null && this != mCurrentIterator) {
                    throw new IllegalStateException("Attempt to close non-current iterator");
                }
                mCurrentIterator = null;
                isClosed = true;
            }
        }

        public boolean hasNext() {
            if (isClosed) throw new IllegalStateException("Iterator has been closed");
            return (mNextRecord != null);
        }

        public SAMRecord next() {
            if (isClosed) throw new IllegalStateException("Iterator has been closed");
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
            final long startCoordinate = mCompressedInputStream.getFilePointer();
            final SAMRecord next = bamRecordCodec.decode();
            final long stopCoordinate = mCompressedInputStream.getFilePointer();

            if(mFileReader != null && next != null)
                next.setFileSource(new SAMFileSource(mFileReader,new BAMFileSpan(new Chunk(startCoordinate,stopCoordinate))));

            return next;
        }

        /**
         * @return The record that will be return by the next call to next()
         */
        protected SAMRecord peek() {
            return mNextRecord;
        }
    }

    /**
     * Prepare to iterate through SAMRecords matching the target interval.
     * @param sequence Desired reference sequence.
     * @param start 1-based start of target interval, inclusive.
     * @param end 1-based end of target interval, inclusive.
     * @param queryType contained, overlapping, or starting-at query.
     */
    private CloseableIterator<SAMRecord> createIndexIterator(final String sequence,
                                                             final int start,
                                                             final int end,
                                                             final QueryType queryType) {
        long[] filePointers = null;

        // Hit the index to determine the chunk boundaries for the required data.
        final SAMFileHeader fileHeader = getFileHeader();
        final int referenceIndex = fileHeader.getSequenceIndex(sequence);
        if (referenceIndex != -1) {
            final BAMIndex fileIndex = getIndex();
            final BAMFileSpan fileSpan = fileIndex.getSpanOverlapping(referenceIndex, start, end);
            filePointers = fileSpan != null ? fileSpan.toCoordinateArray() : null;
        }

        // Create an iterator over the above chunk boundaries.
        final BAMFileIndexIterator iterator = new BAMFileIndexIterator(filePointers);

        // Add some preprocessing filters for edge-case reads that don't fit into this
        // query type.
        return new BAMQueryFilteringIterator(iterator,sequence,start,end,queryType);
    }

    enum QueryType {CONTAINED, OVERLAPPING, STARTING_AT}

    /**
     * Look for BAM index file according to standard naming convention.
     *
     * @param dataFile BAM file name.
     * @return Index file name, or null if not found.
     */
    private static File findIndexFile(final File dataFile) {
        // If input is foo.bam, look for foo.bai
        final String bamExtension = ".bam";
        File indexFile;
        final String fileName = dataFile.getName();
        if (fileName.endsWith(bamExtension)) {
            final String bai = fileName.substring(0, fileName.length() - bamExtension.length()) + BAMIndex.BAMIndexSuffix;
            indexFile = new File(dataFile.getParent(), bai);
            if (indexFile.exists()) {
                return indexFile;
            }
        }

        // If foo.bai doesn't exist look for foo.bam.bai
        indexFile = new File(dataFile.getParent(), dataFile.getName() + ".bai");
        if (indexFile.exists()) {
            return indexFile;
        } else {
            return null;
        }
    }    

    private class BAMFileIndexIterator extends BAMFileIterator {

        private long[] mFilePointers = null;
        private int mFilePointerIndex = 0;
        private long mFilePointerLimit = -1;

        /**
         * Prepare to iterate through SAMRecords stored in the specified compressed blocks at the given offset.
         * @param filePointers the block / offset combination, stored in chunk format.
         */
        BAMFileIndexIterator(final long[] filePointers) {
            super(false);  // delay advance() until after construction
            mFilePointers = filePointers;
            advance();
        }

        SAMRecord getNextRecord()
            throws IOException {
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
            return super.getNextRecord();
        }
    }

    /**
     * A decorating iterator that filters out records that are outside the bounds of the
     * given query parameters.
     */
    private class BAMQueryFilteringIterator implements CloseableIterator<SAMRecord> {
        /**
         * The wrapped iterator.
         */
        private final CloseableIterator<SAMRecord> wrappedIterator;

        /**
         * The next record to be returned.  Will be null if no such record exists.
         */
        private SAMRecord mNextRecord;

        private final int mReferenceIndex;
        private final int mRegionStart;
        private final int mRegionEnd;
        private final QueryType mQueryType;
        private boolean isClosed = false;

        public BAMQueryFilteringIterator(final CloseableIterator<SAMRecord> iterator,final String sequence, final int start, final int end, final QueryType queryType) {
            this.wrappedIterator = iterator;
            final SAMFileHeader fileHeader = getFileHeader();
            mReferenceIndex = fileHeader.getSequenceIndex(sequence);
            mRegionStart = start;
            if (queryType == QueryType.STARTING_AT) {
                mRegionEnd = mRegionStart;
            } else {
                mRegionEnd = (end <= 0) ? Integer.MAX_VALUE : end;
            }
            mQueryType = queryType;
            mNextRecord = advance();
        }

        /**
         * Returns true if a next element exists; false otherwise.
         */
        public boolean hasNext() {
            if (isClosed) throw new IllegalStateException("Iterator has been closed");
            return mNextRecord != null;
        }

        /**
         * Gets the next record from the given iterator.
         * @return The next SAM record in the iterator.
         */
        public SAMRecord next() {
            if(!hasNext())
                throw new NoSuchElementException("BAMQueryFilteringIterator: no next element available");
            final SAMRecord currentRead = mNextRecord;
            mNextRecord = advance();
            return currentRead;
        }

        /**
         * Closes down the existing iterator.
         */
        public void close() {
            if (!isClosed) {
            if (this != mCurrentIterator) {
                throw new IllegalStateException("Attempt to close non-current iterator");
            }
                mCurrentIterator = null;
                isClosed = true;
            }
        }

        /**
         * @throws UnsupportedOperationException always.
         */
        public void remove() {
            throw new UnsupportedOperationException("Not supported: remove");
        }

        SAMRecord advance() {
            while (true) {
                // Pull next record from stream
                if(!wrappedIterator.hasNext())
                    return null;

                final SAMRecord record = wrappedIterator.next();
                // If beyond the end of this reference sequence, end iteration
                final int referenceIndex = record.getReferenceIndex();
                if (referenceIndex != mReferenceIndex) {
                    if (referenceIndex < 0 ||
                        referenceIndex > mReferenceIndex) {
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
                // If read is unmapped but has a coordinate, return it if the coordinate is within
                // the query region, regardless of whether the mapped mate will be returned.
                final int alignmentEnd;
                if (mQueryType == QueryType.STARTING_AT) {
                    alignmentEnd = -1;
                } else {
                    alignmentEnd = (record.getAlignmentEnd() != SAMRecord.NO_ALIGNMENT_START?
                            record.getAlignmentEnd(): alignmentStart);
                }

                if (alignmentStart > mRegionEnd) {
                    // If scanned beyond target region, end iteration
                    return null;
                }
                // Filter for overlap with region
                if (mQueryType == QueryType.CONTAINED) {
                    if (alignmentStart >= mRegionStart && alignmentEnd <= mRegionEnd) {
                        return record;
                    }
                } else if (mQueryType == QueryType.OVERLAPPING) {
                    if (alignmentEnd >= mRegionStart && alignmentStart <= mRegionEnd) {
                        return record;
                    }
                } else {
                    if (alignmentStart == mRegionStart) {
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
