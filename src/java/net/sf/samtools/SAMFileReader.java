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

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.net.URL;


/**
 * Class for reading and querying SAM/BAM files.  Delegates to appropriate concrete implementation.
 */
public class SAMFileReader implements Iterable<SAMRecord>, Closeable {

    private static ValidationStringency defaultValidationStringency = ValidationStringency.DEFAULT_STRINGENCY;

    public static ValidationStringency getDefaultValidationStringency() {
        return defaultValidationStringency;
    }

    /**
     * Set validation stringency for all subsequently-created SAMFileReaders.  This is the only way to
     * change the validation stringency for SAM header.
     */
    public static void setDefaultValidationStringency(final ValidationStringency defaultValidationStringency) {
        SAMFileReader.defaultValidationStringency = defaultValidationStringency;
    }

    private boolean mIsBinary = false;
    private BAMFileIndex mFileIndex = null;
    private ReaderImplementation mReader = null;
    private File samFile = null;

    /**
     * How strict to be when reading a SAM or BAM, beyond bare minimum validation.
     */
    public enum ValidationStringency {
        /**
         * Do the right thing, throw an exception if something looks wrong.
         */
        STRICT,
        /**
         * Emit warnings but keep going if possible.
         */
        LENIENT,
        /**
         * Like LENIENT, only don't emit warning messages.
         */
        SILENT;

        public static final ValidationStringency DEFAULT_STRINGENCY = STRICT;
    }

    /**
     * Internal interface for SAM/BAM file reader implementations.
     * Implemented as an abstract class to enforce better access control.
     */
    static abstract class ReaderImplementation {
        abstract SAMFileHeader getFileHeader();
        abstract CloseableIterator<SAMRecord> getIterator();
        abstract CloseableIterator<SAMRecord> query(String sequence, int start, int end, boolean contained);
        abstract CloseableIterator<SAMRecord> queryAlignmentStart(String sequence, int start);
        abstract public CloseableIterator<SAMRecord> queryUnmapped();
        abstract void close();
        // If true, emit warnings about format errors rather than throwing exceptions;
        abstract void setValidationStringency(final ValidationStringency validationStringency);
        abstract ValidationStringency getValidationStringency();
    }


    /**
     * Prepare to read a SAM or BAM file.  Indexed lookup not allowed because reading from InputStream.
     */
    public SAMFileReader(final InputStream stream) {
        this(stream, false);
    }

    /**
     * Prepare to read a SAM or BAM file.  If the given file is a BAM, and has a companion BAI index file
     * that is named according to the convention, it will be found and opened, and indexed query will be allowed.
     */
    public SAMFileReader(final File file) {
        this(file, null, false);
    }

    /**
     * Prepare to read a SAM or BAM file.  If the given file is a BAM, and an index is present, indexed query
     * will be allowed.
     *
     * @param file SAM or BAM to read.
     * @param indexFile Index file that is companion to BAM, or null if no index file, or if index file
     * should be found automatically.
     */
    public SAMFileReader(final File file, final File indexFile) {
        this(file, indexFile, false);
    }

    /**
     * Read a SAM or BAM file.  Indexed lookup not allowed because reading from InputStream.
     *
     * @param stream input SAM or BAM.
     * @param eagerDecode if true, decode SAM record entirely when reading it.
     */
    public SAMFileReader(final InputStream stream, final boolean eagerDecode) {
        init(stream, eagerDecode, defaultValidationStringency);
    }

    /**
     * Read a SAM or BAM file, possibly with an index file if present.
     * If the given file is a BAM, and an index is present, indexed query will be allowed.
     *
     * @param file SAM or BAM.
     * @param eagerDecode if true, decode SAM record entirely when reading it.
     */
    public SAMFileReader(final File file, final boolean eagerDecode) {
        init(file, null, eagerDecode, defaultValidationStringency);
    }

    /**
     * Read a SAM or BAM file, possibly with an index file. If the given file is a BAM, and an index is present,
     * indexed query will be allowed.
     *
     * @param file SAM or BAM.
     * @param indexFile Location of index file, or null in order to use the default index file (if present).
     * @param eagerDecode eagerDecode if true, decode SAM record entirely when reading it.
     */
    public SAMFileReader(final File file, final File indexFile, final boolean eagerDecode){
        init(file, indexFile, eagerDecode, defaultValidationStringency);
    }

    /**
     * Read a BAM file by http
     * indexed query will be allowed.
     *
     * @param url         BAM.
     * @param indexFile   Location of index file, or null if indexed access not required.
     * @param eagerDecode eagerDecode if true, decode SAM record entirely when reading it.
     */
    public SAMFileReader(final URL url, final File indexFile, final boolean eagerDecode) {
        init(new SeekableBufferedStream(new SeekableHTTPStream(url)),
                indexFile, eagerDecode, defaultValidationStringency);
    }

    /**
     * Read a BAM file via caller-supplied mechanism.  Indexed query will be allowed, but
     * index file must be provided in that case.
     * @param strm BAM -- If the stream is not buffered, caller should wrap in SeekableBufferedStream for
     * better performance.
     * @param indexFile Location of index file, or null indexed access not required.
     * @param eagerDecode if true, decode SAM record entirely when reading it.
     */
    public SAMFileReader(final SeekableStream strm, final File indexFile, final boolean eagerDecode) {
        init(strm, indexFile, eagerDecode, defaultValidationStringency);
    }

    public void close() {
        if (mReader != null) {
            mReader.close();
        }
        if (mFileIndex != null) {
            mFileIndex.close();
        }
        mReader = null;
        mFileIndex = null;
    }

    /**
     * @return True if this is a BAM reader.
     */
    public boolean isBinary() {
        return mIsBinary;
    }

    /**
     * @return true if ths is a BAM file, and has an index
     */
    public boolean hasIndex() {
        return (mFileIndex != null);
    }

    public SAMFileHeader getFileHeader() {
        return mReader.getFileHeader();
    }

    /**
     * Control validation of SAMRecords as they are read from file.
     * In order to control validation stringency for SAM Header, call SAMFileReader.setDefaultValidationStringency
     * before constructing a SAMFileReader.
     */
    public void setValidationStringency(final ValidationStringency validationStringency) {
        mReader.setValidationStringency(validationStringency);
    }

    /**
     * Iterate through file in order.  For a SAMFileReader constructed from an InputStream, and for any SAM file,
     * a 2nd iteration starts where the 1st one left off.  For a BAM constructed from a File, each new iteration
     * starts at the first record.
     * <p/>
     * Only a single open iterator on a SAM or BAM file may be extant at any one time.  If you want to start
     * a second iteration, the first one must be closed first.
     */
    public CloseableIterator<SAMRecord> iterator() {
        return mReader.getIterator();
    }

    /**
     * Iterate over records that match the given interval.  Only valid to call this if hasIndex() == true.
     * <p/>
     * Only a single open iterator on a given SAMFileReader may be extant at any one time.  If you want to start
     * a second iteration, the first one must be closed first.  You can use a second SAMFileReader to iterate
     * in parallel over the same underlying file.
     * <p/>
     * Note that indexed lookup is not perfectly efficient in terms of disk I/O.  I.e. some SAMRecords may be read
     * and then discarded because they do not match the interval of interest.
     * <p/>
     * Note that an unmapped read will be returned by this call if it has a coordinate for the purpose of sorting that
     * is in the query region.
     *
     * @param sequence  Reference sequence of interest.
     * @param start     1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
     * @param end       1-based, inclusive end of interval of interest. Zero implies end of the reference sequence.
     * @param contained If true, each SAMRecord returned is will have its alignment completely contained in the
     *                  interval of interest.  If false, the alignment of the returned SAMRecords need only overlap the interval of interest.
     * @return Iterator over the SAMRecords matching the interval.
     */
    public CloseableIterator<SAMRecord> query(final String sequence, final int start, final int end, final boolean contained) {
        return mReader.query(sequence, start, end, contained);
    }

    /**
     * Iterate over records that overlap the given interval.  Only valid to call this if hasIndex() == true.
     * <p/>
     * Only a single open iterator on a given SAMFileReader may be extant at any one time.  If you want to start
     * a second iteration, the first one must be closed first.
     * <p/>
     * Note that indexed lookup is not perfectly efficient in terms of disk I/O.  I.e. some SAMRecords may be read
     * and then discarded because they do not match the interval of interest.
     * <p/>
     * Note that an unmapped read will be returned by this call if it has a coordinate for the purpose of sorting that
     * is in the query region.
     *
     * @param sequence Reference sequence of interest.
     * @param start    1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
     * @param end      1-based, inclusive end of interval of interest. Zero implies end of the reference sequence.
     * @return Iterator over the SAMRecords overlapping the interval.
     */
    public CloseableIterator<SAMRecord> queryOverlapping(final String sequence, final int start, final int end) {
        return query(sequence, start, end, false);
    }

    /**
     * Iterate over records that are contained in the given interval.  Only valid to call this if hasIndex() == true.
     * <p/>
     * Only a single open iterator on a given SAMFileReader may be extant at any one time.  If you want to start
     * a second iteration, the first one must be closed first.
     * <p/>
     * Note that indexed lookup is not perfectly efficient in terms of disk I/O.  I.e. some SAMRecords may be read
     * and then discarded because they do not match the interval of interest.
     * <p/>
     * Note that an unmapped read will be returned by this call if it has a coordinate for the purpose of sorting that
     * is in the query region.
     *
     * @param sequence Reference sequence of interest.
     * @param start    1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
     * @param end      1-based, inclusive end of interval of interest. Zero implies end of the reference sequence.
     * @return Iterator over the SAMRecords contained in the interval.
     */
    public CloseableIterator<SAMRecord> queryContained(final String sequence, final int start, final int end) {
        return query(sequence, start, end, true);
    }

    public CloseableIterator<SAMRecord> queryUnmapped() {
        return mReader.queryUnmapped();
    }

    /**
     * Iterate over records that map to the given sequence and start at the given position.  Only valid to call this if hasIndex() == true.
     * <p/>
     * Only a single open iterator on a given SAMFileReader may be extant at any one time.  If you want to start
     * a second iteration, the first one must be closed first.
     * <p/>
     * Note that indexed lookup is not perfectly efficient in terms of disk I/O.  I.e. some SAMRecords may be read
     * and then discarded because they do not match the interval of interest.
     * <p/>
     * Note that an unmapped read will be returned by this call if it has a coordinate for the purpose of sorting that
     * matches the arguments.
     *
     * @param sequence Reference sequence of interest.
     * @param start    Alignment start of interest.
     * @return Iterator over the SAMRecords with the given alignment start.
     */
    public CloseableIterator<SAMRecord> queryAlignmentStart(final String sequence, final int start) {
        return mReader.queryAlignmentStart(sequence, start);
    }

    /**
     * Fetch the mate for the given read.  Only valid to call this if hasIndex() == true.
     * This will work whether the mate has a coordinate or not, so long as the given read has correct
     * mate information.  This method iterates over the SAM file, so there may not be an unclosed
     * iterator on the SAM file when this method is called.
     *
     * @param rec Record for which mate is sought.  Must be a paired read.
     * @return rec's mate, or null if it cannot be found.
     */
    public SAMRecord queryMate(final SAMRecord rec) {
        if (!rec.getReadPairedFlag()) {
            throw new IllegalArgumentException("queryMate called for unpaired read.");
        }
        if (rec.getFirstOfPairFlag() == rec.getSecondOfPairFlag()) {
            throw new IllegalArgumentException("SAMRecord must be either first and second of pair, but not both.");
        }
        final boolean firstOfPair = rec.getFirstOfPairFlag();
        final CloseableIterator<SAMRecord> it;
        if (rec.getMateReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            it = queryUnmapped();
        } else {
            it = queryAlignmentStart(rec.getMateReferenceName(), rec.getMateAlignmentStart());
        }
        try {
            SAMRecord mateRec = null;
            while (it.hasNext()) {
                final SAMRecord next = it.next();
                if (!next.getReadPairedFlag()) {
                    if (rec.getReadName().equals(next.getReadName())) {
                        throw new SAMFormatException("Paired and unpaired reads with same name: " + rec.getReadName());
                    }
                    continue;
                }
                if (firstOfPair) {
                    if (next.getFirstOfPairFlag()) continue;
                } else {
                    if (next.getSecondOfPairFlag()) continue;
                }
                if (rec.getReadName().equals(next.getReadName())) {
                    if (mateRec != null) {
                        throw new SAMFormatException("Multiple SAMRecord with read name " + rec.getReadName() +
                                " for " + (firstOfPair ? "second" : "first") + " end.");
                    }
                    mateRec = next;
                }
            }
            return mateRec;
        } finally {
            it.close();
        }
    }

    private void init(final InputStream stream, final boolean eagerDecode, final ValidationStringency validationStringency) {

        try {
            final BufferedInputStream bufferedStream = IOUtil.toBufferedStream(stream);
            if (isBAMFile(bufferedStream)) {
                mIsBinary = true;
                mReader = new BAMFileReader(bufferedStream, eagerDecode, validationStringency);
            } else if (isGzippedSAMFile(bufferedStream)) {
                mIsBinary = false;
                mReader = new SAMTextReader(new GZIPInputStream(bufferedStream), validationStringency);
            } else if (isSAMFile(bufferedStream)) {
                mIsBinary = false;
                mReader = new SAMTextReader(bufferedStream, validationStringency);
            } else {
                throw new SAMFormatException("Unrecognized file format");
            }
            setValidationStringency(validationStringency);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }


    private void init(final SeekableStream strm, final File indexFile, final boolean eagerDecode,
                      final ValidationStringency validationStringency) {

        try {
            // Its too expensive to examine the remote file to determine type.
            // Rely on file extension.
            if (strm.getSource() == null || strm.getSource().toLowerCase().endsWith(".bam")) {
                mIsBinary = true;
                final BAMFileReader reader = new BAMFileReader(strm, eagerDecode, validationStringency);
                mReader = reader;
                if (indexFile != null) {
                    mFileIndex = new BAMFileIndex(indexFile);
                    reader.setFileIndex(mFileIndex);
                }
            } else {
                throw new SAMFormatException("Unrecognized file format: " + strm);
            }
            setValidationStringency(validationStringency);
        }
        catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }


    private void init(final File file, File indexFile, final boolean eagerDecode, final ValidationStringency validationStringency) {
        this.samFile = file;

        try {
            final BufferedInputStream bufferedStream = new BufferedInputStream(new FileInputStream(file));
            if (isBAMFile(bufferedStream)) {
                mIsBinary = true;
                if (!file.isFile()) {
                    // Handle case in which file is a named pipe, e.g. /dev/stdin or created by mkfifo
                    mReader = new BAMFileReader(bufferedStream, eagerDecode, validationStringency);
                } else {
                    bufferedStream.close();
                    final BAMFileReader reader = new BAMFileReader(file, eagerDecode, validationStringency);
                    mReader = reader;
                    if (indexFile == null) {
                        indexFile = findIndexFile(file);
                    }
                    if (indexFile != null) {
                        mFileIndex = new BAMFileIndex(indexFile);
                        reader.setFileIndex(mFileIndex);
                        if (indexFile.lastModified() < file.lastModified()) {
                            System.err.println("WARNING: BAM index file " + indexFile.getAbsolutePath() +
                                    " is older than BAM " + file.getAbsolutePath());
                        }
                    }
                }
            } else if (isGzippedSAMFile(bufferedStream)) {
                mIsBinary = false;
                mReader = new SAMTextReader(new GZIPInputStream(bufferedStream), validationStringency);
            } else if (isSAMFile(bufferedStream)) {
                if (indexFile != null) {
                    bufferedStream.close();
                    throw new RuntimeException("Cannot use index file with textual SAM file");
                }
                mIsBinary = false;
                mReader = new SAMTextReader(bufferedStream, file, validationStringency);
            } else {
                bufferedStream.close();
                throw new SAMFormatException("Unrecognized file format");
            }
            setValidationStringency(validationStringency);
        }
        catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }


    /**
     * Look for BAM index file according to standard naming convention.
     *
     * @param dataFile BAM file name.
     * @return Index file name, or null if not found.
     */
    private File findIndexFile(final File dataFile) {
        // If input is foo.bam, look for foo.bai
        final String bamExtension = ".bam";
        File indexFile;
        final String fileName = dataFile.getName();
        if (fileName.endsWith(bamExtension)) {
            final String bai = fileName.substring(0, fileName.length() - bamExtension.length()) + ".bai";
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

    /**
     * @param stream stream.markSupported() must be true
     * @return true if this looks like a BAM file.
     */
    private boolean isBAMFile(final InputStream stream)
            throws IOException {
        return BlockCompressedInputStream.isValidFile(stream);
    }

    /**
     * Attempts to check whether the file is a gzipped sam file.  Returns true if it
     * is and false otherwise.
     */
    private boolean isGzippedSAMFile(final BufferedInputStream stream) {
        if (!stream.markSupported()) {
            throw new IllegalArgumentException("Cannot test a stream that doesn't support marking.");
        }
        stream.mark(8000);

        try {
            final GZIPInputStream gunzip = new GZIPInputStream(stream);
            final int ch = gunzip.read();
            return true;
        }
        catch (IOException ioe) {
            return false;
        }
        finally {
            try {
                stream.reset();
            }
            catch (IOException ioe) {
                throw new IllegalStateException("Could not reset stream.");
            }
        }
    }

    private boolean isSAMFile(final InputStream stream) {
        // For now, assume every non-binary file is a SAM text file.
        return true;
    }

    @Override
    public String toString() {
        if (this.samFile == null) {
            return getClass().getSimpleName() + "{initialized with stream}";
        } else {
            return getClass().getSimpleName() + "{" + this.samFile.getAbsolutePath() + "}";
        }
    }
}
