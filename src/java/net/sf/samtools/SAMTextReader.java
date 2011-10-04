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


import net.sf.samtools.util.BufferedLineReader;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.io.InputStream;
import java.util.Map;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Internal class for reading SAM text files.
 */
class SAMTextReader extends SAMFileReader.ReaderImplementation {
    // From SAM specification
    private static final int QNAME_COL = 0;
    private static final int FLAG_COL = 1;
    private static final int RNAME_COL = 2;
    private static final int POS_COL = 3;
    private static final int MAPQ_COL = 4;
    private static final int CIGAR_COL = 5;
    private static final int MRNM_COL = 6;
    private static final int MPOS_COL = 7;
    private static final int ISIZE_COL = 8;
    private static final int SEQ_COL = 9;
    private static final int QUAL_COL = 10;

    private static final int NUM_REQUIRED_FIELDS = 11;

    // Read string must contain only these characters
    private static final Pattern VALID_BASES = Pattern.compile("^[acmgrsvtwyhkdbnACMGRSVTWYHKDBN.=]+$");

    private BufferedLineReader mReader;
    private SAMFileHeader mFileHeader = null;
    private String mCurrentLine = null;
    private RecordIterator mIterator = null;
    private File mFile = null;
    private final TextTagCodec tagCodec = new TextTagCodec();
    private SAMFileReader.ValidationStringency validationStringency = SAMFileReader.ValidationStringency.DEFAULT_STRINGENCY;

    /**
     * Add information about the origin (reader and position) to SAM records.
     */
    private SAMFileReader mParentReader;

    /**
     * Prepare to read a SAM text file.
     * @param stream Need not be buffered, as this class provides buffered reading.
     */
    SAMTextReader(final InputStream stream, final SAMFileReader.ValidationStringency validationStringency) {
        mReader = new BufferedLineReader(stream);
        this.validationStringency = validationStringency;
        readHeader();
    }

    /**
     * Prepare to read a SAM text file.
     * @param stream Need not be buffered, as this class provides buffered reading.
     * @param file For error reporting only.
     */
    SAMTextReader(final InputStream stream, final File file, final SAMFileReader.ValidationStringency validationStringency) {
        this(stream, validationStringency);
        mFile = file;
    }

    /**
     * If true, writes the source of every read into the source SAMRecords.
     * @param enabled true to write source information into each SAMRecord.
     */
    void enableFileSource(final SAMFileReader reader, final boolean enabled) {
        this.mParentReader = enabled ? reader : null;
    }

    void enableIndexCaching(final boolean enabled) {
        throw new UnsupportedOperationException("Cannot enable index caching for a SAM text reader");
    }

    void enableIndexMemoryMapping(final boolean enabled) {
        throw new UnsupportedOperationException("Cannot enable index memory mapping for a SAM text reader");
    }

    void enableCrcChecking(final boolean enabled) {
        // Do nothing - this has no meaning for SAM reading
    }

    boolean hasIndex() {
        return false;    
    }

    BAMIndex getIndex() {
        throw new UnsupportedOperationException();
    }

    void close() {
        if (mReader != null) {
            try {
                mReader.close();
            } finally {
                mReader = null;
            }
        }
    }

    SAMFileHeader getFileHeader() {
        return mFileHeader;
    }

    public SAMFileReader.ValidationStringency getValidationStringency() {
        return validationStringency;
    }

    public void setValidationStringency(final SAMFileReader.ValidationStringency stringency) {
        this.validationStringency = stringency;
    }

    /**
     * There can only be one extant iterator on a SAMTextReader at a time.  The previous one must
     * be closed before calling getIterator().  Because the input stream is not seekable, a subsequent
     * call to getIterator() returns an iterator that starts where the last one left off.
     *
     * @return Iterator of SAMRecords in file order.
     */
    CloseableIterator<SAMRecord> getIterator() {
        if (mReader == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (mIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        mIterator = new RecordIterator();
        return mIterator;
    }

    /**
     * Generally loads data at a given point in the file.  Unsupported for SAMTextReaders.
     * @param fileSpan The file span.
     * @return An iterator over the given file span.
     */
    CloseableIterator<SAMRecord> getIterator(final SAMFileSpan fileSpan) {
        throw new UnsupportedOperationException("Cannot directly iterate over regions within SAM text files.");
    }

    /**
     * Generally gets a pointer to the first read in the file.  Unsupported for SAMTextReaders.
     * @return An pointer to the first read in the file.
     */
    SAMFileSpan getFilePointerSpanningReads() {
        throw new UnsupportedOperationException("Cannot retrieve file pointers within SAM text files.");
    }

    /**
     * Unsupported for SAM text files.
     */
    CloseableIterator<SAMRecord> query(final String sequence, final int start, final int end, final boolean contained) {
        throw new UnsupportedOperationException("Cannot query SAM text files");
    }

    /**
     * Unsupported for SAM text files.
     */
    CloseableIterator<SAMRecord> queryAlignmentStart(final String sequence, final int start) {
        throw new UnsupportedOperationException("Cannot query SAM text files");
    }

    public CloseableIterator<SAMRecord> queryUnmapped() {
        throw new UnsupportedOperationException("Cannot query SAM text files");
    }

    private void readHeader() {
        final SAMTextHeaderCodec headerCodec = new SAMTextHeaderCodec();
        headerCodec.setValidationStringency(validationStringency);
        mFileHeader = headerCodec.decode(mReader, (mFile != null? mFile.toString(): null));
        advanceLine();
    }

    private String advanceLine() {
        mCurrentLine = mReader.readLine();
        return mCurrentLine;
    }

    private String makeErrorString(final String reason) {
        String fileMessage = "";
        if (mFile != null) {
            fileMessage = "File " + mFile + "; ";
        }
        return "Error parsing text SAM file. " + reason + "; " + fileMessage +
                "Line " + mReader.getLineNumber() + "\nLine: " + mCurrentLine;
    }

    private RuntimeException reportFatalErrorParsingLine(final String reason) {
        return new SAMFormatException(makeErrorString(reason));
    }

    private void reportErrorParsingLine(final String reason) {
        final String errorMessage = makeErrorString(reason);

        if (validationStringency == SAMFileReader.ValidationStringency.STRICT) {
            throw new SAMFormatException(errorMessage);
        } else if (validationStringency == SAMFileReader.ValidationStringency.LENIENT) {
            System.err.println("Ignoring SAM validation error due to lenient parsing:");
            System.err.println(errorMessage);
        }
    }

    private void reportErrorParsingLine(final Exception e) {
        final String errorMessage = makeErrorString(e.getMessage());
        if (validationStringency == SAMFileReader.ValidationStringency.STRICT) {
            throw new SAMFormatException(errorMessage);
        } else if (validationStringency == SAMFileReader.ValidationStringency.LENIENT) {
            System.err.println("Ignoring SAM validation error due to lenient parsing:");
            System.err.println(errorMessage);
        }
    }

    /**
     * SAMRecord iterator for SAMTextReader
     */
    private class RecordIterator implements CloseableIterator<SAMRecord> {

        /**
         * Allocate this once rather than for every line as a performance optimization.
         * The size is arbitrary -- merely large enough to handle the maximum number
         * of fields we might expect from a reasonable SAM file.
         */
        private final String[] mFields = new String[10000];

        private RecordIterator() {
            if (mReader == null) {
                throw new IllegalStateException("Reader is closed.");
            }
        }

        public void close() {
            SAMTextReader.this.close();
        }

        public boolean hasNext() {
            return mCurrentLine != null;
        }

        public SAMRecord next() {
            if (!hasNext()) {
                throw new IllegalStateException("Cannot call next() on exhausted iterator");
            }
            try {
                return parseLine();
            } finally {
                advanceLine();
            }
        }

        public void remove() {
            throw new UnsupportedOperationException("Not supported: remove");
        }

        int parseInt(final String s, final String fieldName) {
            final int ret;
            try {
                ret = Integer.parseInt(s);
            } catch (NumberFormatException e) {
                throw reportFatalErrorParsingLine("Non-numeric value in " + fieldName + " column");
            }
            return ret;
        }

        void validateReferenceName(final String rname, final String fieldName) {
            if (rname.equals("=")) {
                if (fieldName.equals("MRNM")) {
                    return;
                }
                reportErrorParsingLine("= is not a valid value for " + fieldName + " field.");
            }
            if (getFileHeader().getSequenceDictionary().size() != 0) {
                if (getFileHeader().getSequence(rname) == null) {
                    reportErrorParsingLine(fieldName + " '" + rname + "' not found in any SQ record");
                }
            }
        }

        private SAMRecord parseLine() {
            final int numFields = StringUtil.split(mCurrentLine, mFields, '\t');
            if (numFields < NUM_REQUIRED_FIELDS) {
                throw reportFatalErrorParsingLine("Not enough fields");
            }
            if (numFields == mFields.length) {
                reportErrorParsingLine("Too many fields in SAM text record.");
            }
            for (int i = 0; i < numFields; ++i) {
                if (mFields[i].length() == 0) {
                    reportErrorParsingLine("Empty field at position " + i + " (zero-based)");
                }
            }
            final SAMRecord samRecord = new SAMRecord(mFileHeader);
            samRecord.setValidationStringency(getValidationStringency());
            if(mParentReader != null)
                samRecord.setFileSource(new SAMFileSource(mParentReader,null));
            samRecord.setHeader(mFileHeader);
            samRecord.setReadName(mFields[QNAME_COL]);

            final int flags = parseInt(mFields[FLAG_COL], "FLAG");
            samRecord.setFlags(flags);

            String rname = mFields[RNAME_COL];
            if (!rname.equals("*")) {
                rname = SAMSequenceRecord.truncateSequenceName(rname);
                validateReferenceName(rname, "RNAME");
                samRecord.setReferenceName(rname);
            } else if (!samRecord.getReadUnmappedFlag()) {
                    reportErrorParsingLine("RNAME is not specified but flags indicate mapped");
                }

            final int pos = parseInt(mFields[POS_COL], "POS");
            final int mapq = parseInt(mFields[MAPQ_COL], "MAPQ");
            final String cigar = mFields[CIGAR_COL];
            if (!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(samRecord.getReferenceName())) {
                if (pos == 0) {
                    reportErrorParsingLine("POS must be non-zero if RNAME is specified");
                }
                if (!samRecord.getReadUnmappedFlag() && cigar.equals("*")) {
                    reportErrorParsingLine("CIGAR must not be '*' if RNAME is specified");
                }
            } else {
                if (pos != 0) {
                    reportErrorParsingLine("POS must be zero if RNAME is not specified");
                }
                if (mapq != 0) {
                    reportErrorParsingLine("MAPQ must be zero if RNAME is not specified");
                }
                if (!cigar.equals("*")) {
                    reportErrorParsingLine("CIGAR must be '*' if RNAME is not specified");
                }
            }
            samRecord.setAlignmentStart(pos);
            samRecord.setMappingQuality(mapq);
            samRecord.setCigarString(cigar);

            String mateRName = mFields[MRNM_COL];
            if (mateRName.equals("*")) {
                if (samRecord.getReadPairedFlag() && !samRecord.getMateUnmappedFlag()) {
                    reportErrorParsingLine("MRNM not specified but flags indicate mate mapped");
                }
            }
            else {
                if (!samRecord.getReadPairedFlag()) {
                    reportErrorParsingLine("MRNM specified but flags indicate unpaired");
                }
                if (!"=".equals(mateRName)) {
                    mateRName = SAMSequenceRecord.truncateSequenceName(mateRName);
                }
                validateReferenceName(mateRName, "MRNM");
                if (mateRName.equals("=")) {
                    if (samRecord.getReferenceName() == null) {
                        reportErrorParsingLine("MRNM is '=', but RNAME is not set");
                    }
                    samRecord.setMateReferenceName(samRecord.getReferenceName());
                } else {
                    samRecord.setMateReferenceName(mateRName);
                }
            }

            final int matePos = parseInt(mFields[MPOS_COL], "MPOS");
            final int isize = parseInt(mFields[ISIZE_COL], "ISIZE");
            if (!samRecord.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
                if (matePos == 0) {
                    reportErrorParsingLine("MPOS must be non-zero if MRNM is specified");
                }
            } else {
                if (matePos != 0) {
                    reportErrorParsingLine("MPOS must be zero if MRNM is not specified");
                }
                if (isize != 0) {
                    reportErrorParsingLine("ISIZE must be zero if MRNM is not specified");
                }
            }
            samRecord.setMateAlignmentStart(matePos);
            samRecord.setInferredInsertSize(isize);
            if (!mFields[SEQ_COL].equals("*")) {
                validateReadBases(mFields[SEQ_COL]);
                samRecord.setReadString(mFields[SEQ_COL]);
            } else {
                samRecord.setReadBases(SAMRecord.NULL_SEQUENCE);
            }
            if (!mFields[QUAL_COL].equals("*")) {
                if (samRecord.getReadBases() == SAMRecord.NULL_SEQUENCE) {
                    reportErrorParsingLine("QUAL should not be specified if SEQ is not specified");
                }
                if (samRecord.getReadString().length() != mFields[QUAL_COL].length()) {
                    reportErrorParsingLine("length(QUAL) != length(SEQ)");
                }
                samRecord.setBaseQualityString(mFields[QUAL_COL]);
            } else {
                samRecord.setBaseQualities(SAMRecord.NULL_QUALS);
            }

            for (int i = NUM_REQUIRED_FIELDS; i < numFields; ++i) {
                parseTag(samRecord, mFields[i]);
            }

            final List<SAMValidationError> validationErrors = samRecord.isValid();
            if (validationErrors != null) {
                for (final SAMValidationError errorMessage : validationErrors) {
                    reportErrorParsingLine(errorMessage.getMessage());
                }
            }
            return samRecord;
        }

        private void validateReadBases(final String bases) {
            if (!VALID_BASES.matcher(bases).matches()) {
                reportErrorParsingLine("Invalid character in read bases");
            }
        }

        private void parseTag(final SAMRecord samRecord, final String tag) {
            Map.Entry<String, Object> entry = null;
            try {
                entry = tagCodec.decode(tag);
            } catch (SAMFormatException e) {
                reportErrorParsingLine(e);
            }
            if (entry != null) {
                if (entry.getValue() instanceof TagValueAndUnsignedArrayFlag) {
                    final TagValueAndUnsignedArrayFlag valueAndFlag = (TagValueAndUnsignedArrayFlag) entry.getValue();
                    if (valueAndFlag.isUnsignedArray) {
                        samRecord.setUnsignedArrayAttribute(entry.getKey(), valueAndFlag.value);
                    }
                    else {
                        samRecord.setAttribute(entry.getKey(), valueAndFlag.value);
                    }
                } else {
                    samRecord.setAttribute(entry.getKey(), entry.getValue());
                }
            }
        }
    }
}

