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


import net.sf.samtools.util.AsciiLineReader;
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
class SAMTextReader
    extends SAMFileReader.ReaderImplementation
{
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
    private static final Pattern VALID_BASES = Pattern.compile("^[acgtnACGTN.=]+$");

    private AsciiLineReader mReader;
    private SAMFileHeader mFileHeader = null;
    private String mCurrentLine = null;
    private RecordIterator mIterator = null;
    private File mFile = null;
    private final TextTagCodec tagCodec = new TextTagCodec();
    private SAMFileReader.ValidationStringency validationStringency = SAMFileReader.ValidationStringency.DEFAULT_STRINGENCY;

    /**
     * Prepare to read a SAM text file.
     * @param stream Need not be buffered, as this class provides buffered reading.
     */
    SAMTextReader(final InputStream stream) {
        mReader = new AsciiLineReader(stream);
        readHeader();
    }

    /**
     * Prepare to read a SAM text file.
     * @param stream Need not be buffered, as this class provides buffered reading.
     * @param file For error reporting only.
     */
    SAMTextReader(final InputStream stream, final File file) {
        this(stream);
        mFile = file;
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
     * Unsupported for SAM text files.
     */
    CloseableIterator<SAMRecord> query(final String sequence, final int start, final int end, final boolean contained) {
        throw new UnsupportedOperationException("Cannot query SAM text files");
    }

    private void readHeader() {
        final SAMTextHeaderCodec headerCodec = new SAMTextHeaderCodec();
        mFileHeader = headerCodec.decode(mReader, mFile);
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

        private SAMRecord mCurrentRecord;

        private RecordIterator() {
            assert(mReader != null);
            if (mCurrentLine != null) {
                parseLine();
            }

        }

        public void close() {
            mCurrentRecord = null;
            SAMTextReader.this.close();
        }

        public boolean hasNext() {
            return mCurrentRecord != null;
        }

        public SAMRecord next() {
            if (!hasNext()) {
                throw new IllegalStateException("Cannot call next() on exhausted iterator");
            }
            final SAMRecord ret = mCurrentRecord;
            mCurrentRecord = null;
            advanceLine();
            if (mCurrentLine != null) {
                parseLine();
            }
            return ret;
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
            if (fieldName.equals("MRNM") && rname.equals("=")) {
                return;
            }
            if (getFileHeader().getSequenceDictionary().size() != 0) {
                if (getFileHeader().getSequence(rname) == null) {
                    reportErrorParsingLine(fieldName + " '" + rname + "' not found in any SQ record");
                }
            }
        }

        private void parseLine() {
            final int numFields = StringUtil.split(mCurrentLine, mFields, '\t');
            if (numFields < NUM_REQUIRED_FIELDS) {
                reportErrorParsingLine("Not enough fields");
            }
            if (numFields == mFields.length) {
                reportErrorParsingLine("Too many fields in SAM text record.");
            }
            for (int i = 0; i < numFields; ++i) {
                if (mFields[i].length() == 0) {
                    reportErrorParsingLine("Empty field at position " + i + " (zero-based)");
                }
            }
            mCurrentRecord = new SAMRecord(mFileHeader);
            mCurrentRecord.setValidationStringency(getValidationStringency());
            mCurrentRecord.setHeader(mFileHeader);
            mCurrentRecord.setReadName(mFields[QNAME_COL]);

            final int flags = parseInt(mFields[FLAG_COL], "FLAG");
            mCurrentRecord.setFlags(flags);

            final String rname = mFields[RNAME_COL];
            if (!rname.equals("*")) {
                validateReferenceName(rname, "RNAME");
                mCurrentRecord.setReferenceName(rname);
            } else if (!mCurrentRecord.getReadUnmappedFlag()) {
                    reportErrorParsingLine("RNAME is not specified but flags indicate mapped");
                }

            final int pos = parseInt(mFields[POS_COL], "POS");
            final int mapq = parseInt(mFields[MAPQ_COL], "MAPQ");
            final String cigar = mFields[CIGAR_COL];
            if (!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(mCurrentRecord.getReferenceName())) {
                if (pos == 0) {
                    reportErrorParsingLine("POS must be non-zero if RNAME is specified");
                }
                if (!mCurrentRecord.getReadUnmappedFlag() && cigar.equals("*")) {
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
            mCurrentRecord.setAlignmentStart(pos);
            mCurrentRecord.setMappingQuality(mapq);
            mCurrentRecord.setCigarString(cigar);

            final String mateRName = mFields[MRNM_COL];
            if (mateRName.equals("*")) {
                if (mCurrentRecord.getReadPairedFlag() && !mCurrentRecord.getMateUnmappedFlag()) {
                    reportErrorParsingLine("MRNM not specified but flags indicate mate mapped");
                }
            }
            else {
                if (!mCurrentRecord.getReadPairedFlag()) {
                    reportErrorParsingLine("MRNM specified but flags indicate unpaired");
                }

                validateReferenceName(mateRName, "MRNM");
                if (mateRName.equals("=")) {
                    if (mCurrentRecord.getReferenceName() == null) {
                        reportErrorParsingLine("MRNM is '=', but RNAME is not set");
                    }
                    mCurrentRecord.setMateReferenceName(mCurrentRecord.getReferenceName());
                } else {
                    mCurrentRecord.setMateReferenceName(mateRName);
                }
            }

            final int matePos = parseInt(mFields[MPOS_COL], "MPOS");
            final int isize = parseInt(mFields[ISIZE_COL], "ISIZE");
            if (!mCurrentRecord.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
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
            mCurrentRecord.setMateAlignmentStart(matePos);
            mCurrentRecord.setInferredInsertSize(isize);
            if (!mFields[SEQ_COL].equals("*")) {
                validateReadBases(mFields[SEQ_COL]);
                mCurrentRecord.setReadString(mFields[SEQ_COL]);
            }
            if (!mFields[QUAL_COL].equals("*")) {
                if (mCurrentRecord.getReadString() == null) {
                    reportErrorParsingLine("QUAL should not be specified if SEQ is not specified");
                }
                if (mCurrentRecord.getReadString().length() != mFields[QUAL_COL].length()) {
                    reportErrorParsingLine("length(QUAL) != length(SEQ)");
                }
                mCurrentRecord.setBaseQualityString(mFields[QUAL_COL]);
            }

            for (int i = NUM_REQUIRED_FIELDS; i < numFields; ++i) {
                parseTag(mFields[i]);
            }

            final List<SAMValidationError> validationErrors = mCurrentRecord.isValid();
            if (validationErrors != null) {
                for (final SAMValidationError errorMessage : validationErrors) {
                    reportErrorParsingLine(errorMessage.getMessage());
                }
            }
        }

        private void validateReadBases(final String bases) {
            if (!VALID_BASES.matcher(bases).matches()) {
                reportErrorParsingLine("Invalid character in read bases");
            }
        }

        private void parseTag(final String tag) {
            Map.Entry<String, Object> entry = null;
            try {
                entry = tagCodec.decode(tag);
            } catch (SAMFormatException e) {
                reportErrorParsingLine(e);
            }
            if (entry != null) {
                mCurrentRecord.setAttribute(entry.getKey(), entry.getValue());
            }
        }
    }
}

