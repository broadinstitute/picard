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

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.*;

/**
 * Parser for a SAM text header, and a generator of SAM text header.
 */
public class SAMTextHeaderCodec {
    private static final String HEADER_LINE_START = "@";

    // These attributes are populated when parsing or generating
    private SAMFileHeader mFileHeader;
    private final TextTagCodec mTagCodec = new TextTagCodec();

    // These attributes are populated when parsing text
    private String mCurrentLine;
    private LineReader mReader;
    private File mFile;
    private List<SAMSequenceRecord> sequences;
    private List<SAMReadGroupRecord> readGroups;

    // For error reporting when parsing
    private SAMFileReader.ValidationStringency validationStringency = SAMFileReader.ValidationStringency.SILENT;

    // These attributes are populated when generating text
    private BufferedWriter writer;

    private static final String TAG_KEY_VALUE_SEPARATOR = ":";
    private static final String FIELD_SEPARATOR = "\t";

    /**
     * Reads text SAM header and converts to a SAMFileHeader object.
     * @param reader Where to get header text from.
     * @param file Name of the input file, for error messages.  May be null.
     * @return complete header object.
     */
    public SAMFileHeader decode(final LineReader reader, final File file) {
        mFileHeader = new SAMFileHeader();
        mReader = reader;
        mFile = file;
        sequences = new ArrayList<SAMSequenceRecord>();
        readGroups = new ArrayList<SAMReadGroupRecord>();

        while (advanceLine() != null) {
            final ParsedHeaderLine parsedHeaderLine = new ParsedHeaderLine(mCurrentLine);
            switch (parsedHeaderLine.getHeaderRecordType()) {

                case HD:
                    parseHDLine(parsedHeaderLine);
                    break;
                case PG:
                    parsePGLine(parsedHeaderLine);
                    break;
                case RG:
                    parseRGLine(parsedHeaderLine);
                    break;
                case SQ:
                    parseSQLine(parsedHeaderLine);
                    break;
                default:
                    throw new IllegalStateException("Unrecognized header record type: " +
                            parsedHeaderLine.getHeaderRecordType());
            }
        }
        mFileHeader.setSequenceDictionary(new SAMSequenceDictionary(sequences));
        mFileHeader.setReadGroups(readGroups);
        return mFileHeader;
    }

    private String advanceLine() {
        final int nextChar = mReader.peek();
        if (nextChar != '@') {
            return null;
        }
        mCurrentLine = mReader.readLine();
        return mCurrentLine;
    }

    /**
     * Transfer standard and non-standard tags from text representation to in-memory representation.
     * Standard tags are treated as Strings.  Non-standard tags are typed.
     * @param record attributes get set into this object.
     * @param textAttributes Map of tag type to value.  Some values may be removed by this method.
     */
    private void transferAttributes(final AbstractSAMHeaderRecord record, final Map<String, String> textAttributes) {
        // Transfer standard tags that are of type String
        for (final String standardTag : record.getStandardTags()) {
            final String value = textAttributes.remove(standardTag);
            if (value != null) {
                record.setAttribute(standardTag, value);
            }
        }
        // Transfer any non-standard typed tags
        for (final Map.Entry<String, String> entry : textAttributes.entrySet()) {
            final Object value = mTagCodec.decodeTypeAndValue(entry.getValue());
            record.setAttribute(entry.getKey(), value);
        }

    }

    private void parsePGLine(final ParsedHeaderLine parsedHeaderLine) {
        assert(HeaderRecordType.PG.equals(parsedHeaderLine.getHeaderRecordType()));
        parsedHeaderLine.requireTag(SAMProgramRecord.PROGRAM_GROUP_ID_TAG);
        final SAMProgramRecord programRecord = new SAMProgramRecord(parsedHeaderLine.removeValue(SAMProgramRecord.PROGRAM_GROUP_ID_TAG));

        transferAttributes(programRecord, parsedHeaderLine.mKeyValuePairs);
        mFileHeader.addProgramRecord(programRecord);
    }

    private void parseRGLine(final ParsedHeaderLine parsedHeaderLine) {
        assert(HeaderRecordType.RG.equals(parsedHeaderLine.getHeaderRecordType()));
        parsedHeaderLine.requireTag(SAMReadGroupRecord.READ_GROUP_ID_TAG);
        parsedHeaderLine.requireTag(SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG);
        final SAMReadGroupRecord samReadGroupRecord = new SAMReadGroupRecord(parsedHeaderLine.removeValue(SAMReadGroupRecord.READ_GROUP_ID_TAG));
        transferAttributes(samReadGroupRecord, parsedHeaderLine.mKeyValuePairs);

        // Convert non-String attributes to the appropriate types
        final String predictedMedianInsertSize =
                (String)samReadGroupRecord.getAttribute(SAMReadGroupRecord.PREDICTED_MEDIAN_INSERT_SIZE_TAG);
        if (predictedMedianInsertSize != null) {
            try {
                samReadGroupRecord.setAttribute(SAMReadGroupRecord.PREDICTED_MEDIAN_INSERT_SIZE_TAG,
                    Integer.parseInt(predictedMedianInsertSize));
            } catch (NumberFormatException e) {
                throw new SAMFormatException(SAMReadGroupRecord.PREDICTED_MEDIAN_INSERT_SIZE_TAG +
                        " is not numeric: " + predictedMedianInsertSize, e);
            }
        }

        final String dateRunProduced = (String)samReadGroupRecord.getAttribute(SAMReadGroupRecord.DATE_RUN_PRODUCED_TAG);
        if (dateRunProduced != null) {
            Object date;
            try {
                date = mTagCodec.decodeDate(dateRunProduced);
            } catch (DateParser.InvalidDateException e) {
                switch (validationStringency) {
                    case LENIENT:
                        System.err.println("Ignored error attempting to parse ISO-8601 date string for RG:DT tag: " + dateRunProduced);
                        date = dateRunProduced;
                        break;
                    case SILENT:
                        date = dateRunProduced;
                        break;
                    case STRICT:
                        throw e;
                    default:
                        throw new RuntimeException("Unrecognized validation stringency");
                }
            }
            samReadGroupRecord.setAttribute(SAMReadGroupRecord.DATE_RUN_PRODUCED_TAG, date);
        }

        readGroups.add(samReadGroupRecord);
    }

    private void parseSQLine(final ParsedHeaderLine parsedHeaderLine) {
        assert(HeaderRecordType.SQ.equals(parsedHeaderLine.getHeaderRecordType()));
        parsedHeaderLine.requireTag(SAMSequenceRecord.SEQUENCE_NAME_TAG);
        parsedHeaderLine.requireTag(SAMSequenceRecord.SEQUENCE_LENGTH_TAG);
        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(parsedHeaderLine.removeValue(SAMSequenceRecord.SEQUENCE_NAME_TAG),
                Integer.parseInt(parsedHeaderLine.removeValue(SAMSequenceRecord.SEQUENCE_LENGTH_TAG)));
        transferAttributes(samSequenceRecord, parsedHeaderLine.mKeyValuePairs);
        sequences.add(samSequenceRecord);
    }

    private void parseHDLine(final ParsedHeaderLine parsedHeaderLine) {
        assert(HeaderRecordType.HD.equals(parsedHeaderLine.getHeaderRecordType()));
        parsedHeaderLine.requireTag(SAMFileHeader.VERSION_TAG);
        transferAttributes(mFileHeader, parsedHeaderLine.mKeyValuePairs);
    }

    private RuntimeException reportErrorParsingLine(final String reason) {
        String fileMessage = "";
        if (mFile != null) {
            fileMessage = "File " + mFile + "; ";
        }
        return new SAMFormatException("Error parsing text SAM file. " + reason + "; " + fileMessage +
                "Line " + mReader.getLineNumber() + "\nLine: " + mCurrentLine);
    }

    private enum HeaderRecordType {
        HD, SQ, RG, PG
    }

    private class ParsedHeaderLine {
        private final HeaderRecordType mHeaderRecordType;
        private final Map<String, String> mKeyValuePairs = new HashMap<String, String>();

        ParsedHeaderLine(final String line) {
            assert(line.startsWith(HEADER_LINE_START));
            final String[] fields = line.split(FIELD_SEPARATOR);
            try {
                mHeaderRecordType = HeaderRecordType.valueOf(fields[0].substring(1));
            } catch (IllegalArgumentException e) {
                throw reportErrorParsingLine("Unrecognized header record type");
            }
            for (int i = 1; i < fields.length; ++i) {
                final String[] keyAndValue = fields[i].split(TAG_KEY_VALUE_SEPARATOR, 2);
                if (keyAndValue.length != 2) {
                    throw reportErrorParsingLine("Problem parsing " + HEADER_LINE_START + mHeaderRecordType +
                            " key:value pair");
                }
                mKeyValuePairs.put(keyAndValue[0], keyAndValue[1]);
            }
        }

        void requireTag(final String tag) {
            if (!mKeyValuePairs.containsKey(tag)) {
                throw reportErrorParsingLine(HEADER_LINE_START + mHeaderRecordType + " line missing " + tag + " tag");
            }
        }

        public HeaderRecordType getHeaderRecordType() {
            return mHeaderRecordType;
        }

        boolean containsKey(final String key) {
            return mKeyValuePairs.containsKey(key);
        }

        String getValue(final String key) {
            return mKeyValuePairs.get(key);
        }

        String removeValue(final String key) {
            final String ret = mKeyValuePairs.get(key);
            mKeyValuePairs.remove(key);
            return ret;
        }

    }

    /**
     * Convert SAMFileHeader from in-memory representation to text representation.
     * @param writer where to write the header text.
     * @param header object to be converted to text.
     */
    public void encode(final Writer writer, final SAMFileHeader header) {
        mFileHeader = header;
        this.writer = new BufferedWriter(writer);
        writeHDLine();
        for (final SAMSequenceRecord sequenceRecord: header.getSequenceDictionary().getSequences()) {
            writeSQLine(sequenceRecord);
        }

        for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
            writeRGLine(readGroup);
        }
        for (final SAMProgramRecord programRecord : header.getProgramRecords()) {
            writePGLine(programRecord);
        }
        try {
            this.writer.flush();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private void println(final String s) {
        try {
            writer.append(s);
            writer.append("\n");
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private void writePGLine(final SAMProgramRecord programRecord) {
        if (programRecord == null) {
            return;
        }
        final String[] fields = new String[2 + programRecord.getAttributes().size()];
        fields[0] = HEADER_LINE_START + HeaderRecordType.PG;
        fields[1] = SAMProgramRecord.PROGRAM_GROUP_ID_TAG + TAG_KEY_VALUE_SEPARATOR + programRecord.getProgramGroupId();
        encodeTags(programRecord, fields, 2);
        println(StringUtil.join(FIELD_SEPARATOR, fields));
    }

    private void writeRGLine(final SAMReadGroupRecord readGroup) {
        final String[] fields = new String[2 + readGroup.getAttributes().size()];
        fields[0] = HEADER_LINE_START + HeaderRecordType.RG;
        fields[1] = SAMReadGroupRecord.READ_GROUP_ID_TAG + TAG_KEY_VALUE_SEPARATOR + readGroup.getReadGroupId();
        encodeTags(readGroup, fields, 2);
        println(StringUtil.join(FIELD_SEPARATOR, fields));
    }

    private void writeHDLine() {
        final String[] fields = new String[1 + mFileHeader.getAttributes().size()];
        fields[0] = HEADER_LINE_START + HeaderRecordType.HD;
        encodeTags(mFileHeader, fields, 1);
        println(StringUtil.join(FIELD_SEPARATOR, fields));
    }

    private void writeSQLine(final SAMSequenceRecord sequenceRecord) {
        final int numAttributes =sequenceRecord.getAttributes() != null ? sequenceRecord.getAttributes().size() : 0;
        final String[] fields = new String[3 + numAttributes];
        fields[0] = HEADER_LINE_START + HeaderRecordType.SQ;
        fields[1] = SAMSequenceRecord.SEQUENCE_NAME_TAG + TAG_KEY_VALUE_SEPARATOR + sequenceRecord.getSequenceName();
        fields[2] = SAMSequenceRecord.SEQUENCE_LENGTH_TAG + TAG_KEY_VALUE_SEPARATOR + Integer.toString(sequenceRecord.getSequenceLength());
        encodeTags(sequenceRecord, fields, 3);
        println(StringUtil.join(FIELD_SEPARATOR, fields));
    }

    /**
     * Encode all the attributes in the given object as text
     * @param rec object containing attributes, and knowledge of which are standard tags
     * @param fields where to put the text representation of the tags.  Must be big enough to hold all tags.
     * @param offset where to start putting text tag representations.
     */
    private void encodeTags(final AbstractSAMHeaderRecord rec, final String[] fields, int offset) {
        for (final Map.Entry<String, Object> entry: rec.getAttributes()) {
            final String textTagRepresentation;
            if (rec.getStandardTags().contains(entry.getKey())) {
                textTagRepresentation = mTagCodec.encodeUntypedTag(entry.getKey(), entry.getValue());
            } else {
                textTagRepresentation = mTagCodec.encode(entry.getKey(), entry.getValue());
            }
            fields[offset++] = textTagRepresentation;
        }
    }

    public void setValidationStringency(final SAMFileReader.ValidationStringency validationStringency) {
        if (validationStringency == null) {
            throw new IllegalArgumentException("null validationStringency not allowed");
        }
        this.validationStringency = validationStringency;
    }
}
