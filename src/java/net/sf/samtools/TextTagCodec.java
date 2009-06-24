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

import net.sf.samtools.util.Iso8601Date;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.util.DateParser;

import java.util.Map;
import java.text.DateFormat;
import java.text.ParseException;

/**
 * Converter between SAM text representation of a tag, and in-memory Object representation.
 * Note that this class is not thread-safe, in that some local variables have been made into instance
 * variables in order to reduce object creation, but it should not ever be the case that the same
 * instance is used in multiple threads.
 */
class TextTagCodec {
    private static final int NUM_TAG_FIELDS = 3;

    /**
     * This is really a local variable of decode(), but allocated here to reduce allocations.
     */
    private final String[] fields = new String[NUM_TAG_FIELDS];

    /**
     * This is really a local variable of decodeTypeAndValue(), but allocated here to reduce allocations.
     */
    private final String[] typeAndValueFields = new String[NUM_TAG_FIELDS - 1];

    /**
     * Convert in-memory representation of tag to SAM text representation.
     * @param tagName Two-character tag name.
     * @param value Tag value as approriate Object subclass.
     * @return SAM text String representation, i.e. name:type:value
     */
    String encode(final String tagName, Object value) {
        final StringBuilder sb = new StringBuilder(tagName);
        sb.append(':');
        char tagType = BinaryTagCodec.getTagValueType(value);
        switch (tagType) {
            case 'c':
            case 'C':
            case 's':
            case 'S':
            case 'I':
                tagType = 'i';
        }
        if (tagType == 'H') {
            value = SAMUtils.bytesToHexString((byte[])value);
        } else if (tagType == 'i') {
            final long longVal = ((Number) value).longValue();
            if (longVal > Integer.MAX_VALUE || longVal < Integer.MIN_VALUE) {
                throw new SAMFormatException("Value for tag " + tagName + " cannot be stored in an Integer: " + longVal);
            }
        }
        sb.append(tagType);
        sb.append(':');
        sb.append(value.toString());
        return sb.toString();
    }

    /**
     * Encode a standard tag, which should not have a type field.
     * @param tagName 2-character String.
     * @param value Not necessarily a String.  Some of these are integers but the type is implied by
     * the tagName.  Converted to String with toString().
     * @return Colon-separated text representation suitable for a SAM header, i.e. name:value.
     */
    String encodeUntypedTag(final String tagName, final Object value) {
        final StringBuilder sb = new StringBuilder(tagName);
        sb.append(':');
        sb.append(value.toString());
        return sb.toString();
    }

    /**
     * Convert typed tag in SAM text format (name:type:value) into tag name and Object value representation.
     * @param tag SAM text format name:type:value tag.
     * @return Tag name as 2-character String, and tag value in appropriate class based on tag type.
     */
    Map.Entry<String, Object> decode(final String tag) {
        final int numFields = StringUtil.split(tag, fields, ':');
        if (numFields != TextTagCodec.NUM_TAG_FIELDS) {
            throw new SAMFormatException("Not enough fields in tag '" + tag + "'");
        }
        final String key = fields[0];
        final String type = fields[1];
        final String stringVal = fields[2];
        final Object val = convertStringToObject(type, stringVal);
        return new Map.Entry<String, Object>() {
            public String getKey() {
                return key;
            }

            public Object getValue() {
                return val;
            }

            public Object setValue(final Object o) {
                throw new UnsupportedOperationException();
            }
        };
    }

    /**
     * Similar to decode() method above, but the tag name has already been stripped off.
     * @param typeAndValue type:string-value, or, for backward-compatibility, just string-value.
     * @return Value converted into the appropriate type.
     */
    Object decodeTypeAndValue(final String typeAndValue) {
        final int numFields = StringUtil.split(typeAndValue,  typeAndValueFields, ':');
        if (numFields == 1) {
            // For backward compatibility, if no colon, treat as String type
            return typeAndValue;
        }
        return convertStringToObject(typeAndValueFields[0], typeAndValueFields[1]);
    }

    private Object convertStringToObject(final String type, final String stringVal) {
        final Object val;
        if (type.equals("Z")) {
            val = stringVal;
        } else if (type.equals("A")) {
            if (stringVal.length() != 1) {
                throw new SAMFormatException("Tag of type A should have a single-character value");
            }
            val = stringVal.charAt(0);
        } else if (type.equals("i")) {
            try {
                val = new Integer(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type i should have signed decimal value");
            }
        } else if (type.equals("f")) {
            try {
                val = new Float(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type f should have single-precision floating point value");
            }
        } else if (type.equals("H")) {
            try {
                val = SAMUtils.hexStringToBytes(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type H should have valid hex string with even number of digits");
            }
        } else {
            throw new SAMFormatException("Unrecognized tag type: " + type);
        }
        return val;
    }

    Iso8601Date decodeDate(final String dateStr) {
        try {
            return new Iso8601Date(dateStr);
        } catch (DateParser.InvalidDateException ex) {
            try {
                return new Iso8601Date(DateFormat.getDateTimeInstance().parse(dateStr));
            } catch (ParseException e) {
                throw new DateParser.InvalidDateException("Could not parse as date: " + dateStr, e);
            }
        }
    }
}
