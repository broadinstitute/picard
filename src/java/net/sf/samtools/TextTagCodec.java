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

import java.lang.reflect.Array;
import java.util.Map;
import java.util.Date;
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
            // H should never happen anymore.
            value = StringUtil.bytesToHexString((byte[])value);
        } else if (tagType == 'B') {
            value = getArrayType(value, false) + "," + encodeArrayValue(value);
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

    private char getArrayType(final Object array, final boolean isUnsigned) {
        final char type;
        final Class<?> componentType = array.getClass().getComponentType();
        if (componentType == Float.TYPE) {
            if (isUnsigned) throw new IllegalArgumentException("float array cannot be unsigned");
            return 'f';
        }
        else if (componentType == Byte.TYPE)    type = 'c';
        else if (componentType == Short.TYPE)   type = 's';
        else if (componentType == Integer.TYPE) type = 'i';
        else throw new IllegalArgumentException("Unrecognized array type " + componentType);
        return (isUnsigned? Character.toUpperCase(type): type);
    }

    private String encodeArrayValue(final Object value) {
        final StringBuilder ret = new StringBuilder(Array.get(value, 0).toString());
        final int length = Array.getLength(value);
        for (int i = 1; i < length; ++i) {
            ret.append(",");
            ret.append(Array.get(value, i).toString());
        }
        return ret.toString();

    }

    private long[] widenToUnsigned(final Object array) {
        final Class<?> componentType = array.getClass().getComponentType();
        final long mask;
        if (componentType == Byte.TYPE)    mask = 0xffL;
        else if (componentType == Short.TYPE)   mask = 0xffffL;
        else if (componentType == Integer.TYPE) mask = 0xffffffffL;
        else throw new IllegalArgumentException("Unrecognized unsigned array type " + componentType);
        final long[] ret = new long[Array.getLength(array)];
        for (int i = 0; i < ret.length; ++i) {
            ret[i] = Array.getLong(array, i) & mask;
        }
        return ret;
    }

    String encodeUnsignedArray(final String tagName, final Object array) {
        if (!array.getClass().isArray()) {
            throw new IllegalArgumentException("Non-array passed to encodeUnsignedArray: " + array.getClass());
        }
        final long[] widened = widenToUnsigned(array);
        return tagName + ":B:" + getArrayType(array, true) + "," + encodeArrayValue(widened);
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
     * If value is an unsigned array, then the value is a TagValueAndUnsignedArrayFlag object.
     */
    Map.Entry<String, Object> decode(final String tag) {
        final int numFields = StringUtil.splitConcatenateExcessTokens(tag, fields, ':');
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

    private Object convertStringToObject(final String type, final String stringVal) {
        if (type.equals("Z")) {
            return stringVal;
        } else if (type.equals("A")) {
            if (stringVal.length() != 1) {
                throw new SAMFormatException("Tag of type A should have a single-character value");
            }
            return stringVal.charAt(0);
        } else if (type.equals("i")) {
            try {
                return new Integer(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type i should have signed decimal value");
            }
        } else if (type.equals("f")) {
            try {
                return new Float(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type f should have single-precision floating point value");
            }
        } else if (type.equals("H")) {
            try {
                return StringUtil.hexStringToBytes(stringVal);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Tag of type H should have valid hex string with even number of digits");
            }
        } else if (type.equals("B")) {
            return covertStringArrayToObject(stringVal);
        } else {
            throw new SAMFormatException("Unrecognized tag type: " + type);
        }
    }

    private Object covertStringArrayToObject(final String stringVal) {
        final String[] elementTypeAndValue = new String[2];
        if (StringUtil.splitConcatenateExcessTokens(stringVal, elementTypeAndValue, ',') != 2) {
            throw new SAMFormatException("Tag of type B should have an element type followed by comma");
        }
        if (elementTypeAndValue[0].length() != 1) {
            throw new SAMFormatException("Unrecognized element type for array tag value: " + elementTypeAndValue[0]);
        }
        final char elementType = elementTypeAndValue[0].charAt(0);
        final String[] stringValues = elementTypeAndValue[1].split(",");
        if (stringValues.length == 0) throw new SAMFormatException("Tag of type B should have at least one element");
        if (elementType == 'f') {
            final float[] ret = new float[stringValues.length];
            for (int i = 0; i < stringValues.length; ++i) {
                try {
                    ret[i] = Float.parseFloat(stringValues[i]);
                } catch (NumberFormatException e) {
                    throw new SAMFormatException("Array tag of type f should have single-precision floating point value");
                }
            }
            return ret;
        }
        long mask = Long.MAX_VALUE;
        long minValue = Long.MAX_VALUE;
        long maxValue = Long.MIN_VALUE;
        final boolean isUnsigned = Character.isUpperCase(elementType);
        switch (Character.toLowerCase(elementType)) {
            case 'c':
                if (isUnsigned) {
                    mask = 0xffL;
                } else {
                    minValue = Byte.MIN_VALUE;
                    maxValue = Byte.MAX_VALUE;
                }
                break;
            case 's':
                if (isUnsigned) {
                    mask = 0xffffL;
                } else {
                    minValue = Short.MIN_VALUE;
                    maxValue = Short.MAX_VALUE;
                }
                break;
            case 'i':
                if (isUnsigned) {
                    mask = 0xffffffffL;
                } else {
                    minValue = Integer.MIN_VALUE;
                    maxValue = Integer.MAX_VALUE;
                }
                break;
            default:
                throw new SAMFormatException("Unrecognized array tag element type: " + elementType);
        }
        if (isUnsigned) {
            minValue = 0;
            maxValue = mask;
        }
        final long[] longValues = new long[stringValues.length];
        for (int i = 0; i < stringValues.length; ++i) {
            final long longValue;
            try {
                longValue = Long.parseLong(stringValues[i]);
            } catch (NumberFormatException e) {
                throw new SAMFormatException("Array tag of type " + elementType + " should have integral value");
            }
            if (longValue < minValue || longValue > maxValue) {
                throw new SAMFormatException("Value for element of array tag of type " + elementType +
                " is out of allowed range: " + longValue);
            }
            longValues[i] = longValue;
        }

        switch (Character.toLowerCase(elementType)) {
            case 'c': {
                final byte[] array = new byte[longValues.length];
                for (int i = 0; i < longValues.length; ++i) array[i] = (byte)longValues[i];
                if (isUnsigned) return new TagValueAndUnsignedArrayFlag(array, true);
                else return array;
            }
            case 's': {
                final short[] array = new short[longValues.length];
                for (int i = 0; i < longValues.length; ++i) array[i] = (short)longValues[i];
                if (isUnsigned) return new TagValueAndUnsignedArrayFlag(array, true);
                else return array;
            }
            case 'i':{
                final int[] array = new int[longValues.length];
                for (int i = 0; i < longValues.length; ++i) array[i] = (int)longValues[i];
                if (isUnsigned) return new TagValueAndUnsignedArrayFlag(array, true);
                else return array;
            }
            default:
                throw new SAMFormatException("Unrecognized array tag element type: " + elementType);
        }
    }

    Iso8601Date decodeDate(final String dateStr) {
        try {
            return new Iso8601Date(dateStr);
        } catch (DateParser.InvalidDateException ex) {
            try {
                return new Iso8601Date(DateFormat.getDateTimeInstance().parse(dateStr));
            } catch (ParseException e) {
                try {
                    return new Iso8601Date(new Date(dateStr));
                } catch (Exception e1) {
                    throw new DateParser.InvalidDateException("Could not parse as date: " + dateStr, e);
                }
            }
        }
    }
}
