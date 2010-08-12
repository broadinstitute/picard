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
import net.sf.samtools.util.StringUtil;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.List;

/**
 * Converter between disk and in-memory representation of a SAMRecord tag.
 */
class BinaryTagCodec {
    // Size in bytes of the fixed part of the disk representation of a tag,
    // i.e. the number of bytes occupied by the tag name and tag type fields.
    private static final int FIXED_TAG_SIZE = 3;

    // Integers are stored in the smallest size that will hold them.
    private static final long MAX_INT = Integer.MAX_VALUE;
    private static final long MAX_UINT = MAX_INT * 2 + 1;
    private static final long MAX_SHORT = Short.MAX_VALUE;
    private static final long MAX_USHORT = MAX_SHORT * 2 + 1;
    private static final long MAX_BYTE = Byte.MAX_VALUE;
    private static final long MAX_UBYTE = MAX_BYTE * 2 + 1;

    // Source or sink for disk representation.
    final BinaryCodec binaryCodec;

    /**
     * For writing tags.
     * For reading tags, a BinaryCodec is not used.  See readTags() below.
     * @param binaryCodec where to write the file rep of the tags
     */
    BinaryTagCodec(final BinaryCodec binaryCodec) {
        this.binaryCodec = binaryCodec;
    }

    /**
     * @param attributeValue In-memory representation of a tag value.
     * @return Size in bytes to store the value on disk.
     */
    private static int getBinaryValueSize(final Object attributeValue) {
        switch (getTagValueType(attributeValue)) {
            case 'Z':
                return ((String)attributeValue).length() + 1;
            case 'A':
                return 1;
            case 'I':
            case 'i':
                return 4;
            case 's':
            case 'S':
                return 2;
            case 'c':
            case 'C':
                return 1;
            case 'f':
                return 4;
            case 'H':
                final byte[] byteArray = (byte[])attributeValue;
                return byteArray.length * 2 + 1;
            default:
                throw new IllegalArgumentException("When writing BAM, unrecognized tag type " +
                        attributeValue.getClass().getName());
        }
    }

    /**
     * @param value In-memory representation of a tag value.
     * @return Size in bytes to store the tag name, tag type and tag value on disk.
     */
    static int getTagSize(final Object value) {
        return FIXED_TAG_SIZE + getBinaryValueSize(value);
    }

    /**
     * @param value In-memory representation of a tag value.
     * @return One-character disk representation of tag type.
     */
    static char getTagValueType(final Object value) {
        if (value instanceof String) {
            return 'Z';
        } else if (value instanceof Character) {
            return 'A';
        } else if (value instanceof Float) {
            return 'f';
        } else if (value instanceof Number) {
            if (!(value instanceof Byte || value instanceof Short || value instanceof Integer || value instanceof Long)) {
                throw new IllegalArgumentException("Unrecognized tag type " + value.getClass().getName());
            }
            return getIntegerType(((Number)value).longValue());
        } else if (value instanceof byte[]) {
            return 'H';
        } else {
            throw new IllegalArgumentException("When writing BAM, unrecognized tag type " +
                    value.getClass().getName());
        }
    }

    /**
     * @param val Integer tag value.
     * @return Tag type corresponding to the smallest integer type that will hold the given value.
     */
    static private char getIntegerType(final long val) {
        if (val > MAX_UINT) {
            throw new IllegalArgumentException("Integer attribute value too large to be encoded in BAM");
        }
        if (val > MAX_INT) {
            return 'I';
        }
        if (val > MAX_USHORT) {
            return 'i';
        }
        if (val > MAX_SHORT) {
            return 'S';
        }
        if (val > MAX_UBYTE) {
            return 's';
        }
        if (val > MAX_BYTE) {
            return 'C';
        }
        if (val >= Byte.MIN_VALUE) {
            return 'c';
        }
        if (val >= Short.MIN_VALUE) {
            return 's';
        }
        if (val >= Integer.MIN_VALUE) {
            return 'i';
        }
        throw new IllegalArgumentException("Integer attribute value too negative to be encoded in BAM");
    }

    /**
     * Write the given tag name and value to disk.
     */
    void writeTag(final short tag, final Object value) {
        binaryCodec.writeShort(tag);
        final char tagValueType = getTagValueType(value);
        binaryCodec.writeByte(tagValueType);

        switch (tagValueType) {
            case 'Z':
                binaryCodec.writeString((String)value, false, true);
                break;
            case 'A':
                binaryCodec.writeByte(((Character)value));
                break;
            case 'I':
                binaryCodec.writeUInt((Long)value);
                break;
            case 'i':
                binaryCodec.writeInt(((Number)value).intValue());
                break;
            case 's':
                binaryCodec.writeShort(((Number)value).shortValue());
                break;
            case 'S':
                binaryCodec.writeUShort(((Number)value).intValue());
                break;
            case 'c':
                binaryCodec.writeByte(((Number)value).byteValue());
                break;
            case 'C':
                binaryCodec.writeUByte(((Integer)value).shortValue());
                break;
            case 'f':
                binaryCodec.writeFloat((Float)value);
                break;
            case 'H':
                final byte[] byteArray = (byte[])value;
                binaryCodec.writeString(StringUtil.bytesToHexString(byteArray), false, true);
                break;
            default:
                throw new IllegalArgumentException("When writing BAM, unrecognized tag type " +
                        value.getClass().getName());
        }
    }

    /**
     * Convert tags from little-endian disk representation to in-memory representation.
     * @param tagCollection Append converted tags to this collection.
     * @param binaryRep Byte buffer containing file representation of tags.
     * @param offset Where in binaryRep tags start.
     * @param length How many bytes in binaryRep are tag storage.
     */
    static SAMBinaryTagAndValue readTags(SAMBinaryTagAndValue tagCollection, final byte[] binaryRep, final int offset,
                         final int length, final SAMFileReader.ValidationStringency validationStringency) {
        final ByteBuffer byteBuffer = ByteBuffer.wrap(binaryRep, offset, length);
        byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        while (byteBuffer.hasRemaining()) {
            final short tag = byteBuffer.getShort();
            final byte tagType = byteBuffer.get();
            final Object value = readValue(tagType, byteBuffer, validationStringency);
            final SAMBinaryTagAndValue newHead = new SAMBinaryTagAndValue(tag, value);
            newHead.setNext(tagCollection);
            tagCollection = newHead;
        }
        return tagCollection;
    }

    /**
     * Read value of specified type.
     * @param tagType What type to read.
     * @param byteBuffer Little-ending byte buffer to read value from.
     * @return Value in in-memory Object form.
     */
    private static  Object readValue(final byte tagType, final ByteBuffer byteBuffer,
                                     final SAMFileReader.ValidationStringency validationStringency) {
        switch (tagType) {
            case 'Z':
                return readNullTerminatedString(byteBuffer);
            case 'A':
                return (char)byteBuffer.get();
            case 'I':
                final long val = byteBuffer.getInt() & 0xffffffffL;
                if (val <= Integer.MAX_VALUE) {
                    return (int)val;
                }
                SAMUtils.processValidationError(new SAMValidationError(SAMValidationError.Type.TAG_VALUE_TOO_LARGE,
                        "Tag value " + val + " too large to store as signed integer.", null), validationStringency);
                // convert to unsigned int stored in a long
                return val;
            case 'i':
                return byteBuffer.getInt();
            case 's':
                return (int)byteBuffer.getShort();
            case 'S':
                // Convert to unsigned short stored in an int
                return byteBuffer.getShort() & 0xffff;
            case 'c':
                return (int)byteBuffer.get();
            case 'C':
                // Convert to unsigned byte stored in an int
                return (int)byteBuffer.get() & 0xff;
            case 'f':
                return byteBuffer.getFloat();
            case 'H':
                final String hexRep = readNullTerminatedString(byteBuffer);
                return StringUtil.hexStringToBytes(hexRep);
            default:
                throw new SAMFormatException("Unrecognized tag type: " + (char)tagType);
        }
    }

    private static String readNullTerminatedString(final ByteBuffer byteBuffer) {
        // Count the number of bytes in the string
        byteBuffer.mark();
        final int startPosition = byteBuffer.position();
        while (byteBuffer.get() != 0) {}
        final int endPosition = byteBuffer.position();

        // Don't count null terminator
        final byte[] buf = new byte[endPosition - startPosition - 1];
        // Go back to the start of the string and read out the bytes
        byteBuffer.reset();
        byteBuffer.get(buf);
        // Skip over the null terminator
        byteBuffer.get();
        return StringUtil.bytesToString(buf);
    }

}
