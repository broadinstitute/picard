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

import java.lang.reflect.Array;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Converter between disk and in-memory representation of a SAMRecord tag.
 */
class BinaryTagCodec {
    // Size in bytes of the fixed part of the disk representation of a tag,
    // i.e. the number of bytes occupied by the tag name and tag type fields.
    private static final int FIXED_TAG_SIZE = 3;

    // Size in bytes of the fixed part of the value of a binary array,
    // i.e. the number of bytes occupied by the array type and the array length.
    private static final int FIXED_BINARY_ARRAY_TAG_SIZE = 5;

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
            case 'B':
                final int numElements = Array.getLength(attributeValue);
                final int elementSize;
                if(attributeValue instanceof byte[]) {
                    elementSize = 1;
                } else if(attributeValue instanceof short[]) {
                    elementSize = 2;
                } else if(attributeValue instanceof int[]) {
                    elementSize = 4;
                } else if(attributeValue instanceof float[]) {
                    elementSize = 4;
                } else {
                    throw new IllegalArgumentException("Unsupported array type: " + attributeValue.getClass());
                }
                return numElements * elementSize + FIXED_BINARY_ARRAY_TAG_SIZE;
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
        } /*
           Note that H tag type is never written anymore, because B style is more compact.
           else if (value instanceof byte[]) {
            return 'H';
           }
          */
        else if (value instanceof byte[] || value instanceof short[] || value instanceof int[] || value instanceof float[]) {
            return 'B';
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
    void writeTag(final short tag, final Object value, final boolean isUnsignedArray) {
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
            /*
            Writing H is no longer supported
            case 'H':
                final byte[] byteArray = (byte[])value;
                binaryCodec.writeString(StringUtil.bytesToHexString(byteArray), false, true);
                break;
             */
            case 'B':
                writeArray(value, isUnsignedArray);
                break;
            default:
                throw new IllegalArgumentException("When writing BAM, unrecognized tag type " +
                        value.getClass().getName());
        }
    }

    private void writeArray(final Object value, final boolean isUnsignedArray) {
        if (value instanceof byte[]) {
            binaryCodec.writeByte(isUnsignedArray? 'C': 'c');
            final byte[] array = (byte[]) value;
            binaryCodec.writeInt(array.length);
            for (final byte element: array) binaryCodec.writeByte(element);

        } else if (value instanceof short[]) {
            binaryCodec.writeByte(isUnsignedArray? 'S': 's');
            final short[] array = (short[]) value;
            binaryCodec.writeInt(array.length);
            for (final short element: array) binaryCodec.writeShort(element);

        } else if (value instanceof int[]) {
            binaryCodec.writeByte(isUnsignedArray? 'I': 'i');
            final int[] array = (int[]) value;
            binaryCodec.writeInt(array.length);
            for (final int element: array) binaryCodec.writeInt(element);

        } else if (value instanceof float[]) {
            binaryCodec.writeByte('f');
            final float[] array = (float[]) value;
            binaryCodec.writeInt(array.length);
            for (final float element: array) binaryCodec.writeFloat(element);

        } else throw new SAMException("Unrecognized array value type: " + value.getClass());
    }

    /**
     * Convert tags from little-endian disk representation to in-memory representation.
     * @param binaryRep Byte buffer containing file representation of tags.
     * @param offset Where in binaryRep tags start.
     * @param length How many bytes in binaryRep are tag storage.
     */
    static SAMBinaryTagAndValue readTags(final byte[] binaryRep, final int offset,
                                         final int length, final SAMFileReader.ValidationStringency validationStringency) {
        final ByteBuffer byteBuffer = ByteBuffer.wrap(binaryRep, offset, length);
        byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        SAMBinaryTagAndValue head = null;
        SAMBinaryTagAndValue tail = null;

        while (byteBuffer.hasRemaining()) {
            final short tag = byteBuffer.getShort();
            final byte tagType = byteBuffer.get();
            final SAMBinaryTagAndValue tmp;
            if (tagType != 'B') {
                tmp = new SAMBinaryTagAndValue(tag, readSingleValue(tagType, byteBuffer, validationStringency));
            } else {
                final TagValueAndUnsignedArrayFlag valueAndFlag = readArray(byteBuffer, validationStringency);
                if (valueAndFlag.isUnsignedArray) tmp = new SAMBinaryTagAndUnsignedArrayValue(tag, valueAndFlag.value);
                else tmp = new SAMBinaryTagAndValue(tag, valueAndFlag.value);
            }

            // If samjdk wrote the BAM then the attributes will be in lowest->highest tag order, to inserting at the
            // head each time will be very inefficient. To fix that we check here to see if the tag should go right on
            // the tail and if so stick it there, else insert it through the head.
            if (head == null) {
                head = tmp;
                tail = tmp;
            }
            else if (tmp.tag > tail.tag) {
                tail.insert(tmp);
                tail = tmp;
            }
            else {
                head = head.insert(tmp);
            }
        }

        return head;
    }

    /**
     * Read value of specified non-array type.
     * @param tagType What type to read.
     * @param byteBuffer Little-ending byte buffer to read value from.
     * @return Value in in-memory Object form.
     */
    private static  Object readSingleValue(final byte tagType, final ByteBuffer byteBuffer,
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




    /**
     * Read value of specified type.
     * @param byteBuffer Little-ending byte buffer to read value from.
     * @return CVO containing the value in in-memory Object form, and a flag indicating whether it is unsigned or not.
     */
    private static TagValueAndUnsignedArrayFlag readArray(final ByteBuffer byteBuffer,
                                                          final SAMFileReader.ValidationStringency validationStringency) {
        final byte arrayType = byteBuffer.get();
        final boolean isUnsigned = Character.isUpperCase(arrayType);
        final int length = byteBuffer.getInt();
        final Object value;
        switch (Character.toLowerCase(arrayType)) {
            case 'c': {
                final byte[] array = new byte[length];
                value = array;
                byteBuffer.get(array);
                break;
            }
            case 's': {
                final short[] array = new short[length];
                value = array;
                for (int i = 0; i < length; ++i) {
                    array[i] = byteBuffer.getShort();
                }
                break;
            }

            case 'i': {
                final int[] array = new int[length];
                value = array;
                for (int i = 0; i < length; ++i) {
                    array[i] = byteBuffer.getInt();
                }
                break;
            }

            case 'f': {
                final float[] array = new float[length];
                value = array;
                for (int i = 0; i < length; ++i) {
                    array[i] = byteBuffer.getFloat();
                }
                break;
            }

            default:
                throw new SAMFormatException("Unrecognized tag array type: " + (char)arrayType);
        }
        return new TagValueAndUnsignedArrayFlag(value, isUnsigned);
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
