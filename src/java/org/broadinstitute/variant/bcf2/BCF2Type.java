/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.bcf2;

import com.google.java.contract.Requires;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.EnumSet;

/**
 * BCF2 types and associated information
 *
 * @author depristo
 * @since 05/12
 */
public enum BCF2Type {
    // the actual values themselves
    MISSING(0, 0, 0x00) {
        @Override public int read(final InputStream in) throws IOException {
            throw new IllegalArgumentException("Cannot read MISSING type");
        }
        @Override public void write(final int value, final OutputStream out) throws IOException {
            throw new IllegalArgumentException("Cannot write MISSING type");
        }
    },

    INT8 (1, 1, 0xFFFFFF80,        -127,        127) {
        @Override
        public int read(final InputStream in) throws IOException {
            return BCF2Utils.readByte(in);
        }

        @Override
        public void write(final int value, final OutputStream out) throws IOException {
            out.write(0xFF & value);   // TODO -- do we need this operation?
        }
    },

    INT16(2, 2, 0xFFFF8000,      -32767,      32767) {
        @Override
        public int read(final InputStream in) throws IOException {
            final int b2 = BCF2Utils.readByte(in) & 0xFF;
            final int b1 = BCF2Utils.readByte(in) & 0xFF;
            return (short)((b1 << 8) | b2);
        }

        @Override
        public void write(final int value, final OutputStream out) throws IOException {
            // TODO -- optimization -- should we put this in a local buffer?
            out.write((0x00FF & value));
            out.write((0xFF00 & value) >> 8);
        }
    },

    INT32(3, 4, 0x80000000, -2147483647, 2147483647) {
        @Override
        public int read(final InputStream in) throws IOException {
            final int b4 = BCF2Utils.readByte(in) & 0xFF;
            final int b3 = BCF2Utils.readByte(in) & 0xFF;
            final int b2 = BCF2Utils.readByte(in) & 0xFF;
            final int b1 = BCF2Utils.readByte(in) & 0xFF;
            return (int)(b1 << 24 | b2 << 16 | b3 << 8 | b4);
        }

        @Override
        public void write(final int value, final OutputStream out) throws IOException {
            out.write((0x000000FF & value));
            out.write((0x0000FF00 & value) >> 8);
            out.write((0x00FF0000 & value) >> 16);
            out.write((0xFF000000 & value) >> 24);
        }
    },

    FLOAT(5, 4, 0x7F800001) {
        @Override
        public int read(final InputStream in) throws IOException {
            return INT32.read(in);
        }

        @Override
        public void write(final int value, final OutputStream out) throws IOException {
            INT32.write(value, out);
        }
    },

    CHAR (7, 1, 0x00000000) {
        @Override
        public int read(final InputStream in) throws IOException {
            return INT8.read(in);
        }

        @Override
        public void write(final int value, final OutputStream out) throws IOException {
            INT8.write(value, out);
        }
    };

    private final int id;
    private final Object missingJavaValue;
    private final int missingBytes;
    private final int sizeInBytes;
    private final long minValue, maxValue;

    BCF2Type(final int id, final int sizeInBytes, final int missingBytes) {
        this(id, sizeInBytes, missingBytes, 0, 0);
    }

    BCF2Type(final int id, final int sizeInBytes, final int missingBytes, final long minValue, final long maxValue) {
        this.id = id;
        this.sizeInBytes = sizeInBytes;
        this.missingJavaValue = null;
        this.missingBytes = missingBytes;
        this.minValue = minValue;
        this.maxValue = maxValue;
    }

    /**
     * How many bytes are used to represent this type on disk?
     * @return
     */
    public int getSizeInBytes() {
        return sizeInBytes;
    }

    /**
     * The ID according to the BCF2 specification
     * @return
     */
    public int getID() { return id; }

    /**
     * Can we encode value v in this type, according to its declared range.
     *
     * Only makes sense for integer values
     *
     * @param v
     * @return
     */
    @Requires("this.isIntegerType()")
    public final boolean withinRange(final long v) { return v >= minValue && v <= maxValue; }

    /**
     * Return the java object (aka null) that is used to represent a missing value for this
     * type in Java
     *
     * @return
     */
    public Object getMissingJavaValue() { return missingJavaValue; }

    /**
     * The bytes (encoded as an int) that are used to represent a missing value
     * for this type in BCF2
     *
     * @return
     */
    public int getMissingBytes() { return missingBytes; }

    /**
     * An enum set of the types that might represent Integer values
     */
    private final static EnumSet<BCF2Type> INTEGERS = EnumSet.of(INT8, INT16, INT32);

    /**
     * @return true if this BCF2Type corresponds to the magic "MISSING" type (0x00)
     */
    public boolean isMissingType() {
        return this == MISSING;
    }

    public boolean isIntegerType() {
        return INTEGERS.contains(this);
    }

    /**
     * Read a value from in stream of this BCF2 type as an int [32 bit] collection of bits
     *
     * For intX and char values this is just the int / byte value of the underlying data represented as a 32 bit int
     * For a char the result must be converted to a char by (char)(byte)(0x0F & value)
     * For doubles it's necessary to convert subsequently this value to a double via Double.bitsToDouble()
     *
     * @param in
     * @return
     * @throws IOException
     */
    @Requires("in != null")
    public int read(final InputStream in) throws IOException {
        throw new IllegalArgumentException("Not implemented");
    }

    @Requires("out != null")
    public void write(final int value, final OutputStream out) throws IOException {
        throw new IllegalArgumentException("Not implemented");
    }
}
