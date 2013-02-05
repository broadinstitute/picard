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

package org.broadinstitute.variant.variantcontext.writer;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.variant.bcf2.BCF2Type;
import org.broadinstitute.variant.bcf2.BCF2Utils;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * See #BCFWriter for documentation on this classes role in encoding BCF2 files
 *
 * @author Mark DePristo
 * @since 06/12
 */
public final class BCF2Encoder {
    // TODO -- increase default size?
    public static final int WRITE_BUFFER_INITIAL_SIZE = 16384;
    private ByteArrayOutputStream encodeStream = new ByteArrayOutputStream(WRITE_BUFFER_INITIAL_SIZE);

    // --------------------------------------------------------------------------------
    //
    // Functions to return the data being encoded here
    //
    // --------------------------------------------------------------------------------

    @Ensures("result != null")
    public byte[] getRecordBytes() {
        byte[] bytes = encodeStream.toByteArray();
        encodeStream.reset();
        return bytes;
    }

    // --------------------------------------------------------------------------------
    //
    // Writing typed values (have type byte)
    //
    // --------------------------------------------------------------------------------

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeTypedMissing(final BCF2Type type) throws IOException {
        encodeType(0, type);
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeTyped(final Object value, final BCF2Type type) throws IOException {
        if ( value == null )
            encodeTypedMissing(type);
        else {
            switch ( type ) {
                case INT8:
                case INT16:
                case INT32: encodeTypedInt((Integer)value, type); break;
                case FLOAT: encodeTypedFloat((Double) value); break;
                case CHAR:  encodeTypedString((String) value); break;
                default:    throw new IllegalArgumentException("Illegal type encountered " + type);
            }
        }
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeTypedInt(final int v) throws IOException {
        final BCF2Type type = BCF2Utils.determineIntegerType(v);
        encodeTypedInt(v, type);
    }

    @Requires("type.isIntegerType()")
    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeTypedInt(final int v, final BCF2Type type) throws IOException {
        encodeType(1, type);
        encodeRawInt(v, type);
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeTypedString(final String s) throws IOException {
        encodeTypedString(s.getBytes());
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeTypedString(final byte[] s) throws IOException {
        if ( s == null )
            encodeType(0, BCF2Type.CHAR);
        else {
            encodeType(s.length, BCF2Type.CHAR);
            for ( int i = 0; i < s.length; i++ ) {
                encodeRawChar(s[i]);
            }
        }
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeTypedFloat(final double d) throws IOException {
        encodeType(1, BCF2Type.FLOAT);
        encodeRawFloat(d);
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeTyped(List<? extends Object> v, final BCF2Type type) throws IOException {
        if ( type == BCF2Type.CHAR && v.size() != 0 ) {
            final String s = BCF2Utils.collapseStringList((List<String>) v);
            v = stringToBytes(s);
        }

        encodeType(v.size(), type);
        encodeRawValues(v, type);
    }

    // --------------------------------------------------------------------------------
    //
    // Writing raw values (don't have a type byte)
    //
    // --------------------------------------------------------------------------------

    public final <T extends Object> void encodeRawValues(final Collection<T> v, final BCF2Type type) throws IOException {
        for ( final T v1 : v ) {
            encodeRawValue(v1, type);
        }
    }

    public final <T extends Object> void encodeRawValue(final T value, final BCF2Type type) throws IOException {
        try {
            if ( value == type.getMissingJavaValue() )
                encodeRawMissingValue(type);
            else {
                switch (type) {
                    case INT8:
                    case INT16:
                    case INT32: encodeRawBytes((Integer) value, type); break;
                    case FLOAT: encodeRawFloat((Double) value); break;
                    case CHAR:  encodeRawChar((Byte) value); break;
                    default:    throw new IllegalArgumentException("Illegal type encountered " + type);
                }
            }
        } catch ( ClassCastException e ) {
            throw new ClassCastException("BUG: invalid type cast to " + type + " from " + value);
        }
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeRawMissingValue(final BCF2Type type) throws IOException {
        encodeRawBytes(type.getMissingBytes(), type);
    }

    @Requires("size >= 0")
    public final void encodeRawMissingValues(final int size, final BCF2Type type) throws IOException {
        for ( int i = 0; i < size; i++ )
            encodeRawMissingValue(type);
    }

    // --------------------------------------------------------------------------------
    //
    // low-level encoders
    //
    // --------------------------------------------------------------------------------

    public final void encodeRawChar(final byte c) throws IOException {
        encodeStream.write(c);
    }

    public final void encodeRawFloat(final double value) throws IOException {
        encodeRawBytes(Float.floatToIntBits((float) value), BCF2Type.FLOAT);
    }

    @Requires("size >= 0")
    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeType(final int size, final BCF2Type type) throws IOException {
        if ( size <= BCF2Utils.MAX_INLINE_ELEMENTS ) {
            final int typeByte = BCF2Utils.encodeTypeDescriptor(size, type);
            encodeStream.write(typeByte);
        } else {
            final int typeByte = BCF2Utils.encodeTypeDescriptor(BCF2Utils.OVERFLOW_ELEMENT_MARKER, type);
            encodeStream.write(typeByte);
            // write in the overflow size
            encodeTypedInt(size);
        }
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeRawInt(final int value, final BCF2Type type) throws IOException {
        type.write(value, encodeStream);
    }

    @Ensures("encodeStream.size() > old(encodeStream.size())")
    public final void encodeRawBytes(final int value, final BCF2Type type) throws IOException {
        type.write(value, encodeStream);
    }

    // --------------------------------------------------------------------------------
    //
    // utility functions
    //
    // --------------------------------------------------------------------------------

    @Requires({"s != null", "sizeToWrite >= 0"})
    public void encodeRawString(final String s, final int sizeToWrite) throws IOException {
        final byte[] bytes = s.getBytes();
        for ( int i = 0; i < sizeToWrite; i++ )
            if ( i < bytes.length )
                encodeRawChar(bytes[i]);
            else
                encodeRawMissingValue(BCF2Type.CHAR);
    }

    /**
     * Totally generic encoder that examines o, determines the best way to encode it, and encodes it
     *
     * This method is incredibly slow, but it's only used for UnitTests so it doesn't matter
     *
     * @param o
     * @return
     */
    @Requires("o != null")
    public final BCF2Type encode(final Object o) throws IOException {
        if ( o == null ) throw new IllegalArgumentException("Generic encode cannot deal with null values");

        if ( o instanceof List ) {
            final BCF2Type type = determineBCFType(((List) o).get(0));
            encodeTyped((List) o, type);
            return type;
        } else {
            final BCF2Type type = determineBCFType(o);
            encodeTyped(o, type);
            return type;
        }
    }

    @Requires("arg != null")
    private final BCF2Type determineBCFType(final Object arg) {
        final Object toType = arg instanceof List ? ((List)arg).get(0) : arg;

        if ( toType instanceof Integer )
            return BCF2Utils.determineIntegerType((Integer) toType);
        else if ( toType instanceof String )
            return BCF2Type.CHAR;
        else if ( toType instanceof Double )
            return BCF2Type.FLOAT;
        else
            throw new IllegalArgumentException("No native encoding for Object of type " + arg.getClass().getSimpleName());
    }

    private final List<Byte> stringToBytes(final String v) throws IOException {
        if ( v == null || v.equals("") )
            return Collections.emptyList();
        else {
            // TODO -- this needs to be optimized away for efficiency
            final byte[] bytes = v.getBytes();
            final List<Byte> l = new ArrayList<Byte>(bytes.length);
            for ( int i = 0; i < bytes.length; i++) l.add(bytes[i]);
            return l;
        }
    }
}