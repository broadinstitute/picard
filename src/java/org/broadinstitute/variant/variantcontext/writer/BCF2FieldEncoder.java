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
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.broadinstitute.variant.bcf2.BCF2Type;
import org.broadinstitute.variant.bcf2.BCF2Utils;
import org.broadinstitute.variant.vcf.VCFCompoundHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * See #BCFWriter for documentation on this classes role in encoding BCF2 files
 *
 * @author Mark DePristo
 * @since 06/12
 */
@Invariant({
        "headerLine != null",
        "dictionaryOffsetType.isIntegerType()",
        "dictionaryOffset >= 0"
})
public abstract class BCF2FieldEncoder {
    /**
     * The header line describing the field we will encode values of
     */
    final VCFCompoundHeaderLine headerLine;

    /**
     * The BCF2 type we'll use to encoder this field, if it can be determined statically.
     * If not, this variable must be null
     */
    final BCF2Type staticType;

    /**
     * The integer offset into the strings map of the BCF2 file corresponding to this
     * field.
     */
    final int dictionaryOffset;

    /**
     * The integer type we use to encode our dictionary offset in the BCF2 file
     */
    final BCF2Type dictionaryOffsetType;

    // ----------------------------------------------------------------------
    //
    // Constructor
    //
    // ----------------------------------------------------------------------

    @Requires({"headerLine != null", "dict != null"})
    private BCF2FieldEncoder(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict, final BCF2Type staticType) {
        this.headerLine = headerLine;
        this.staticType = staticType;

        final Integer offset = dict.get(getField());
        if ( offset == null ) throw new IllegalStateException("Format error: could not find string " + getField() + " in header as required by BCF");
        this.dictionaryOffset = offset;
        dictionaryOffsetType = BCF2Utils.determineIntegerType(offset);
    }

    // ----------------------------------------------------------------------
    //
    // Basic accessors
    //
    // ----------------------------------------------------------------------

    @Ensures("result != null")
    public final String getField() { return headerLine.getID(); }

    /**
     * Write the field key (dictionary offset and type) into the BCF2Encoder stream
     *
     * @param encoder where we write our dictionary offset
     * @throws IOException
     */
    @Requires("encoder != null")
    public final void writeFieldKey(final BCF2Encoder encoder) throws IOException {
        encoder.encodeTypedInt(dictionaryOffset, dictionaryOffsetType);
    }

    @Override
    public String toString() {
        return "BCF2FieldEncoder for " + getField() + " with count " + getCountType() + " encoded with " + getClass().getSimpleName();
    }

    // ----------------------------------------------------------------------
    //
    // methods to determine the number of encoded elements
    //
    // ----------------------------------------------------------------------

    @Ensures("result != null")
    protected final VCFHeaderLineCount getCountType() {
        return headerLine.getCountType();
    }

    /**
     * True if this field has a constant, fixed number of elements (such as 1 for an atomic integer)
     *
     * @return
     */
    @Ensures("result != (hasValueDeterminedNumElements() || hasContextDeterminedNumElements())")
    public boolean hasConstantNumElements() {
        return getCountType() == VCFHeaderLineCount.INTEGER;
    }

    /**
     * True if the only way to determine how many elements this field contains is by
     * inspecting the actual value directly, such as when the number of elements
     * is a variable length list per site or per genotype.
     * @return
     */
    @Ensures("result != (hasConstantNumElements() || hasContextDeterminedNumElements())")
    public boolean hasValueDeterminedNumElements() {
        return getCountType() == VCFHeaderLineCount.UNBOUNDED;
    }

    /**
     * True if this field has a non-fixed number of elements that depends only on the properties
     * of the current VariantContext, such as one value per Allele or per genotype configuration.
     *
     * @return
     */
    @Ensures("result != (hasValueDeterminedNumElements() || hasConstantNumElements())")
    public boolean hasContextDeterminedNumElements() {
        return ! hasConstantNumElements() && ! hasValueDeterminedNumElements();
    }

    /**
     * Get the number of elements, assuming this field has a constant number of elements.
     * @return
     */
    @Requires("hasConstantNumElements()")
    @Ensures("result >= 0")
    public int numElements() {
        return headerLine.getCount();
    }

    /**
     * Get the number of elements by looking at the actual value provided
     * @return
     */
    @Requires("hasValueDeterminedNumElements()")
    @Ensures("result >= 0")
    public int numElements(final Object value) {
        return numElementsFromValue(value);
    }

    /**
     * Get the number of elements, assuming this field has context-determined number of elements.
     * @return
     */
    @Requires("hasContextDeterminedNumElements()")
    @Ensures("result >= 0")
    public int numElements(final VariantContext vc) {
        return headerLine.getCount(vc);
    }

    /**
     * A convenience access for the number of elements, returning
     * the number of encoded elements, either from the fixed number
     * it has, from the VC, or from the value itself.
     * @param vc
     * @param value
     * @return
     */
    @Ensures("result >= 0")
    public final int numElements(final VariantContext vc, final Object value) {
        if ( hasConstantNumElements() ) return numElements();
        else if ( hasContextDeterminedNumElements() ) return numElements(vc);
        else return numElements(value);
    }

    /**
     * Given a value, return the number of elements we will encode for it.
     *
     * Assumes the value is encoded as a List
     *
     * @param value
     * @return
     */
    @Requires("hasValueDeterminedNumElements()")
    @Ensures("result >= 0")
    protected int numElementsFromValue(final Object value) {
        if ( value == null ) return 0;
        else if ( value instanceof List ) return ((List) value).size();
        else return 1;
    }

    // ----------------------------------------------------------------------
    //
    // methods to determine the BCF2 type of the encoded values
    //
    // ----------------------------------------------------------------------

    /**
     * Is the BCF2 type of this field static, or does it have to be determine from
     * the actual field value itself?
     * @return
     */
    @Ensures("result || isDynamicallyTyped()")
    public final boolean isStaticallyTyped() { return ! isDynamicallyTyped(); }

    /**
     * Is the BCF2 type of this field static, or does it have to be determine from
     * the actual field value itself?
     * @return
     */
    @Ensures("result || isStaticallyTyped()")
    public final boolean isDynamicallyTyped() { return staticType == null; }

    /**
     * Get the BCF2 type for this field, either from the static type of the
     * field itself or by inspecting the value itself.
     *
     * @return
     */
    public final BCF2Type getType(final Object value) {
        return isDynamicallyTyped() ? getDynamicType(value) : getStaticType();
    }

    @Requires("isStaticallyTyped()")
    @Ensures("result != null")
    public final BCF2Type getStaticType() {
        return staticType;
    }

    @Requires("isDynamicallyTyped()")
    @Ensures("result != null")
    public BCF2Type getDynamicType(final Object value) {
        throw new IllegalStateException("BUG: cannot get dynamic type for statically typed BCF2 field " + getField());
    }

    // ----------------------------------------------------------------------
    //
    // methods to encode values, including the key abstract method
    //
    // ----------------------------------------------------------------------

    /**
     * Key abstract method that should encode a value of the given type into the encoder.
     *
     * Value will be of a type appropriate to the underlying encoder.  If the genotype field is represented as
     * an int[], this will be value, and the encoder needs to handle encoding all of the values in the int[].
     *
     * The argument should be used, not the getType() method in the superclass as an outer loop might have
     * decided a more general type (int16) to use, even through this encoder could have been done with int8.
     *
     * If minValues > 0, then encodeValue must write in at least minValues items from value.  If value is atomic,
     * this means that minValues - 1 MISSING values should be added to the encoder.  If minValues is a collection
     * type (int[]) then minValues - values.length should be added.  This argument is intended to handle padding
     * of values in genotype fields.
     *
     * @param encoder
     * @param value
     * @param type
     * @param minValues
     * @throws IOException
     */
    @Requires({"encoder != null", "isDynamicallyTyped() || type == getStaticType()", "minValues >= 0"})
    public abstract void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException;

    // ----------------------------------------------------------------------
    //
    // Subclass to encode Strings
    //
    // ----------------------------------------------------------------------

    public static class StringOrCharacter extends BCF2FieldEncoder {
        public StringOrCharacter(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, BCF2Type.CHAR);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            final String s = javaStringToBCF2String(value);
            encoder.encodeRawString(s, Math.max(s.length(), minValues));
        }

        //
        // Regardless of what the header says, BCF2 strings and characters are always encoded
        // as arrays of CHAR type, which has a variable number of elements depending on the
        // exact string being encoded
        //
        @Override public boolean hasConstantNumElements()          { return false; }
        @Override public boolean hasContextDeterminedNumElements() { return false; }
        @Override public boolean hasValueDeterminedNumElements()   { return true; }
        @Override protected int numElementsFromValue(final Object value) {
            return value == null ? 0 : javaStringToBCF2String(value).length();
        }

        /**
         * Recode the incoming object to a String, compacting it into a
         * BCF2 string if the value is a list.
         *
         * @param value a String or List<String> to encode, or null
         * @return a non-null string to encode
         */
        @Ensures("result != null")
        private String javaStringToBCF2String(final Object value) {
            if ( value == null )
                return "";
            else if (value instanceof List) {
                final List<String> l = (List<String>)value;
                if ( l.isEmpty() ) return "";
                else return BCF2Utils.collapseStringList(l);
            } else
                return (String)value;
        }
    }

    // ----------------------------------------------------------------------
    //
    // Subclass to encode FLAG
    //
    // ----------------------------------------------------------------------

    public static class Flag extends BCF2FieldEncoder {
        public Flag(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, BCF2Type.INT8);
            if ( ! headerLine.isFixedCount() || headerLine.getCount() != 0 )
                throw new IllegalStateException("Flag encoder only supports atomic flags for field " + getField());
        }

        @Override
        public int numElements() {
            return 1; // the header says 0 but we will write 1 value
        }

        @Override
        @Requires({"minValues <= 1", "value != null", "value instanceof Boolean", "((Boolean)value) == true"})
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            encoder.encodeRawBytes(1, getStaticType());
        }
    }

    // ----------------------------------------------------------------------
    //
    // Subclass to encode FLOAT
    //
    // ----------------------------------------------------------------------

    public static class Float extends BCF2FieldEncoder {
        final boolean isAtomic;

        public Float(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, BCF2Type.FLOAT);
            isAtomic = hasConstantNumElements() && numElements() == 1;
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            int count = 0;
            // TODO -- can be restructured to avoid toList operation
            if ( isAtomic ) {
                // fast path for fields with 1 fixed float value
                if ( value != null ) {
                    encoder.encodeRawFloat((Double)value);
                    count++;
                }
            } else {
                // handle generic case
                final List<Double> doubles = toList(Double.class, value);
                for ( final Double d : doubles ) {
                    if ( d != null ) { // necessary because .,. => [null, null] in VC
                        encoder.encodeRawFloat(d);
                        count++;
                    }
                }
            }
            for ( ; count < minValues; count++ ) encoder.encodeRawMissingValue(type);
        }
    }

    // ----------------------------------------------------------------------
    //
    // Subclass to encode int[]
    //
    // ----------------------------------------------------------------------

    public static class IntArray extends BCF2FieldEncoder {
        public IntArray(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, null);
        }

        @Override
        protected int numElementsFromValue(final Object value) {
            return value == null ? 0 : ((int[])value).length;
        }

        @Override
        public BCF2Type getDynamicType(final Object value) {
            return value == null ? BCF2Type.INT8 : BCF2Utils.determineIntegerType((int[])value);
        }

        @Requires("value == null || ((int[])value).length <= minValues")
        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            int count = 0;
            if ( value != null ) {
                for ( final int i : (int[])value ) {
                    encoder.encodeRawInt(i, type);
                    count++;
                }
            }
            for ( ; count < minValues; count++ ) encoder.encodeRawMissingValue(type);
        }
    }

    // ----------------------------------------------------------------------
    //
    // Subclass to encode List<Integer>
    //
    // ----------------------------------------------------------------------

    /**
     * Specialized int encoder for atomic (non-list) integers
     */
    public static class AtomicInt extends BCF2FieldEncoder {
        public AtomicInt(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, null);
        }

        @Override
        public BCF2Type getDynamicType(final Object value) {
            return value == null ? BCF2Type.INT8 : BCF2Utils.determineIntegerType((Integer)value);
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            int count = 0;
            if ( value != null ) {
                encoder.encodeRawInt((Integer)value, type);
                count++;
            }
            for ( ; count < minValues; count++ ) encoder.encodeRawMissingValue(type);
        }
    }

    public static class GenericInts extends BCF2FieldEncoder {
        public GenericInts(final VCFCompoundHeaderLine headerLine, final Map<String, Integer> dict ) {
            super(headerLine, dict, null);
        }

        @Override
        public BCF2Type getDynamicType(final Object value) {
            return value == null ? BCF2Type.INT8 : BCF2Utils.determineIntegerType(toList(Integer.class, value));
        }

        @Override
        public void encodeValue(final BCF2Encoder encoder, final Object value, final BCF2Type type, final int minValues) throws IOException {
            int count = 0;
            for ( final Integer i : toList(Integer.class, value) ) {
                if ( i != null ) { // necessary because .,. => [null, null] in VC
                    encoder.encodeRawInt(i, type);
                    count++;
                }
            }
            for ( ; count < minValues; count++ ) encoder.encodeRawMissingValue(type);
        }
    }


    // ----------------------------------------------------------------------
    //
    // Helper methods
    //
    // ----------------------------------------------------------------------

    /**
     * Helper function that takes an object and returns a list representation
     * of it:
     *
     * o == null => []
     * o is a list => o
     * else => [o]
     *
     * @param o
     * @return
     */
    private final static <T> List<T> toList(final Class<T> c, final Object o) {
        if ( o == null ) return Collections.emptyList();
        else if ( o instanceof List ) return (List<T>)o;
        else return Collections.singletonList((T)o);
    }
}
