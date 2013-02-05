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

// the imports for unit testing.
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.variantcontext.writer.BCF2Encoder;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


public class BCF2EncoderDecoderUnitTest extends VariantBaseTest {
    private final double FLOAT_TOLERANCE = 1e-6;
    final List<BCF2TypedValue> primitives = new ArrayList<BCF2TypedValue>();
    final List<BCF2TypedValue> basicTypes = new ArrayList<BCF2TypedValue>();
    final List<BCF2TypedValue> forCombinations = new ArrayList<BCF2TypedValue>();

    @BeforeSuite
    public void before() {
        basicTypes.add(new BCF2TypedValue(1, BCF2Type.INT8));
        basicTypes.add(new BCF2TypedValue(1000, BCF2Type.INT16));
        basicTypes.add(new BCF2TypedValue(1000000, BCF2Type.INT32));
        basicTypes.add(new BCF2TypedValue(1.2345e6, BCF2Type.FLOAT));
        basicTypes.add(new BCF2TypedValue("A", BCF2Type.CHAR));

        // small ints
        primitives.add(new BCF2TypedValue(0, BCF2Type.INT8));
        primitives.add(new BCF2TypedValue(10, BCF2Type.INT8));
        primitives.add(new BCF2TypedValue(-1, BCF2Type.INT8));
        primitives.add(new BCF2TypedValue(100, BCF2Type.INT8));
        primitives.add(new BCF2TypedValue(-100, BCF2Type.INT8));
        primitives.add(new BCF2TypedValue(-127, BCF2Type.INT8));    // last value in range
        primitives.add(new BCF2TypedValue( 127, BCF2Type.INT8));    // last value in range

        // medium ints
        primitives.add(new BCF2TypedValue(-1000, BCF2Type.INT16));
        primitives.add(new BCF2TypedValue(1000, BCF2Type.INT16));
        primitives.add(new BCF2TypedValue(-128, BCF2Type.INT16));    // first value in range
        primitives.add(new BCF2TypedValue( 128, BCF2Type.INT16));    // first value in range
        primitives.add(new BCF2TypedValue(-32767, BCF2Type.INT16)); // last value in range
        primitives.add(new BCF2TypedValue( 32767, BCF2Type.INT16)); // last value in range

        // larger ints
        primitives.add(new BCF2TypedValue(-32768, BCF2Type.INT32)); // first value in range
        primitives.add(new BCF2TypedValue( 32768, BCF2Type.INT32)); // first value in range
        primitives.add(new BCF2TypedValue(-100000, BCF2Type.INT32));
        primitives.add(new BCF2TypedValue(100000, BCF2Type.INT32));
        primitives.add(new BCF2TypedValue(-2147483647, BCF2Type.INT32));
        primitives.add(new BCF2TypedValue(2147483647, BCF2Type.INT32));

        // floats
        primitives.add(new BCF2TypedValue(0.0, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-0.0, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(1.0, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-1.0, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(1.1, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-1.1, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(5.0 / 3.0, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-5.0 / 3.0, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(1.23e3, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(1.23e6, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(1.23e9, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(1.23e12, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(1.23e15, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-1.23e3, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-1.23e6, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-1.23e9, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-1.23e12, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(-1.23e15, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(Float.MIN_VALUE, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(Float.MAX_VALUE, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(Double.NEGATIVE_INFINITY, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(Double.POSITIVE_INFINITY, BCF2Type.FLOAT));
        primitives.add(new BCF2TypedValue(Double.NaN, BCF2Type.FLOAT));

        // strings
        //primitives.add(new BCF2TypedValue("", BCFType.CHAR)); <- will be null (which is right)
        primitives.add(new BCF2TypedValue("S", BCF2Type.CHAR));
        primitives.add(new BCF2TypedValue("S2", BCF2Type.CHAR));
        primitives.add(new BCF2TypedValue("12345678910", BCF2Type.CHAR));
        primitives.add(new BCF2TypedValue("ABCDEFGHIJKLMNOPQRSTUVWXYZ", BCF2Type.CHAR));
        primitives.add(new BCF2TypedValue("ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ", BCF2Type.CHAR));

        // missing values
        for ( BCF2Type type : BCF2Type.values() ) {
            primitives.add(new BCF2TypedValue(null, type));
        }

        forCombinations.add(new BCF2TypedValue(10, BCF2Type.INT8));
        forCombinations.add(new BCF2TypedValue(100, BCF2Type.INT8));
        forCombinations.add(new BCF2TypedValue(-100, BCF2Type.INT8));
        forCombinations.add(new BCF2TypedValue(-128, BCF2Type.INT16));    // first value in range
        forCombinations.add(new BCF2TypedValue( 128, BCF2Type.INT16));    // first value in range
        forCombinations.add(new BCF2TypedValue(-100000, BCF2Type.INT32));
        forCombinations.add(new BCF2TypedValue(100000, BCF2Type.INT32));
        forCombinations.add(new BCF2TypedValue(0.0, BCF2Type.FLOAT));
        forCombinations.add(new BCF2TypedValue(1.23e6, BCF2Type.FLOAT));
        forCombinations.add(new BCF2TypedValue(-1.23e6, BCF2Type.FLOAT));
        forCombinations.add(new BCF2TypedValue("S", BCF2Type.CHAR));
        forCombinations.add(new BCF2TypedValue("ABCDEFGHIJKLMNOPQRSTUVWXYZ", BCF2Type.CHAR));
        forCombinations.add(new BCF2TypedValue("ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ", BCF2Type.CHAR));

        // missing values
        for ( BCF2Type type : BCF2Type.values() ) {
            forCombinations.add(new BCF2TypedValue(null, type));
        }
    }

    // --------------------------------------------------------------------------------
    //
    // merge case Provider
    //
    // --------------------------------------------------------------------------------

    private class BCF2TypedValue {
        final BCF2Type type;
        final Object value;

        private BCF2TypedValue(final int value, final BCF2Type type) {
            this(new Integer(value), type);
        }

        private BCF2TypedValue(final double value, final BCF2Type type) {
            this(new Double(value), type);
        }

        private BCF2TypedValue(final Object value, final BCF2Type type) {
            this.type = type;
            this.value = value;
        }

        public boolean isMissing() { return value == null; }

        @Override
        public String toString() {
            return String.format("%s of %s", value, type);
        }
    }

    // -----------------------------------------------------------------
    //
    // Test encoding of basic types
    //
    // -----------------------------------------------------------------

    @DataProvider(name = "BCF2EncodingTestProviderBasicTypes")
    public Object[][] BCF2EncodingTestProviderBasicTypes() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( BCF2TypedValue tv : basicTypes )
            tests.add(new Object[]{Arrays.asList(tv)});
        return tests.toArray(new Object[][]{});
    }

    private interface EncodeMe {
        public void encode(final BCF2Encoder encoder, final BCF2TypedValue tv) throws IOException;
    }


    @Test(dataProvider = "BCF2EncodingTestProviderBasicTypes")
    public void testBCF2BasicTypesWithStaticCalls(final List<BCF2TypedValue> toEncode) throws IOException {
        testBCF2BasicTypesWithEncodeMe(toEncode,
                new EncodeMe() {
                    @Override
                    public void encode(final BCF2Encoder encoder, final BCF2TypedValue tv) throws IOException {
                        switch ( tv.type ) {
                            case INT8:
                            case INT16:
                            case INT32:
                                encoder.encodeTypedInt((Integer)tv.value, tv.type);
                                break;
                            case FLOAT:
                                encoder.encodeTypedFloat((Double)tv.value);
                                break;
                            case CHAR:
                                encoder.encodeTypedString((String)tv.value);
                                break;
                        }
                    }
                });
    }

    @Test(dataProvider = "BCF2EncodingTestProviderBasicTypes")
    public void testBCF2BasicTypesWithObjectType(final List<BCF2TypedValue> toEncode) throws IOException {
        testBCF2BasicTypesWithEncodeMe(toEncode,
                new EncodeMe() {
                    @Override
                    public void encode(final BCF2Encoder encoder, final BCF2TypedValue tv) throws IOException {
                        encoder.encodeTyped(tv.value, tv.type);
                    }
                });
    }

    @Test(dataProvider = "BCF2EncodingTestProviderBasicTypes")
    public void testBCF2BasicTypesWithObjectNoType(final List<BCF2TypedValue> toEncode) throws IOException {
        testBCF2BasicTypesWithEncodeMe(toEncode,
                new EncodeMe() {
                    @Override
                    public void encode(final BCF2Encoder encoder, final BCF2TypedValue tv) throws IOException {
                        encoder.encode(tv.value);
                    }
                });
    }

    public void testBCF2BasicTypesWithEncodeMe(final List<BCF2TypedValue> toEncode, final EncodeMe func) throws IOException {
        for ( final BCF2TypedValue tv : toEncode ) {
            BCF2Encoder encoder = new BCF2Encoder();
            func.encode(encoder, tv);

            BCF2Decoder decoder = new BCF2Decoder(encoder.getRecordBytes());
            final Object decoded = decoder.decodeTypedValue();

            Assert.assertNotNull(decoded);
            Assert.assertFalse(decoded instanceof List);
            myAssertEquals(tv, decoded);
        }
    }

    @Test(dataProvider = "BCF2EncodingTestProviderBasicTypes")
    public void testBCF2EncodingVectors(final List<BCF2TypedValue> toEncode) throws IOException {
        for ( final BCF2TypedValue tv : toEncode ) {
            for ( final int length : Arrays.asList(2, 5, 10, 15, 20, 25) ) {
                BCF2Encoder encoder = new BCF2Encoder();
                List<Object> expected = Collections.nCopies(length, tv.value);
                encoder.encodeTyped(expected, tv.type);

                BCF2Decoder decoder = new BCF2Decoder(encoder.getRecordBytes());
                final Object decoded = decoder.decodeTypedValue();

                Assert.assertTrue(decoded instanceof List);
                final List<Object> decodedList = (List<Object>)decoded;
                Assert.assertEquals(decodedList.size(), expected.size());
                for ( Object decodedValue : decodedList )
                    myAssertEquals(tv, decodedValue);
            }
        }
    }

    @DataProvider(name = "BCF2EncodingTestProviderSingletons")
    public Object[][] BCF2EncodingTestProviderSingletons() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( BCF2TypedValue tv : primitives )
            tests.add(new Object[]{Arrays.asList(tv)});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BCF2EncodingTestProviderSingletons")
    public void testBCF2EncodingSingletons(final List<BCF2TypedValue> toEncode) throws IOException {
        final byte[] record = encodeRecord(toEncode);
        decodeRecord(toEncode, record);
    }

    // -----------------------------------------------------------------
    //
    // Test encoding of vectors
    //
    // -----------------------------------------------------------------

    @DataProvider(name = "BCF2EncodingTestProviderSequences")
    public Object[][] BCF2EncodingTestProviderSequences() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( BCF2TypedValue tv1 : forCombinations )
            for ( BCF2TypedValue tv2 : forCombinations )
                for ( BCF2TypedValue tv3 : forCombinations )
                    tests.add(new Object[]{Arrays.asList(tv1, tv2, tv3)});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BCF2EncodingTestProviderBasicTypes")
    public void testBCF2EncodingVectorsWithMissing(final List<BCF2TypedValue> toEncode) throws IOException {
        for ( final BCF2TypedValue tv : toEncode ) {
            if ( tv.type != BCF2Type.CHAR ) {
                for ( final int length : Arrays.asList(2, 5, 10, 15, 20, 25) ) {
                    final byte td = BCF2Utils.encodeTypeDescriptor(1, tv.type);

                    final BCF2Encoder encoder = new BCF2Encoder();
                    for ( int i = 0; i < length; i++ ) {
                        encoder.encodeRawValue(i % 2 == 0 ? null : tv.value, tv.type);
                    }

                    final BCF2Decoder decoder = new BCF2Decoder(encoder.getRecordBytes());

                    for ( int i = 0; i < length; i++ ) {
                        final Object decoded = decoder.decodeTypedValue(td);
                        myAssertEquals(i % 2 == 0 ? new BCF2TypedValue(null, tv.type) : tv, decoded);
                    }
                }
            }
        }
    }

    @Test(dataProvider = "BCF2EncodingTestProviderSequences", dependsOnMethods = "testBCF2EncodingSingletons")
    public void testBCF2EncodingTestProviderSequences(final List<BCF2TypedValue> toEncode) throws IOException {
        final byte[] record = encodeRecord(toEncode);
        decodeRecord(toEncode, record);
    }

    // -----------------------------------------------------------------
    //
    // Test strings and lists of strings
    //
    // -----------------------------------------------------------------

    @DataProvider(name = "ListOfStrings")
    public Object[][] listOfStringsProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();
        tests.add(new Object[]{Arrays.asList("s1", "s2"), ",s1,s2"});
        tests.add(new Object[]{Arrays.asList("s1", "s2", "s3"), ",s1,s2,s3"});
        tests.add(new Object[]{Arrays.asList("s1", "s2", "s3", "s4"), ",s1,s2,s3,s4"});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ListOfStrings")
    public void testEncodingListOfString(List<String> strings, String expected) throws IOException {
        final String collapsed = BCF2Utils.collapseStringList(strings);
        Assert.assertEquals(collapsed, expected);
        Assert.assertEquals(BCF2Utils.explodeStringList(collapsed), strings);
    }

    // -----------------------------------------------------------------
    //
    // Tests to determine the best type of arrays of integers
    //
    // -----------------------------------------------------------------

    @DataProvider(name = "BestIntTypeTests")
    public Object[][] BestIntTypeTests() {
        List<Object[]> tests = new ArrayList<Object[]>();
        tests.add(new Object[]{Arrays.asList(1), BCF2Type.INT8});
        tests.add(new Object[]{Arrays.asList(1, 10), BCF2Type.INT8});
        tests.add(new Object[]{Arrays.asList(1, 10, 100), BCF2Type.INT8});
        tests.add(new Object[]{Arrays.asList(1, -1), BCF2Type.INT8});
        tests.add(new Object[]{Arrays.asList(1, 1000), BCF2Type.INT16});
        tests.add(new Object[]{Arrays.asList(1, 1000, 10), BCF2Type.INT16});
        tests.add(new Object[]{Arrays.asList(1, 1000, 100), BCF2Type.INT16});
        tests.add(new Object[]{Arrays.asList(1000), BCF2Type.INT16});
        tests.add(new Object[]{Arrays.asList(100000), BCF2Type.INT32});
        tests.add(new Object[]{Arrays.asList(100000, 10), BCF2Type.INT32});
        tests.add(new Object[]{Arrays.asList(100000, 100), BCF2Type.INT32});
        tests.add(new Object[]{Arrays.asList(100000, 1, -10), BCF2Type.INT32});
        tests.add(new Object[]{Arrays.asList(-100000, 1, -10), BCF2Type.INT32});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BestIntTypeTests")
    public void determineBestEncoding(final List<Integer> ints, final BCF2Type expectedType) throws IOException {
        BCF2Encoder encoder = new BCF2Encoder();
        Assert.assertEquals(BCF2Utils.determineIntegerType(ints), expectedType);
        Assert.assertEquals(BCF2Utils.determineIntegerType(toPrimitive(ints.toArray(new Integer[0]))), expectedType);
    }

    private static int[] toPrimitive ( final Integer[] array ) {
        if ( array == null ) {
            return null;
        }
        else if ( array.length == 0 ) {
            return new int[0];
        }

        final int[] result = new int[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i].intValue();
        }
        return result;
    }

    // -----------------------------------------------------------------
    //
    // Tests managing and skipping multiple blocks
    //
    // -----------------------------------------------------------------

    @Test(dataProvider = "BCF2EncodingTestProviderSequences", dependsOnMethods = "testBCF2EncodingTestProviderSequences")
    public void testReadAndSkipWithMultipleBlocks(final List<BCF2TypedValue> block) throws IOException {
        testReadAndSkipWithMultipleBlocks(block, forCombinations);
        testReadAndSkipWithMultipleBlocks(forCombinations, block);
    }

    public void testReadAndSkipWithMultipleBlocks(final List<BCF2TypedValue> block1, final List<BCF2TypedValue> block2) throws IOException {
        final byte[] record1 = encodeRecord(block1);
        final byte[] record2 = encodeRecord(block2);

        // each record is individually good
        decodeRecord(block1, record1);
        decodeRecord(block2, record2);

        BCF2Decoder decoder = new BCF2Decoder();

        // test setting
        decoder.setRecordBytes(record1);
        decodeRecord(block1, decoder);
        decoder.setRecordBytes(record2);
        decodeRecord(block2, decoder);

        // test combining the streams
        final byte[] combined = combineRecords(record1, record2);
        final List<BCF2TypedValue> combinedObjects = new ArrayList<BCF2TypedValue>(block1);
        combinedObjects.addAll(block2);

        // the combined bytes is the same as the combined objects
        InputStream stream = new ByteArrayInputStream(combined);
        decoder.readNextBlock(record1.length, stream);
        decodeRecord(block1, decoder);
        decoder.readNextBlock(record2.length, stream);
        decodeRecord(block2, decoder);

        // skipping the first block allows us to read the second block directly
        stream = new ByteArrayInputStream(combined);
        decoder.skipNextBlock(record1.length, stream);
        decoder.readNextBlock(record2.length, stream);
        decodeRecord(block2, decoder);
    }

    // -----------------------------------------------------------------
    //
    // Test encoding / decoding arrays of ints
    //
    // This checks that we can encode and decode correctly with the
    // low-level decodeIntArray function arrays of values.  This
    // has to be pretty comprehensive as decodeIntArray is a highly optimized
    // piece of code with lots of edge cases.  The values we are encoding
    // don't really matter -- just that the values come back as expected.
    //
    // -----------------------------------------------------------------

    @DataProvider(name = "IntArrays")
    public Object[][] makeIntArrays() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( int nValues : Arrays.asList(0, 1, 2, 5, 10, 100) ) {
            for ( int nPad : Arrays.asList(0, 1, 2, 5, 10, 100) ) {
                int nElements = nValues + nPad;

                List<Integer> values = new ArrayList<Integer>(nElements);

                // add nValues from 0 to nValues - 1
                for ( int i = 0; i < nValues; i++ )
                    values.add(i);

                // add nPad nulls
                for ( int i = 0; i < nPad; i++ )
                    values.add(null);

                tests.add(new Object[]{values});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "IntArrays")
    public void testIntArrays(final List<Integer> ints) throws IOException {
        final BCF2Encoder encoder = new BCF2Encoder();
        encoder.encodeTyped(ints, BCF2Type.INT16);

        final BCF2Decoder decoder = new BCF2Decoder(encoder.getRecordBytes());

        final byte typeDescriptor = decoder.readTypeDescriptor();

        // read the int[] with the low-level version
        final int size = decoder.decodeNumberOfElements(typeDescriptor);
        final int[] decoded = decoder.decodeIntArray(typeDescriptor, size);

        if ( isMissing(ints) ) {
            // we expect that the result is null in this case
            Assert.assertNull(decoded, "Encoded all missing values -- expected null");
        } else {
            // we expect at least some values to come back
            Assert.assertTrue(decoded.length > 0, "Must have at least 1 element for non-null encoded data");

            // check corresponding values
            for ( int i = 0; i < ints.size(); i++ ) {
                final Integer expected = ints.get(i);

                if ( expected == null ) {
                    Assert.assertTrue(decoded.length <= i, "we expect decoded to be truncated for missing values");
                } else {
                    Assert.assertTrue(decoded.length > i, "we expected at least " + i + " values in decoded array");
                    Assert.assertEquals(decoded[i], (int)expected);
                }
            }
        }
    }

    // -----------------------------------------------------------------
    //
    // Helper routines
    //
    // -----------------------------------------------------------------

    private final byte[] combineRecords(final byte[] record1, final byte[] record2) throws IOException {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        baos.write(record1);
        baos.write(record2);
        return baos.toByteArray();
    }

    private final byte[] encodeRecord(final List<BCF2TypedValue> toEncode) throws IOException {
        BCF2Encoder encoder = new BCF2Encoder();

        for ( final BCF2TypedValue tv : toEncode ) {
            if ( tv.isMissing() )
                encoder.encodeTypedMissing(tv.type);
            else {
                final BCF2Type encodedType = encoder.encode(tv.value);
                if ( tv.type != null ) // only if we have an expectation
                    Assert.assertEquals(encodedType, tv.type);
            }
        }

        // check output
        final byte[] record = encoder.getRecordBytes();
        Assert.assertNotNull(record);
        Assert.assertTrue(record.length > 0);
        return record;
    }

    private final void decodeRecord(final List<BCF2TypedValue> toEncode, final byte[] record) throws IOException {
        decodeRecord(toEncode, new BCF2Decoder(record));
    }

    private final void decodeRecord(final List<BCF2TypedValue> toEncode, final BCF2Decoder decoder) throws IOException {
        for ( final BCF2TypedValue tv : toEncode ) {
            Assert.assertFalse(decoder.blockIsFullyDecoded());
            final Object decoded = decoder.decodeTypedValue();

            myAssertEquals(tv, decoded);
        }

        Assert.assertTrue(decoder.blockIsFullyDecoded());
    }

    private final void myAssertEquals(final BCF2TypedValue tv, final Object decoded) {
        if ( tv.value == null ) { // special needs for instanceof double
            Assert.assertEquals(decoded, tv.value);
        } else if ( tv.type == BCF2Type.FLOAT ) { // need tolerance for floats, and they aren't null
            Assert.assertTrue(decoded instanceof Double);

            final double valueFloat = (Double)tv.value;
            final double decodedFloat = (Double)decoded;

            VariantBaseTest.assertEqualsDoubleSmart(decodedFloat, valueFloat, FLOAT_TOLERANCE);
        } else
            Assert.assertEquals(decoded, tv.value);
    }

    private final boolean isMissing(final List<Integer> values) {
        if ( values != null )
            for ( Integer value : values )
                if ( value != null )
                    return false;
        return true;
    }
}