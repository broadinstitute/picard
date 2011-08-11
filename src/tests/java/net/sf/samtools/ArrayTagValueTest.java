/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

import net.sf.picard.io.IoUtil;
import net.sf.picard.sam.SamFormatConverter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;

/**
 * Test 'B' array tag type.
 */
public class ArrayTagValueTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/samtools");
    private static final SortedSet<String> arrayTags = new TreeSet<String>(Arrays.asList(
            "X0:B:c,0",
            "X1:B:c,-128,127,-55,3",
            "X2:B:C,12",
            "X3:B:C,0,255,129,128,127",
            "X4:B:s,-23",
            "X5:B:s,32767,-32768,16384,-16384,0",
            "X6:B:S,65535",
            "X7:B:S,0,65535,13,14,15,22",
            "X8:B:i,0",
            "X9:B:i,2147483647,-2147483648,49,-1345678",
            "Y0:B:I,4000000000",
            "Y1:B:I,0,4294967295,3000000000"
            ));

    // Floats must be done separately because of ambiguity in string representation.
    private static final SortedMap<String, String> floatArrayTags = new TreeMap<String, String>();
    static {
        floatArrayTags.put("F0", "B:f,0");
        floatArrayTags.put("F1", "B:f,0,-0,.0,-.0,-0.0,3.4028234663852886E38f,-3.4028234663852886E38f,1.401298464324817E-45f,-1.401298464324817E-45f");

    }

    /**
     * Convert hand-crafted SAM text to BAM and then back to SAM, in order to test reading and writing of both SAM and BAM.
     */
    @Test
    public void basicTest() throws Exception {

        // Write a SAM file not using SAM-JDK
        final String inputSamWithNoTags = IoUtil.readFully(new FileInputStream(new File(TEST_DATA_DIR, "array_tags.sam")));
        final File inputSam = File.createTempFile("arrayTagValueTest.", ".sam");
        inputSam.deleteOnExit();
        final FileWriter os = new FileWriter(inputSam);
        os.write(inputSamWithNoTags);
        for (final String tag: arrayTags) {
            os.write("\t");
            os.write(tag);
        }
        for (final Map.Entry<String, String> floatTag: floatArrayTags.entrySet()) {
            os.write("\t");
            os.write(floatTag.getKey() + ":" + floatTag.getValue());
        }
        os.close();

        // Convert the SAM to BAM, in order to test SAM reading and BAM writing.
        final File bam = File.createTempFile("arrayTagValueTest.", ".bam");
        bam.deleteOnExit();

        final String[] samToBamArgs = {"INPUT="+inputSam.getAbsolutePath(),
                                 "OUTPUT="+bam.getAbsolutePath()};
        Assert.assertEquals(new SamFormatConverter().instanceMain(samToBamArgs), 0);

        // Convert the BAM to SAM, in order to test BAM reading and SAM writing.
        final File outputSam = File.createTempFile("arrayTagValueTest.convertedFromBam.", ".sam");
        outputSam.deleteOnExit();
        final String[] bamToSamArgs = {"INPUT="+bam.getAbsolutePath(),
                                 "OUTPUT="+outputSam.getAbsolutePath()};
        Assert.assertEquals(new SamFormatConverter().instanceMain(bamToSamArgs), 0);

        // Read the SAM file not using SAM-JDK, in order to validate that the tags are as expected.
        final BufferedReader reader = new BufferedReader(new FileReader(outputSam));
        // Skip header lines.
        reader.readLine(); reader.readLine();
        final String samTextLine = reader.readLine();

        // Get the single SAM record, and pull off the tags.
        final String[] fields = samTextLine.split("\t");
        final SortedSet<String> actualTags = new TreeSet<String>();
        final SortedMap<String, String> actualFloatTags = new TreeMap<String, String>();
        for (int i = 11; i < fields.length; ++i) {
            if (fields[i].startsWith("F")) {
                final String[] tagAndValue = fields[i].split(":", 2);
                actualFloatTags.put(tagAndValue[0], tagAndValue[1]);
            } else {
                actualTags.add(fields[i]);
            }
        }
        Assert.assertEquals(actualTags, arrayTags);

        // Floats must be done separately because of ambiguity in string representation.
        Assert.assertEquals(actualFloatTags.size(), floatArrayTags.size());

        for (final String tag: floatArrayTags.keySet()) {
            assertFloatTagValuesEqual(actualFloatTags.get(tag), floatArrayTags.get(tag));
        }
    }

    private void assertFloatTagValuesEqual(final String actualFloatStr, final String expectedFloatStr) {
        final float[] actualFloats = convertFloatStr(actualFloatStr);
        final float[] expectedFloats = convertFloatStr(expectedFloatStr);
        Assert.assertTrue(Arrays.equals(actualFloats, expectedFloats));
    }

    private float[] convertFloatStr(final String floatStr) {
        final String[] fields = floatStr.split(":")[1].split(",");
        final float[] ret = new float[fields.length - 1];
        for (int i = 1; i < fields.length; ++i) {
            ret[i-1] = Float.parseFloat(fields[i]);
        }
        return ret;
    }

    private static final String HEX_TAG = "X0";
    private static final String HEX_TAG_STRING = ":H:00ff7FF7";
    private static final byte[] EXPECTED_HEX_ARRAY = {0, (byte)0xff, 0x7f, (byte)0xf7};
    /**
     * Old-style hex arrays (type H) are converted to new-style byte arrays when read.
     * @throws Exception
     */
    @Test
    public void hexArrayFromSamTest() throws Exception {
        // Write a SAM file not using SAM-JDK
        final String inputSamWithNoTags = IoUtil.readFully(new FileInputStream(new File(TEST_DATA_DIR, "array_tags.sam")));
        final File inputSam = File.createTempFile("arrayTagValueTest.", ".sam");
        inputSam.deleteOnExit();
        final FileWriter os = new FileWriter(inputSam);
        os.write(inputSamWithNoTags);
        os.write("\t" + HEX_TAG + HEX_TAG_STRING);
        os.close();
        final SAMFileReader samReader = new SAMFileReader(inputSam);
        final SAMRecord samRecord = samReader.iterator().next();
        Assert.assertTrue(Arrays.equals((byte[])samRecord.getAttribute(HEX_TAG), EXPECTED_HEX_ARRAY));
        Assert.assertTrue(Arrays.equals(samRecord.getByteArrayAttribute(HEX_TAG), EXPECTED_HEX_ARRAY));
        Assert.assertTrue(Arrays.equals(samRecord.getSignedByteArrayAttribute(HEX_TAG), EXPECTED_HEX_ARRAY));
    }

    /**
     * Old-style hex arrays (type H) are converted to new-style byte arrays when read.  Confirm that a BAM containing
     * an old-style hex array is converted properly.
     * @throws Exception
     */
    @Test
    public void hexArrayFromBamTest() throws Exception {
        final SAMFileReader samReader = new SAMFileReader(new File(TEST_DATA_DIR, "recordWithHexArrayAttribute.bam"));
        final SAMRecord samRecord = samReader.iterator().next();
        Assert.assertTrue(Arrays.equals((byte[])samRecord.getAttribute(HEX_TAG), EXPECTED_HEX_ARRAY));
        Assert.assertTrue(Arrays.equals(samRecord.getByteArrayAttribute(HEX_TAG), EXPECTED_HEX_ARRAY));
        Assert.assertTrue(Arrays.equals(samRecord.getSignedByteArrayAttribute(HEX_TAG), EXPECTED_HEX_ARRAY));
    }
}
