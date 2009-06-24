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

import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;

/**
 * Confirm that integer tag types are stored and retrieved properly.
 *
 * @author alecw@broadinstitute.org
 */
public class SAMIntegerTagTest {
    private static final String BYTE_TAG = "BY";
    private static final String SHORT_TAG = "SH";
    private static final String INTEGER_TAG = "IN";
    private static final String UNSIGNED_INTEGER_TAG = "UI";
    private static final String STRING_TAG = "ST";

    @Test
    public void testBAM() throws Exception {
        final SAMRecord rec = writeAndReadSamRecord("bam");
        Assert.assertTrue(rec.getAttribute(BYTE_TAG) instanceof Byte);
        Assert.assertEquals(((Number)rec.getAttribute(BYTE_TAG)).intValue(), 1);
        Assert.assertTrue(rec.getAttribute(SHORT_TAG) instanceof Short);
        Assert.assertEquals(((Number)rec.getAttribute(SHORT_TAG)).intValue(), 1);
        Assert.assertTrue(rec.getAttribute(INTEGER_TAG) instanceof Integer);
        Assert.assertEquals(((Number)rec.getAttribute(INTEGER_TAG)).intValue(), 1);
    }

    @Test
    public void testSAM() throws Exception {
        final SAMRecord rec = writeAndReadSamRecord("sam");
        Assert.assertTrue(rec.getAttribute(BYTE_TAG) instanceof Integer);
        Assert.assertEquals(((Number)rec.getAttribute(BYTE_TAG)).intValue(), 1);
        Assert.assertTrue(rec.getAttribute(SHORT_TAG) instanceof Integer);
        Assert.assertEquals(((Number)rec.getAttribute(SHORT_TAG)).intValue(), 1);
        Assert.assertTrue(rec.getAttribute(INTEGER_TAG) instanceof Integer);
        Assert.assertEquals(((Number)rec.getAttribute(INTEGER_TAG)).intValue(), 1);
    }

    @Test
    public void testUnsignedIntegerBAM() throws Exception {
        SAMRecord rec = createSamRecord();
        final long val = 1l + Integer.MAX_VALUE;
        rec.setAttribute(UNSIGNED_INTEGER_TAG, val);
        rec = writeAndReadSamRecord("bam", rec);
        Assert.assertTrue(rec.getAttribute(UNSIGNED_INTEGER_TAG) instanceof Long);
        Assert.assertEquals(rec.getAttribute(UNSIGNED_INTEGER_TAG), val);
    }

    /**
     * Cannot store unsigned int in SAM text format.
     */
    @Test(expectedExceptions = SAMFormatException.class)
    public void testUnsignedIntegerSAM() throws Exception {
        final SAMRecord rec = createSamRecord();
        final long val = 1l + Integer.MAX_VALUE;
        rec.setAttribute(UNSIGNED_INTEGER_TAG, val);
        writeAndReadSamRecord("sam", rec);
        Assert.fail("Exception should have been thrown.");
    }

    @Test
    public void testGetTypedAttributeMethods() throws Exception {
        final SAMRecord rec = writeAndReadSamRecord("bam");
        Assert.assertEquals(rec.getByteAttribute(INTEGER_TAG).intValue(), 1);
        Assert.assertEquals(rec.getShortAttribute(INTEGER_TAG).intValue(), 1);
        Assert.assertEquals(rec.getIntegerAttribute(INTEGER_TAG).intValue(), 1);
    }

    /**
     * Should be an exception if a typed attribute call is made for the wrong type.
     */
    @Test(expectedExceptions = RuntimeException.class)
    public void testGetTypedAttributeForWrongType() throws Exception {
        final SAMRecord rec = createSamRecord();
        rec.setAttribute(STRING_TAG, "Hello, World!");
        writeAndReadSamRecord("bam", rec);
        rec.getIntegerAttribute(STRING_TAG);
        Assert.fail("Exception should have been thrown.");
    }

    /**
     * Should be an exception if a typed attribute call is made for a value that cannot
     * be coerced into the correct type.
     * This test is a little lame because a RuntimeException could be thrown for some other reason.
     */
    @Test(expectedExceptions = RuntimeException.class)
    public void testGetTypedAttributeOverflow() throws Exception {
        final SAMRecord rec = createSamRecord();
        rec.setAttribute(INTEGER_TAG, Integer.MAX_VALUE);
        writeAndReadSamRecord("bam", rec);
        rec.getShortAttribute(INTEGER_TAG);
        Assert.fail("Exception should have been thrown.");
    }

    /**
     * Should be an exception if a typed attribute call is made for a value that cannot
     * be coerced into the correct type.
     * This test is a little lame because a RuntimeException could be thrown for some other reason.
     */
    @Test(expectedExceptions = RuntimeException.class)
    public void testGetTypedAttributeUnerflow() throws Exception {
        final SAMRecord rec = createSamRecord();
        rec.setAttribute(INTEGER_TAG, Integer.MIN_VALUE);
        writeAndReadSamRecord("bam", rec);
        rec.getShortAttribute(INTEGER_TAG);
        Assert.fail("Exception should have been thrown.");
    }

    /**
     * Create a SAMRecord with integer tags of various sizes, write to a file, and read it back.
     * @param format "sam" or "bam".
     * @return The record after having being read from file.
     */
    private SAMRecord writeAndReadSamRecord(final String format) throws IOException {
        SAMRecord rec = createSamRecord();
        rec.setAttribute(BYTE_TAG, (byte)1);
        rec.setAttribute(SHORT_TAG, (short)1);
        rec.setAttribute(INTEGER_TAG, 1);
        rec = writeAndReadSamRecord(format, rec);
        return rec;
    }

    /**
     * Write a SAMRecord to a SAM file in the given format, and read it back.
     * @param format "sam" or "bam".
     * @param rec The record to write.
     * @return The same record, after having being written and read back.
     */
    private SAMRecord writeAndReadSamRecord(final String format, SAMRecord rec) throws IOException {
        final File bamFile = File.createTempFile("test.", "." + format);
        bamFile.deleteOnExit();
        final SAMFileWriter bamWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(rec.getHeader(), false, bamFile);
        bamWriter.addAlignment(rec);
        bamWriter.close();
        final SAMFileReader reader = new SAMFileReader(bamFile);
        rec = reader.iterator().next();
        reader.close();
        return rec;
    }

    private SAMRecord createSamRecord() {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.unsorted);
        builder.addFrag("readA", 20, 140, false);
        return builder.iterator().next();
    }
}
