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

import org.testng.Assert;

/**
 * Misc methods for SAM-related unit tests.  These are in the src tree rather than the tests tree
 * so that they will be included in sam.jar, and therefore can be used by tests outside of net.sf.samtools.
 * These methods use org.testng.Assert methods.
 */
public class SAMTestUtil {
    /**
     * Basic sanity check for a pair of SAMRecords.
     */
    public void assertPairValid(final SAMRecord firstEnd, final SAMRecord secondEnd) {
        Assert.assertEquals(firstEnd.getReadName(), secondEnd.getReadName());
        Assert.assertTrue(firstEnd.getFirstOfPairFlag());
        Assert.assertTrue(secondEnd.getSecondOfPairFlag());
        Assert.assertFalse(secondEnd.getFirstOfPairFlag());
        Assert.assertFalse(firstEnd.getSecondOfPairFlag());
        if (!firstEnd.getReadUnmappedFlag() && !secondEnd.getReadUnmappedFlag()) {
            Assert.assertNotSame(firstEnd.getReadNegativeStrandFlag(),
                    secondEnd.getReadNegativeStrandFlag());
        }
    }

    /**
     * Basic sanity check for a SAMRecord.
     */
    public void assertReadValid(final SAMRecord read) {
        Assert.assertEquals(read.getReadBases().length, read.getBaseQualities().length);
        // Note that it is possible to have an unmapped read that has a coordinate
        if (read.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
            Assert.assertEquals(read.getAlignmentStart(), SAMRecord.NO_ALIGNMENT_START);
            Assert.assertTrue(read.getReadUnmappedFlag());
        } else {
            Assert.assertNotSame(read.getAlignmentStart(), SAMRecord.NO_ALIGNMENT_START);
        }
        if (read.getReadUnmappedFlag()) {
            Assert.assertEquals(read.getMappingQuality(), SAMRecord.NO_MAPPING_QUALITY);
            Assert.assertEquals(read.getCigar().getCigarElements().size(), 0);
        } else {
            Assert.assertNotSame(read.getCigar().getCigarElements(), 0);
        }
        if (read.getReadPairedFlag()) {
            if (read.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
                Assert.assertEquals(read.getMateAlignmentStart(), SAMRecord.NO_ALIGNMENT_START);
                Assert.assertTrue(read.getMateUnmappedFlag());
            } else {
                // Even if the mate is unmapped, if it has a reference name, it should have a position.
                Assert.assertNotSame(read.getMateAlignmentStart(), SAMRecord.NO_ALIGNMENT_START);
            }
            if (read.getReadUnmappedFlag() || read.getMateUnmappedFlag() ||
                    !read.getReferenceName().equals(read.getMateReferenceName())) {
                Assert.assertEquals(read.getInferredInsertSize(), 0);
            } else {
                Assert.assertNotSame(read.getInferredInsertSize(), 0);
            }
            if (!read.getReadUnmappedFlag() && !read.getMateUnmappedFlag()) {
                Assert.assertNotSame(read.getReadNegativeStrandFlag(), read.getMateNegativeStrandFlag(),
                        read.getReadName());
            }

        } else {
            Assert.assertEquals(read.getInferredInsertSize(), 0);
        }
    }
}
