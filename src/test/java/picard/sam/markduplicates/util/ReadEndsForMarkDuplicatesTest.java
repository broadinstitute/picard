/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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

package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Tests for {@link ReadEndsForMarkDuplicates} and {@link ReadEndsForMarkDuplicatesWithBarcodes}: the copy
 * constructors and {@code clone()}, which MarkDuplicates uses to start a pair from a fragment end, and
 * {@link ReadEndsForMarkDuplicates#isSingleEndReadKnownFragment(SAMRecord)}.
 */
public class ReadEndsForMarkDuplicatesTest {

    /** Every field set to a distinct value that differs from its default. */
    static ReadEndsForMarkDuplicates populate(final ReadEndsForMarkDuplicates ends) {
        ends.libraryId = 3;
        ends.orientation = ReadEnds.RF;
        ends.read1ReferenceIndex = 4;
        ends.read1Coordinate = 1005;
        ends.read2ReferenceIndex = 6;
        ends.read2Coordinate = 2007;
        ends.readGroup = 8;
        ends.tile = 1109;
        ends.x = 10010;
        ends.y = 20011;
        ends.orientationForOpticalDuplicates = ReadEnds.RR;
        ends.score = 1213;
        ends.read1IndexInFile = 3000000014L;
        ends.read2IndexInFile = 4000000015L;
        return ends;
    }

    static void assertSameFields(final ReadEndsForMarkDuplicates actual, final ReadEndsForMarkDuplicates expected) {
        Assert.assertEquals(actual.libraryId, expected.libraryId, "libraryId");
        Assert.assertEquals(actual.orientation, expected.orientation, "orientation");
        Assert.assertEquals(actual.read1ReferenceIndex, expected.read1ReferenceIndex, "read1ReferenceIndex");
        Assert.assertEquals(actual.read1Coordinate, expected.read1Coordinate, "read1Coordinate");
        Assert.assertEquals(actual.read2ReferenceIndex, expected.read2ReferenceIndex, "read2ReferenceIndex");
        Assert.assertEquals(actual.read2Coordinate, expected.read2Coordinate, "read2Coordinate");
        Assert.assertEquals(actual.readGroup, expected.readGroup, "readGroup");
        Assert.assertEquals(actual.tile, expected.tile, "tile");
        Assert.assertEquals(actual.x, expected.x, "x");
        Assert.assertEquals(actual.y, expected.y, "y");
        Assert.assertEquals(actual.orientationForOpticalDuplicates, expected.orientationForOpticalDuplicates, "orientationForOpticalDuplicates");
        Assert.assertEquals(actual.score, expected.score, "score");
        Assert.assertEquals(actual.read1IndexInFile, expected.read1IndexInFile, "read1IndexInFile");
        Assert.assertEquals(actual.read2IndexInFile, expected.read2IndexInFile, "read2IndexInFile");
    }

    @Test
    public void testCopyConstructorCopiesFields() {
        final ReadEndsForMarkDuplicates original = populate(new ReadEndsForMarkDuplicates());
        final ReadEndsForMarkDuplicates copy = new ReadEndsForMarkDuplicates(original);
        Assert.assertNotSame(copy, original);
        assertSameFields(copy, original);
        Assert.assertFalse(copy.isOpticalDuplicate);
    }

    /**
     * MarkDuplicates sends a fragment end to the fragment sort and stores its clone for the mate, then edits the
     * clone, so the two must not share state.
     */
    @Test
    public void testCloneIsIndependent() {
        final ReadEndsForMarkDuplicates original = populate(new ReadEndsForMarkDuplicates());
        final ReadEndsForMarkDuplicates clone = original.clone();
        Assert.assertNotSame(clone, original);
        assertSameFields(clone, original);

        clone.read2Coordinate = 1;
        clone.read2IndexInFile = 2;
        clone.orientation = ReadEnds.FF;
        clone.score = 3;
        clone.read1Coordinate = 4;
        Assert.assertEquals(original.read2Coordinate, 2007);
        Assert.assertEquals(original.read2IndexInFile, 4000000015L);
        Assert.assertEquals(original.orientation, ReadEnds.RF);
        Assert.assertEquals(original.score, 1213);
        Assert.assertEquals(original.read1Coordinate, 1005);
    }

    @Test
    public void testBarcodeDefaults() {
        final ReadEndsForMarkDuplicatesWithBarcodes ends = new ReadEndsForMarkDuplicatesWithBarcodes();
        Assert.assertEquals(ends.barcode, 0);
        Assert.assertEquals(ends.readOneBarcode, 0);
        Assert.assertEquals(ends.readTwoBarcode, 0);
        Assert.assertEquals(ends.read1Coordinate, -1);
    }

    @Test
    public void testBarcodeCopyFromPlainEndsLeavesBarcodesUnset() {
        final ReadEndsForMarkDuplicates original = populate(new ReadEndsForMarkDuplicates());
        final ReadEndsForMarkDuplicatesWithBarcodes copy = new ReadEndsForMarkDuplicatesWithBarcodes(original);
        assertSameFields(copy, original);
        Assert.assertEquals(copy.barcode, 0);
        Assert.assertEquals(copy.readOneBarcode, 0);
        Assert.assertEquals(copy.readTwoBarcode, 0);
    }

    @Test
    public void testBarcodeCopyAndClone() {
        final ReadEndsForMarkDuplicatesWithBarcodes original = (ReadEndsForMarkDuplicatesWithBarcodes) populate(new ReadEndsForMarkDuplicatesWithBarcodes());
        original.barcode = Integer.MIN_VALUE;
        original.readOneBarcode = -7;
        original.readTwoBarcode = Integer.MAX_VALUE;

        final ReadEndsForMarkDuplicatesWithBarcodes copy = new ReadEndsForMarkDuplicatesWithBarcodes(original);
        assertSameFields(copy, original);
        Assert.assertEquals(copy.barcode, Integer.MIN_VALUE);
        Assert.assertEquals(copy.readOneBarcode, -7);
        Assert.assertEquals(copy.readTwoBarcode, Integer.MAX_VALUE);

        final ReadEndsForMarkDuplicatesWithBarcodes clone = original.clone();
        Assert.assertNotSame(clone, original);
        assertSameFields(clone, original);
        Assert.assertEquals(clone.barcode, Integer.MIN_VALUE);
        Assert.assertEquals(clone.readOneBarcode, -7);
        Assert.assertEquals(clone.readTwoBarcode, Integer.MAX_VALUE);

        clone.readOneBarcode = 5;
        Assert.assertEquals(original.readOneBarcode, -7);

        // Called through the base type, as MarkDuplicates does; its barcode comparator reads the clone's barcodes.
        final ReadEndsForMarkDuplicatesWithBarcodes baseClone =
                (ReadEndsForMarkDuplicatesWithBarcodes) ((ReadEndsForMarkDuplicates) original).clone();
        assertSameFields(baseClone, original);
        Assert.assertEquals(baseClone.barcode, Integer.MIN_VALUE);
        Assert.assertEquals(baseClone.readOneBarcode, -7);
        Assert.assertEquals(baseClone.readTwoBarcode, Integer.MAX_VALUE);
    }

    private static final String CLIPPING_TAG = "tm";

    private static final SAMFileHeader HEADER = new SAMFileHeader();

    static {
        HEADER.addSequence(new SAMSequenceRecord("chr1", 10000));
    }

    private static SAMRecord fragment(final String cigar, final boolean negativeStrand) {
        final SAMRecord rec = new SAMRecord(HEADER);
        rec.setReadName("read");
        rec.setReferenceIndex(0);
        rec.setAlignmentStart(100);
        rec.setCigarString(cigar);
        rec.setReadNegativeStrandFlag(negativeStrand);
        return rec;
    }

    @DataProvider(name = "knownFragments")
    public Object[][] knownFragments() {
        final SAMRecord unmapped = fragment("50M", false);
        unmapped.setReadUnmappedFlag(true);
        unmapped.setAttribute(CLIPPING_TAG, "A");

        final SAMRecord paired = fragment("45M5S", false);
        paired.setReadPairedFlag(true);
        paired.setAttribute(CLIPPING_TAG, "A");

        final SAMRecord adapterTag = fragment("50M", false);
        adapterTag.setAttribute(CLIPPING_TAG, "A");
        final SAMRecord adapterAndQualityTag = fragment("50M", true);
        adapterAndQualityTag.setAttribute(CLIPPING_TAG, "QA");
        final SAMRecord qualityTag = fragment("50M", false);
        qualityTag.setAttribute(CLIPPING_TAG, "Q");
        final SAMRecord zTag = fragment("50M", false);
        zTag.setAttribute(CLIPPING_TAG, "Z");

        return new Object[][]{
                {"unmapped", unmapped, false},
                {"paired", paired, false},
                // A clipped 3' end (the end of a forward read, the start of a reverse one) marks a read that ran past its fragment.
                {"forward, unclipped", fragment("50M", false), false},
                {"forward, 3' soft clip", fragment("45M5S", false), true},
                {"forward, 3' hard clip", fragment("45M5H", false), true},
                {"forward, 5' soft clip", fragment("5S45M", false), false},
                {"reverse, unclipped", fragment("50M", true), false},
                {"reverse, 3' soft clip", fragment("5S45M", true), true},
                {"reverse, 3' hard clip", fragment("5H45M", true), true},
                {"reverse, 5' soft clip", fragment("45M5S", true), false},
                {"adapter clipped", adapterTag, true},
                {"adapter and quality clipped", adapterAndQualityTag, true},
                {"quality clipped", qualityTag, false},
                {"tm without A", zTag, false}
        };
    }

    @Test(dataProvider = "knownFragments")
    public void testIsSingleEndReadKnownFragment(final String description, final SAMRecord rec, final boolean expected) {
        Assert.assertEquals(ReadEndsForMarkDuplicates.isSingleEndReadKnownFragment(rec), expected, description);
    }
}
