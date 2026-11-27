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

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SortingCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import static picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesCodecTest.ORIENTATIONS;
import static picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesCodecTest.Values;
import static picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesCodecTest.assertDecoded;
import static picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesCodecTest.encode;
import static picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesCodecTest.fill;
import static picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesCodecTest.roundTrip;

/**
 * Tests for {@link ReadEndsForMarkDuplicatesWithBarcodesCodec}: round trips of the base record and the three
 * barcodes at extreme values for every orientation, decoding to the barcoded type, the end of a stream, truncated
 * records, rejection of records without barcodes, {@code clone()}, and spilling through a {@link SortingCollection}.
 */
public class ReadEndsForMarkDuplicatesWithBarcodesCodecTest {

    private static ReadEndsForMarkDuplicatesWithBarcodes barcoded(final Values values, final byte orientation) {
        final ReadEndsForMarkDuplicatesWithBarcodes ends = (ReadEndsForMarkDuplicatesWithBarcodes) fill(new ReadEndsForMarkDuplicatesWithBarcodes(), values, orientation);
        switch (values) {
            case MINIMUM:
                ends.barcode = ends.readOneBarcode = ends.readTwoBarcode = Integer.MIN_VALUE;
                break;
            case MAXIMUM:
                ends.barcode = ends.readOneBarcode = ends.readTwoBarcode = Integer.MAX_VALUE;
                break;
            case MINUS_ONE:
                ends.barcode = ends.readOneBarcode = ends.readTwoBarcode = -1;
                break;
            case ZERO:
                ends.barcode = ends.readOneBarcode = ends.readTwoBarcode = 0;
                break;
            case MIXED:
                ends.barcode = 0x41424344;
                ends.readOneBarcode = 0x45464748;
                ends.readTwoBarcode = 0x494A4B4C;
                break;
        }
        return ends;
    }

    private static void assertBarcodesDecoded(final ReadEndsForMarkDuplicates actual, final ReadEndsForMarkDuplicatesWithBarcodes expected, final boolean read2Kept) {
        assertDecoded(actual, expected, read2Kept);
        final ReadEndsForMarkDuplicatesWithBarcodes withBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) actual;
        Assert.assertEquals(withBarcodes.barcode, expected.barcode, "barcode");
        Assert.assertEquals(withBarcodes.readOneBarcode, expected.readOneBarcode, "readOneBarcode");
        Assert.assertEquals(withBarcodes.readTwoBarcode, expected.readTwoBarcode, "readTwoBarcode");
    }

    @DataProvider(name = "roundTrips")
    public Object[][] roundTrips() {
        final List<Object[]> rows = new ArrayList<>();
        for (final byte orientation : ORIENTATIONS) {
            for (final Values values : Values.values()) {
                rows.add(new Object[]{orientation, values});
            }
        }
        return rows.toArray(new Object[0][]);
    }

    @Test(dataProvider = "roundTrips")
    public void testRoundTrip(final byte orientation, final Values values) {
        final ReadEndsForMarkDuplicatesWithBarcodes ends = barcoded(values, orientation);
        final ReadEndsForMarkDuplicatesWithBarcodesCodec codec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        assertBarcodesDecoded(roundTrip(codec, (ReadEndsForMarkDuplicatesCodec) codec.clone(), ends), ends, orientation >= ReadEnds.FF);
    }

    /** Forced serialization set on the codec instance itself writes and reads the mate fields of fragments. */
    @Test(dataProvider = "roundTrips")
    public void testForcedRoundTripOnOneInstance(final byte orientation, final Values values) {
        final ReadEndsForMarkDuplicatesWithBarcodes ends = barcoded(values, orientation);
        final ReadEndsForMarkDuplicatesWithBarcodesCodec codec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        codec.setForceSerialization();
        assertBarcodesDecoded(roundTrip(codec, codec, ends), ends, true);
    }

    @Test
    public void testManyRecordsInOneStream() {
        final List<ReadEndsForMarkDuplicatesWithBarcodes> records = new ArrayList<>();
        for (int i = 0; i < ORIENTATIONS.length; i++) {
            final ReadEndsForMarkDuplicatesWithBarcodes ends = barcoded(Values.MIXED, ORIENTATIONS[(i * 5) % ORIENTATIONS.length]);
            ends.barcode = i;
            ends.readOneBarcode = 10 * i;
            ends.readTwoBarcode = -i;
            records.add(ends);
        }
        final ReadEndsForMarkDuplicatesWithBarcodesCodec codec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        final byte[] bytes = encode(codec, records.toArray(new ReadEndsForMarkDuplicates[0]));
        codec.setInputStream(new ByteArrayInputStream(bytes));
        for (final ReadEndsForMarkDuplicatesWithBarcodes expected : records) {
            assertBarcodesDecoded(codec.decode(), expected, expected.orientation >= ReadEnds.FF);
        }
        Assert.assertNull(codec.decode());
    }

    @Test
    public void testDecodeEmptyStreamReturnsNull() {
        final ReadEndsForMarkDuplicatesWithBarcodesCodec codec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        codec.setInputStream(new ByteArrayInputStream(new byte[0]));
        Assert.assertNull(codec.decode());
    }

    /** Whether to keep half of the encoded record or all of it but its last byte. */
    @DataProvider(name = "truncations")
    public Object[][] truncations() {
        return new Object[][]{{true}, {false}};
    }

    @Test(dataProvider = "truncations")
    public void testDecodeTruncatedRecordThrows(final boolean keepHalf) {
        final ReadEndsForMarkDuplicatesWithBarcodesCodec codec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        final byte[] bytes = encode(codec, barcoded(Values.MIXED, ReadEnds.R));
        codec.setInputStream(new ByteArrayInputStream(Arrays.copyOf(bytes, keepHalf ? bytes.length / 2 : bytes.length - 1)));
        final PicardException e = Assert.expectThrows(PicardException.class, codec::decode);
        Assert.assertTrue(e.getCause() instanceof EOFException, String.valueOf(e.getCause()));
    }

    @Test
    public void testEncodeRejectsEndsWithoutBarcodes() {
        final ReadEndsForMarkDuplicatesWithBarcodesCodec codec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        codec.setOutputStream(new ByteArrayOutputStream());
        final PicardException e = Assert.expectThrows(PicardException.class,
                () -> codec.encode(fill(new ReadEndsForMarkDuplicates(), Values.MIXED, ReadEnds.FR)));
        Assert.assertTrue(e.getMessage().contains("ReadEndsForMarkDuplicatesWithBarcodes"), e.getMessage());
    }

    @Test
    public void testEncodeFailureIsWrapped() {
        final ReadEndsForMarkDuplicatesWithBarcodesCodec codec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        codec.setOutputStream(new ReadEndsForMarkDuplicatesCodecTest.FailingOutputStream());
        final PicardException e = Assert.expectThrows(PicardException.class, () -> codec.encode(barcoded(Values.MIXED, ReadEnds.F)));
        Assert.assertTrue(e.getMessage().contains("writing"), e.getMessage());
        Assert.assertTrue(e.getCause() instanceof IOException);
    }

    @Test
    public void testClone() {
        final ReadEndsForMarkDuplicatesWithBarcodesCodec codec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        final SortingCollection.Codec<ReadEndsForMarkDuplicates> clone = codec.clone();
        Assert.assertNotSame(clone, codec);

        final ReadEndsForMarkDuplicatesWithBarcodes ends = barcoded(Values.MIXED, ReadEnds.RF);
        assertBarcodesDecoded(roundTrip((ReadEndsForMarkDuplicatesCodec) clone, codec, ends), ends, true);
    }

    @Test
    public void testSortingCollectionSpillRoundTrip() {
        final File tmpDir = IOUtil.createTempDir("barcodeCodecSpill").toFile();
        try {
            final SortingCollection<ReadEndsForMarkDuplicates> sorter = SortingCollection.newInstance(ReadEndsForMarkDuplicates.class,
                    new ReadEndsForMarkDuplicatesWithBarcodesCodec(), Comparator.comparingLong(e -> e.read1IndexInFile), 2, tmpDir);
            final List<ReadEndsForMarkDuplicatesWithBarcodes> expected = new ArrayList<>();
            for (int i = 0; i < 9; i++) {
                final ReadEndsForMarkDuplicatesWithBarcodes ends = barcoded(Values.MIXED, ORIENTATIONS[i % ORIENTATIONS.length]);
                ends.read1IndexInFile = i;
                ends.barcode = 100 - i;
                ends.readOneBarcode = i;
                ends.readTwoBarcode = i * i;
                expected.add(ends);
            }
            for (int i = expected.size() - 1; i >= 0; i--) {
                sorter.add(expected.get(i));
            }
            sorter.doneAdding();

            int i = 0;
            try (final CloseableIterator<ReadEndsForMarkDuplicates> it = sorter.iterator()) {
                while (it.hasNext()) {
                    final ReadEndsForMarkDuplicatesWithBarcodes exp = expected.get(i++);
                    assertBarcodesDecoded(it.next(), exp, exp.orientation >= ReadEnds.FF);
                }
            }
            Assert.assertEquals(i, expected.size());
            sorter.cleanup();
        } finally {
            IOUtil.deleteDirectoryTree(tmpDir);
        }
    }
}
