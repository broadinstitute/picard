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
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Tests for {@link ReadEndsForMarkDuplicatesCodec}: round trips of every field at its extreme values for every
 * orientation byte, where the mate's coordinate and file index of a fragment are kept only when serialization is
 * forced; several records in one stream; decoding at and before the end of a stream; write failures;
 * {@code clone()}; and spilling through a {@link SortingCollection}, which reads its spill files with clones of the
 * codec.
 */
public class ReadEndsForMarkDuplicatesCodecTest {

    static final byte[] ORIENTATIONS = {ReadEnds.F, ReadEnds.R, ReadEnds.FF, ReadEnds.FR, ReadEnds.RR, ReadEnds.RF};

    enum Values {MINIMUM, MAXIMUM, MINUS_ONE, ZERO, MIXED}

    static ReadEndsForMarkDuplicates fill(final ReadEndsForMarkDuplicates ends, final Values values, final byte orientation) {
        ends.orientation = orientation;
        switch (values) {
            case MINIMUM:
                ends.score = Short.MIN_VALUE;
                ends.libraryId = Short.MIN_VALUE;
                ends.read1ReferenceIndex = Integer.MIN_VALUE;
                ends.read1Coordinate = Integer.MIN_VALUE;
                ends.read1IndexInFile = Long.MIN_VALUE;
                ends.read2ReferenceIndex = Integer.MIN_VALUE;
                ends.read2Coordinate = Integer.MIN_VALUE;
                ends.read2IndexInFile = Long.MIN_VALUE;
                ends.readGroup = Short.MIN_VALUE;
                ends.tile = Short.MIN_VALUE;
                ends.x = Integer.MIN_VALUE;
                ends.y = Integer.MIN_VALUE;
                ends.orientationForOpticalDuplicates = Byte.MIN_VALUE;
                ends.duplicateSetSize = Integer.MIN_VALUE;
                break;
            case MAXIMUM:
                ends.score = Short.MAX_VALUE;
                ends.libraryId = Short.MAX_VALUE;
                ends.read1ReferenceIndex = Integer.MAX_VALUE;
                ends.read1Coordinate = Integer.MAX_VALUE;
                ends.read1IndexInFile = Long.MAX_VALUE;
                ends.read2ReferenceIndex = Integer.MAX_VALUE;
                ends.read2Coordinate = Integer.MAX_VALUE;
                ends.read2IndexInFile = Long.MAX_VALUE;
                ends.readGroup = Short.MAX_VALUE;
                ends.tile = Short.MAX_VALUE;
                ends.x = Integer.MAX_VALUE;
                ends.y = Integer.MAX_VALUE;
                ends.orientationForOpticalDuplicates = Byte.MAX_VALUE;
                ends.duplicateSetSize = Integer.MAX_VALUE;
                break;
            case MINUS_ONE:
                ends.score = -1;
                ends.libraryId = -1;
                ends.read1ReferenceIndex = -1;
                ends.read1Coordinate = -1;
                ends.read1IndexInFile = -1;
                ends.read2ReferenceIndex = -1;
                ends.read2Coordinate = -1;
                ends.read2IndexInFile = -1;
                ends.readGroup = -1;
                ends.tile = -1;
                ends.x = -1;
                ends.y = -1;
                ends.orientationForOpticalDuplicates = -1;
                ends.duplicateSetSize = -1;
                break;
            case ZERO:
                ends.score = 0;
                ends.libraryId = 0;
                ends.read1ReferenceIndex = 0;
                ends.read1Coordinate = 0;
                ends.read1IndexInFile = 0;
                ends.read2ReferenceIndex = 0;
                ends.read2Coordinate = 0;
                ends.read2IndexInFile = 0;
                ends.readGroup = 0;
                ends.tile = 0;
                ends.x = 0;
                ends.y = 0;
                ends.orientationForOpticalDuplicates = 0;
                ends.duplicateSetSize = 0;
                break;
            case MIXED:
                // Distinct values in every byte position catch swapped or misaligned fields.
                ends.score = 0x0102;
                ends.libraryId = 0x0304;
                ends.read1ReferenceIndex = 0x05060708;
                ends.read1Coordinate = 0x090A0B0C;
                ends.read1IndexInFile = 0x0D0E0F1011121314L;
                ends.read2ReferenceIndex = 0x15161718;
                ends.read2Coordinate = 0x191A1B1C;
                ends.read2IndexInFile = 0x1D1E1F2021222324L;
                ends.readGroup = 0x2526;
                ends.tile = 0x2728;
                ends.x = 0x292A2B2C;
                ends.y = 0x2D2E2F30;
                ends.orientationForOpticalDuplicates = 0x31;
                ends.duplicateSetSize = 0x32333435;
                break;
        }
        return ends;
    }

    static byte[] encode(final ReadEndsForMarkDuplicatesCodec codec, final ReadEndsForMarkDuplicates... records) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        codec.setOutputStream(bytes);
        for (final ReadEndsForMarkDuplicates record : records) {
            codec.encode(record);
        }
        return bytes.toByteArray();
    }

    /** Encodes one record with one codec and decodes it with another, as SortingCollection does with clones. */
    static ReadEndsForMarkDuplicates roundTrip(final ReadEndsForMarkDuplicatesCodec encoder, final ReadEndsForMarkDuplicatesCodec decoder,
                                               final ReadEndsForMarkDuplicates ends) {
        decoder.setInputStream(new ByteArrayInputStream(encode(encoder, ends)));
        final ReadEndsForMarkDuplicates decoded = decoder.decode();
        Assert.assertNull(decoder.decode());
        return decoded;
    }

    /** Checks every field the codec keeps; the read2 fields of a fragment are kept only when serialization is forced. */
    static void assertDecoded(final ReadEndsForMarkDuplicates actual, final ReadEndsForMarkDuplicates expected, final boolean read2Kept) {
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.score, expected.score, "score");
        Assert.assertEquals(actual.libraryId, expected.libraryId, "libraryId");
        Assert.assertEquals(actual.orientation, expected.orientation, "orientation");
        Assert.assertEquals(actual.read1ReferenceIndex, expected.read1ReferenceIndex, "read1ReferenceIndex");
        Assert.assertEquals(actual.read1Coordinate, expected.read1Coordinate, "read1Coordinate");
        Assert.assertEquals(actual.read1IndexInFile, expected.read1IndexInFile, "read1IndexInFile");
        Assert.assertEquals(actual.read2ReferenceIndex, expected.read2ReferenceIndex, "read2ReferenceIndex");
        if (read2Kept) {
            Assert.assertEquals(actual.read2Coordinate, expected.read2Coordinate, "read2Coordinate");
            Assert.assertEquals(actual.read2IndexInFile, expected.read2IndexInFile, "read2IndexInFile");
        }
        Assert.assertEquals(actual.readGroup, expected.readGroup, "readGroup");
        Assert.assertEquals(actual.tile, expected.tile, "tile");
        Assert.assertEquals(actual.x, expected.x, "x");
        Assert.assertEquals(actual.y, expected.y, "y");
        Assert.assertEquals(actual.orientationForOpticalDuplicates, expected.orientationForOpticalDuplicates, "orientationForOpticalDuplicates");
    }

    @DataProvider(name = "roundTrips")
    public Object[][] roundTrips() {
        final List<Object[]> rows = new ArrayList<>();
        for (final byte orientation : ORIENTATIONS) {
            for (final Values values : Values.values()) {
                for (final boolean force : new boolean[]{false, true}) {
                    rows.add(new Object[]{orientation, values, force});
                }
            }
        }
        return rows.toArray(new Object[0][]);
    }

    @Test(dataProvider = "roundTrips")
    public void testRoundTrip(final byte orientation, final Values values, final boolean force) {
        final ReadEndsForMarkDuplicates ends = fill(new ReadEndsForMarkDuplicates(), values, orientation);
        final ReadEndsForMarkDuplicatesCodec codec = new ReadEndsForMarkDuplicatesCodec();
        if (force) codec.setForceSerialization();

        // Decode with a clone, as SortingCollection does for its spill files.
        final ReadEndsForMarkDuplicates decoded = roundTrip(codec, (ReadEndsForMarkDuplicatesCodec) codec.clone(), ends);
        assertDecoded(decoded, ends, force || orientation >= ReadEnds.FF);
        Assert.assertEquals(decoded.duplicateSetSize, ends.duplicateSetSize);
    }

    @Test
    public void testManyRecordsInOneStream() {
        final List<ReadEndsForMarkDuplicates> records = new ArrayList<>();
        int i = 0;
        for (final byte orientation : new byte[]{ReadEnds.F, ReadEnds.FR, ReadEnds.R, ReadEnds.RF, ReadEnds.RR, ReadEnds.F, ReadEnds.FF}) {
            final ReadEndsForMarkDuplicates ends = fill(new ReadEndsForMarkDuplicates(), Values.MIXED, orientation);
            ends.read1Coordinate = 1000 + i;
            ends.read2Coordinate = 2000 + i;
            ends.read1IndexInFile = i;
            i++;
            records.add(ends);
        }

        final ReadEndsForMarkDuplicatesCodec codec = new ReadEndsForMarkDuplicatesCodec();
        final byte[] bytes = encode(codec, records.toArray(new ReadEndsForMarkDuplicates[0]));
        codec.setInputStream(new ByteArrayInputStream(bytes));
        for (final ReadEndsForMarkDuplicates expected : records) {
            assertDecoded(codec.decode(), expected, expected.orientation >= ReadEnds.FF);
        }
        Assert.assertNull(codec.decode());
        Assert.assertNull(codec.decode());
    }

    @Test
    public void testDecodeEmptyStreamReturnsNull() {
        final ReadEndsForMarkDuplicatesCodec codec = new ReadEndsForMarkDuplicatesCodec();
        codec.setInputStream(new ByteArrayInputStream(new byte[0]));
        Assert.assertNull(codec.decode());
    }

    /** An orientation, and whether to keep half of the encoded record or all of it but its last byte. */
    @DataProvider(name = "truncations")
    public Object[][] truncations() {
        return new Object[][]{
                {ReadEnds.F, true},
                {ReadEnds.F, false},
                {ReadEnds.FR, true},
                {ReadEnds.FR, false}
        };
    }

    /** A stream that ends inside a record is an error, not the end of the data. */
    @Test(dataProvider = "truncations")
    public void testDecodeTruncatedRecordThrows(final byte orientation, final boolean keepHalf) {
        final ReadEndsForMarkDuplicatesCodec codec = new ReadEndsForMarkDuplicatesCodec();
        final byte[] bytes = encode(codec, fill(new ReadEndsForMarkDuplicates(), Values.MIXED, orientation));
        codec.setInputStream(new ByteArrayInputStream(Arrays.copyOf(bytes, keepHalf ? bytes.length / 2 : bytes.length - 1)));
        final PicardException e = Assert.expectThrows(PicardException.class, codec::decode);
        Assert.assertTrue(e.getCause() instanceof EOFException, String.valueOf(e.getCause()));
    }

    @Test
    public void testEncodeFailureIsWrapped() {
        final ReadEndsForMarkDuplicatesCodec codec = new ReadEndsForMarkDuplicatesCodec();
        codec.setOutputStream(new FailingOutputStream());
        final PicardException e = Assert.expectThrows(PicardException.class,
                () -> codec.encode(fill(new ReadEndsForMarkDuplicates(), Values.MIXED, ReadEnds.FR)));
        Assert.assertTrue(e.getMessage().contains("writing"), e.getMessage());
        Assert.assertTrue(e.getCause() instanceof IOException);
    }

    @Test
    public void testStreamAccessorsWrapTheGivenStreams() throws IOException {
        final ReadEndsForMarkDuplicatesCodec codec = new ReadEndsForMarkDuplicatesCodec();
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        codec.setOutputStream(bytes);
        final ReadEndsForMarkDuplicates ends = fill(new ReadEndsForMarkDuplicates(), Values.MIXED, ReadEnds.RR);

        // The disk-based read ends map interleaves a key with each record this way.
        codec.getOutputStream().writeUTF("key1");
        codec.encode(ends);
        codec.getOutputStream().writeInt(42);

        codec.setInputStream(new ByteArrayInputStream(bytes.toByteArray()));
        Assert.assertEquals(codec.getInputStream().readUTF(), "key1");
        assertDecoded(codec.decode(), ends, true);
        Assert.assertEquals(codec.getInputStream().readInt(), 42);
        Assert.assertNull(codec.decode());
    }

    @Test
    public void testCloneKeepsForceSerialization() {
        final ReadEndsForMarkDuplicates fragment = fill(new ReadEndsForMarkDuplicates(), Values.MIXED, ReadEnds.F);

        final ReadEndsForMarkDuplicatesCodec plain = new ReadEndsForMarkDuplicatesCodec();
        final SortingCollection.Codec<ReadEndsForMarkDuplicates> plainClone = plain.clone();
        Assert.assertNotSame(plainClone, plain);
        assertDecoded(roundTrip(plain, (ReadEndsForMarkDuplicatesCodec) plainClone, fragment), fragment, false);

        final ReadEndsForMarkDuplicatesCodec forced = new ReadEndsForMarkDuplicatesCodec();
        forced.setForceSerialization();
        final ReadEndsForMarkDuplicatesCodec forcedClone = (ReadEndsForMarkDuplicatesCodec) forced.clone();
        assertDecoded(roundTrip(forced, forcedClone, fragment), fragment, true);

        // The clone's setting is its own.
        forced.unsetForceSerialization();
        assertDecoded(roundTrip(forcedClone, (ReadEndsForMarkDuplicatesCodec) forcedClone.clone(), fragment), fragment, true);
        assertDecoded(roundTrip(forced, (ReadEndsForMarkDuplicatesCodec) forced.clone(), fragment), fragment, false);
    }

    @Test
    public void testUnsetForceSerialization() {
        final ReadEndsForMarkDuplicates fragment = fill(new ReadEndsForMarkDuplicates(), Values.MIXED, ReadEnds.R);
        final ReadEndsForMarkDuplicatesCodec codec = new ReadEndsForMarkDuplicatesCodec();
        codec.setForceSerialization();
        codec.setForceSerialization();
        assertDecoded(roundTrip(codec, codec, fragment), fragment, true);
        codec.unsetForceSerialization();
        assertDecoded(roundTrip(codec, codec, fragment), fragment, false);
    }

    @DataProvider(name = "spills")
    public Object[][] spills() {
        return new Object[][]{{false}, {true}};
    }

    /** Spills to disk every three records and reads back through clones of the codec. */
    @Test(dataProvider = "spills")
    public void testSortingCollectionSpillRoundTrip(final boolean force) {
        final ReadEndsForMarkDuplicatesCodec codec = new ReadEndsForMarkDuplicatesCodec();
        if (force) codec.setForceSerialization();
        final File tmpDir = IOUtil.createTempDir("readEndsCodecSpill").toFile();
        try {
            assertSpillRoundTrip(codec, tmpDir, force);
        } finally {
            IOUtil.deleteDirectoryTree(tmpDir);
        }
    }

    private static void assertSpillRoundTrip(final ReadEndsForMarkDuplicatesCodec codec, final File tmpDir, final boolean force) {
        final SortingCollection<ReadEndsForMarkDuplicates> sorter = SortingCollection.newInstance(ReadEndsForMarkDuplicates.class,
                codec, Comparator.comparingLong(e -> e.read1IndexInFile), 3, tmpDir);
        final List<ReadEndsForMarkDuplicates> expected = new ArrayList<>();
        for (int i = 0; i < 11; i++) {
            final ReadEndsForMarkDuplicates ends = fill(new ReadEndsForMarkDuplicates(), Values.MIXED, ORIENTATIONS[i % ORIENTATIONS.length]);
            ends.read1IndexInFile = i;
            ends.read1Coordinate = 500 - i;
            ends.read2Coordinate = 700 + i;
            ends.score = (short) (i * 3);
            expected.add(ends);
        }
        // added in reverse so that every spill file holds out-of-order records
        for (int i = expected.size() - 1; i >= 0; i--) {
            sorter.add(expected.get(i));
        }
        sorter.doneAdding();

        int i = 0;
        try (final CloseableIterator<ReadEndsForMarkDuplicates> it = sorter.iterator()) {
            while (it.hasNext()) {
                final ReadEndsForMarkDuplicates actual = it.next();
                final ReadEndsForMarkDuplicates exp = expected.get(i++);
                assertDecoded(actual, exp, force || exp.orientation >= ReadEnds.FF);
                Assert.assertEquals(actual.duplicateSetSize, exp.duplicateSetSize);
            }
        }
        Assert.assertEquals(i, expected.size());
        sorter.cleanup();
    }

    static final class FailingOutputStream extends OutputStream {
        @Override
        public void write(final int b) throws IOException {
            throw new IOException("disk full");
        }
    }
}
