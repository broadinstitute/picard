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
import picard.sam.util.RepresentativeReadIndexer;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Tests for {@link RepresentativeReadIndexerCodec}: round trips of the read index, set size and representative read
 * index at extreme values, several records in one stream, the end of a stream, truncated records, write failures,
 * {@code clone()} and spilling through a {@link SortingCollection}.
 */
public class RepresentativeReadIndexerCodecTest {

    private static RepresentativeReadIndexer indexer(final int readIndex, final int setSize, final int representativeIndex) {
        final RepresentativeReadIndexer rri = new RepresentativeReadIndexer();
        rri.readIndexInFile = readIndex;
        rri.setSize = setSize;
        rri.representativeReadIndexInFile = representativeIndex;
        return rri;
    }

    private static byte[] encode(final RepresentativeReadIndexerCodec codec, final RepresentativeReadIndexer... records) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        codec.setOutputStream(bytes);
        for (final RepresentativeReadIndexer rri : records) {
            codec.encode(rri);
        }
        return bytes.toByteArray();
    }

    private static void assertIndexer(final RepresentativeReadIndexer actual, final RepresentativeReadIndexer expected) {
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.readIndexInFile, expected.readIndexInFile, "readIndexInFile");
        Assert.assertEquals(actual.setSize, expected.setSize, "setSize");
        Assert.assertEquals(actual.representativeReadIndexInFile, expected.representativeReadIndexInFile, "representativeReadIndexInFile");
    }

    @DataProvider(name = "values")
    public Object[][] values() {
        return new Object[][]{
                {0, 0, 0},
                {-1, -1, -1},
                {Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE},
                {Integer.MIN_VALUE, Integer.MIN_VALUE, Integer.MIN_VALUE},
                {Integer.MAX_VALUE, 2, Integer.MIN_VALUE},
                {0x01020304, 0x05060708, 0x090A0B0C}
        };
    }

    @Test(dataProvider = "values")
    public void testRoundTrip(final int readIndex, final int setSize, final int representativeIndex) {
        final RepresentativeReadIndexer rri = indexer(readIndex, setSize, representativeIndex);
        final RepresentativeReadIndexerCodec codec = new RepresentativeReadIndexerCodec();
        final byte[] bytes = encode(codec, rri);

        final RepresentativeReadIndexerCodec decoder = (RepresentativeReadIndexerCodec) codec.clone();
        decoder.setInputStream(new ByteArrayInputStream(bytes));
        assertIndexer(decoder.decode(), rri);
        Assert.assertNull(decoder.decode());
    }

    @Test
    public void testManyRecordsInOneStream() {
        final RepresentativeReadIndexer[] records = {indexer(5, 3, 2), indexer(6, 1, 6), indexer(7, 3, 2), indexer(Integer.MAX_VALUE, 2, 0)};
        final RepresentativeReadIndexerCodec codec = new RepresentativeReadIndexerCodec();
        final byte[] bytes = encode(codec, records);
        codec.setInputStream(new ByteArrayInputStream(bytes));
        for (final RepresentativeReadIndexer expected : records) {
            assertIndexer(codec.decode(), expected);
        }
        Assert.assertNull(codec.decode());
        Assert.assertNull(codec.decode());
    }

    @Test
    public void testDecodeEmptyStreamReturnsNull() {
        final RepresentativeReadIndexerCodec codec = new RepresentativeReadIndexerCodec();
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
        final RepresentativeReadIndexerCodec codec = new RepresentativeReadIndexerCodec();
        final byte[] bytes = encode(codec, indexer(1, 2, 3));
        codec.setInputStream(new ByteArrayInputStream(Arrays.copyOf(bytes, keepHalf ? bytes.length / 2 : bytes.length - 1)));
        final PicardException e = Assert.expectThrows(PicardException.class, codec::decode);
        Assert.assertTrue(e.getCause() instanceof EOFException, String.valueOf(e.getCause()));
    }

    @Test
    public void testEncodeFailureIsWrapped() {
        final RepresentativeReadIndexerCodec codec = new RepresentativeReadIndexerCodec();
        codec.setOutputStream(new ReadEndsForMarkDuplicatesCodecTest.FailingOutputStream());
        final PicardException e = Assert.expectThrows(PicardException.class, () -> codec.encode(indexer(1, 2, 3)));
        Assert.assertTrue(e.getCause() instanceof IOException);
    }

    @Test
    public void testStreamAccessors() throws IOException {
        final RepresentativeReadIndexerCodec codec = new RepresentativeReadIndexerCodec();
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        codec.setOutputStream(bytes);
        codec.getOutputStream().writeInt(9);
        codec.encode(indexer(1, 2, 3));

        codec.setInputStream(new ByteArrayInputStream(bytes.toByteArray()));
        Assert.assertEquals(codec.getInputStream().readInt(), 9);
        assertIndexer(codec.decode(), indexer(1, 2, 3));
    }

    @Test
    public void testClone() {
        final RepresentativeReadIndexerCodec codec = new RepresentativeReadIndexerCodec();
        final SortingCollection.Codec<RepresentativeReadIndexer> clone = codec.clone();
        Assert.assertNotSame(clone, codec);
        final byte[] bytes = encode((RepresentativeReadIndexerCodec) clone, indexer(4, 2, 3));
        codec.setInputStream(new ByteArrayInputStream(bytes));
        assertIndexer(codec.decode(), indexer(4, 2, 3));
        Assert.assertNull(codec.decode());
    }

    /** MarkDuplicates sorts these by read index and spills them to disk. */
    @Test
    public void testSortingCollectionSpillRoundTrip() {
        final File tmpDir = IOUtil.createTempDir("representativeReadIndexerSpill").toFile();
        try {
            final SortingCollection<RepresentativeReadIndexer> sorter = SortingCollection.newInstance(RepresentativeReadIndexer.class,
                    new RepresentativeReadIndexerCodec(), Comparator.comparing(rri -> rri.readIndexInFile), 2, tmpDir);
            final int[] readIndices = {9, 3, 7, 0, 5, 1, 8};
            for (final int readIndex : readIndices) {
                sorter.add(indexer(readIndex, readIndex + 1, readIndex / 2));
            }
            sorter.doneAdding();

            final int[] sorted = readIndices.clone();
            Arrays.sort(sorted);
            int i = 0;
            try (final CloseableIterator<RepresentativeReadIndexer> it = sorter.iterator()) {
                while (it.hasNext()) {
                    final int readIndex = sorted[i++];
                    assertIndexer(it.next(), indexer(readIndex, readIndex + 1, readIndex / 2));
                }
            }
            Assert.assertEquals(i, sorted.length);
            sorter.cleanup();
        } finally {
            IOUtil.deleteDirectoryTree(tmpDir);
        }
    }
}
