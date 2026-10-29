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

package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SortingCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadCodec;
import picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequence;
import picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequenceWithBarcodes;
import picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadWithBarcodesCodec;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.Arrays;

/**
 * Tests of the pieces of {@link EstimateLibraryComplexity} below the command line: the PairedReadSequence classes,
 * their spill-to-disk codecs, the split of reads into hash parts, barcode values, the read quality filter, and
 * command-line validation.
 */
public class EstimateLibraryComplexityPairedReadTest {

    private static PairedReadSequence sequence(final short readGroup, final short tile, final int x, final int y,
                                               final String read1, final String read2) {
        final PairedReadSequence prs = new PairedReadSequence();
        prs.setReadGroup(readGroup);
        prs.setTile(tile);
        prs.setX(x);
        prs.setY(y);
        prs.read1 = read1.getBytes();
        prs.read2 = read2.getBytes();
        return prs;
    }

    private static byte[] encode(final SortingCollection.Codec<PairedReadSequence> codec, final PairedReadSequence... values) {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        codec.setOutputStream(bytes);
        for (final PairedReadSequence value : values) {
            codec.encode(value);
        }
        return bytes.toByteArray();
    }

    @Test
    public void testPairedReadCodecRoundTrip() {
        final PairedReadSequence first = sequence((short) 3, (short) 1101, 70_000, 123_456, "ACGTN", "GG");
        final PairedReadSequence second = sequence((short) -1, (short) -1, -1, -1, "T", "CA");
        final byte[] encoded = encode(new PairedReadCodec(), first, second);

        final PairedReadCodec decoder = new PairedReadCodec();
        decoder.setInputStream(new ByteArrayInputStream(encoded));
        for (final PairedReadSequence original : new PairedReadSequence[]{first, second}) {
            final PairedReadSequence decoded = decoder.decode();
            Assert.assertEquals(decoded.getReadGroup(), original.getReadGroup());
            Assert.assertEquals(decoded.getTile(), original.getTile());
            Assert.assertEquals(decoded.getX(), original.getX());
            Assert.assertEquals(decoded.getY(), original.getY());
            Assert.assertEquals(decoded.read1, original.read1);
            Assert.assertEquals(decoded.read2, original.read2);
        }
        Assert.assertNull(decoder.decode());
    }

    @Test
    public void testPairedReadWithBarcodesCodecRoundTrip() {
        final PairedReadSequenceWithBarcodes value = new PairedReadSequenceWithBarcodes(
                sequence((short) 7, (short) 2202, 4, 5, "CCA", "TTTT"));
        value.barcode = "AAAA".hashCode();
        value.readOneBarcode = -42;
        value.readTwoBarcode = Integer.MAX_VALUE;
        final byte[] encoded = encode(new PairedReadWithBarcodesCodec(), value);

        final PairedReadWithBarcodesCodec decoder = new PairedReadWithBarcodesCodec();
        decoder.setInputStream(new ByteArrayInputStream(encoded));
        final PairedReadSequence decoded = decoder.decode();
        Assert.assertTrue(decoded instanceof PairedReadSequenceWithBarcodes);
        final PairedReadSequenceWithBarcodes withBarcodes = (PairedReadSequenceWithBarcodes) decoded;
        Assert.assertEquals(withBarcodes.getReadGroup(), 7);
        Assert.assertEquals(withBarcodes.getTile(), 2202);
        Assert.assertEquals(withBarcodes.getX(), 4);
        Assert.assertEquals(withBarcodes.getY(), 5);
        Assert.assertEquals(new String(withBarcodes.read1), "CCA");
        Assert.assertEquals(new String(withBarcodes.read2), "TTTT");
        Assert.assertEquals(withBarcodes.barcode, "AAAA".hashCode());
        Assert.assertEquals(withBarcodes.readOneBarcode, -42);
        Assert.assertEquals(withBarcodes.readTwoBarcode, Integer.MAX_VALUE);
        Assert.assertNull(decoder.decode());
    }

    @Test
    public void testBarcodesCodecRejectsSequencesWithoutBarcodes() {
        final PairedReadWithBarcodesCodec codec = new PairedReadWithBarcodesCodec();
        codec.setOutputStream(new ByteArrayOutputStream());
        try {
            codec.encode(sequence((short) 0, (short) 1, 1, 1, "A", "C"));
            Assert.fail("expected a PicardException");
        } catch (final PicardException e) {
            Assert.assertTrue(e.getMessage().contains("not a PairedReadSequenceWithBarcodes"), e.getMessage());
        }
    }

    /** The factory codec and the clones that SortingCollection reads its spill files with agree on the format. */
    @Test
    public void testCodecFactoriesAndClones() {
        final PairedReadSequence plainValue = sequence((short) 4, (short) 1201, 6, 7, "GATTACA", "CAT");
        final SortingCollection.Codec<PairedReadSequence> factory = PairedReadSequence.getCodec();
        final SortingCollection.Codec<PairedReadSequence> plainClone = factory.clone();
        Assert.assertNotSame(plainClone, factory);
        plainClone.setInputStream(new ByteArrayInputStream(encode(factory, plainValue)));
        final PairedReadSequence plainDecoded = plainClone.decode();
        Assert.assertEquals(plainDecoded.getReadGroup(), 4);
        Assert.assertEquals(plainDecoded.getTile(), 1201);
        Assert.assertEquals(plainDecoded.getX(), 6);
        Assert.assertEquals(plainDecoded.getY(), 7);
        Assert.assertEquals(new String(plainDecoded.read1), "GATTACA");
        Assert.assertEquals(new String(plainDecoded.read2), "CAT");
        Assert.assertNull(plainClone.decode());

        final PairedReadSequenceWithBarcodes barcodedValue = new PairedReadSequenceWithBarcodes(plainValue);
        barcodedValue.barcode = 11;
        barcodedValue.readOneBarcode = -12;
        barcodedValue.readTwoBarcode = 13;
        final PairedReadWithBarcodesCodec barcodes = new PairedReadWithBarcodesCodec();
        final SortingCollection.Codec<PairedReadSequence> barcodesClone = barcodes.clone();
        Assert.assertNotSame(barcodesClone, barcodes);
        barcodesClone.setInputStream(new ByteArrayInputStream(encode(barcodes, barcodedValue)));
        final PairedReadSequenceWithBarcodes barcodedDecoded = (PairedReadSequenceWithBarcodes) barcodesClone.decode();
        Assert.assertEquals(barcodedDecoded.getTile(), 1201);
        Assert.assertEquals(new String(barcodedDecoded.read1), "GATTACA");
        Assert.assertEquals(barcodedDecoded.barcode, 11);
        Assert.assertEquals(barcodedDecoded.readOneBarcode, -12);
        Assert.assertEquals(barcodedDecoded.readTwoBarcode, 13);
        Assert.assertNull(barcodesClone.decode());
    }

    @Test
    public void testCopyIntoSequenceWithBarcodes() {
        final PairedReadSequence source = sequence((short) 2, (short) 1102, 30, 40, "ACGT", "TTAA");
        source.qualityOk = false;
        source.setLibraryId((short) 9);

        final PairedReadSequenceWithBarcodes copy = new PairedReadSequenceWithBarcodes(source);
        Assert.assertEquals(copy.getReadGroup(), 2);
        Assert.assertEquals(copy.getTile(), 1102);
        Assert.assertEquals(copy.getX(), 30);
        Assert.assertEquals(copy.getY(), 40);
        Assert.assertFalse(copy.qualityOk);
        Assert.assertEquals(copy.getLibraryId(), 9);
        Assert.assertEquals(copy.read1, source.read1);
        Assert.assertEquals(copy.read2, source.read2);

        // the bases are copies, not shared arrays
        source.read1[0] = 'N';
        source.read2[0] = 'N';
        Assert.assertEquals(new String(copy.read1), "ACGT");
        Assert.assertEquals(new String(copy.read2), "TTAA");
    }

    @Test
    public void testCopyOfNullIsRejected() {
        try {
            new PairedReadSequenceWithBarcodes(null);
            Assert.fail("expected a PicardException");
        } catch (final PicardException e) {
            Assert.assertTrue(e.getMessage().contains("val was null"), e.getMessage());
        }
    }

    /**
     * A new pair has not failed the quality filter, which each read can only clear, has no location until its name
     * parses, and has a read group outside the header's until one is set, which keeps a pair without one apart from
     * read-grouped pairs when optical duplicates are counted.
     */
    @Test
    public void testNewSequenceDefaults() {
        final PairedReadSequence prs = new PairedReadSequence();
        Assert.assertEquals(prs.getReadGroup(), -1);
        Assert.assertTrue(prs.qualityOk);
        Assert.assertFalse(prs.hasLocation());
    }

    @DataProvider(name = "hashParts")
    public Object[][] hashParts() {
        // 5 parts after skipping 1 base, up to base 15: base i (1 <= i < 15) belongs to part (i - 1) % 5
        return new Object[][]{
                {0, new boolean[]{false, false, false, false, false}},
                {1, new boolean[]{true, false, false, false, false}},
                {6, new boolean[]{true, false, false, false, false}},
                {12, new boolean[]{false, true, false, false, false}},
                {14, new boolean[]{false, false, false, true, false}},
                {5, new boolean[]{false, false, false, false, true}},
                {15, new boolean[]{false, false, false, false, false}},
                {19, new boolean[]{false, false, false, false, false}},
        };
    }

    @Test(dataProvider = "hashParts")
    public void testHashPartsCoverInterleavedBases(final int changedBase, final boolean[] partChanged) {
        final String bases = "ACGTACGTACGTACGTACGT";
        final char[] changed = bases.toCharArray();
        changed[changedBase] = changed[changedBase] == 'A' ? 'C' : 'A';

        final PairedReadSequence reference = sequence((short) 0, (short) 1, 1, 1, bases, bases);
        final PairedReadSequence read1Changed = sequence((short) 0, (short) 1, 1, 1, new String(changed), bases);
        final PairedReadSequence read2Changed = sequence((short) 0, (short) 1, 1, 1, bases, new String(changed));
        for (final PairedReadSequence prs : Arrays.asList(reference, read1Changed, read2Changed)) {
            prs.initHashes(5, 1, 15);
            Assert.assertEquals(prs.hashes1.length, 5);
            Assert.assertEquals(prs.hashes2.length, 5);
        }
        for (int part = 0; part < 5; part++) {
            Assert.assertEquals(read1Changed.hashes1[part] != reference.hashes1[part], partChanged[part], "read 1 part " + part);
            Assert.assertEquals(read1Changed.hashes2[part], reference.hashes2[part], "read 2 part " + part);
            Assert.assertEquals(read2Changed.hashes2[part] != reference.hashes2[part], partChanged[part], "read 2 part " + part);
            Assert.assertEquals(read2Changed.hashes1[part], reference.hashes1[part], "read 1 part " + part);
        }
    }

    @Test
    public void testReadBarcodeValues() {
        final SAMRecord rec = new SAMRecord(new SAMFileHeader());
        rec.setAttribute("BC", "ACGT");
        rec.setAttribute("RX", "");

        Assert.assertEquals(EstimateLibraryComplexity.getReadBarcodeValue(rec, null), 0);
        Assert.assertEquals(EstimateLibraryComplexity.getReadBarcodeValue(rec, "BX"), 0);
        Assert.assertEquals(EstimateLibraryComplexity.getReadBarcodeValue(rec, "BC"), "ACGT".hashCode());
        Assert.assertEquals(EstimateLibraryComplexity.getReadBarcodeValue(rec, "RX"), "".hashCode());

        final EstimateLibraryComplexity elc = new EstimateLibraryComplexity();
        Assert.assertEquals(elc.getBarcodeValue(rec), 0);
        elc.BARCODE_TAG = "BC";
        Assert.assertEquals(elc.getBarcodeValue(rec), "ACGT".hashCode());
    }

    private static byte[] quals(final int... values) {
        final byte[] out = new byte[values.length];
        for (int i = 0; i < values.length; i++) {
            out[i] = (byte) values[i];
        }
        return out;
    }

    private static byte[] flat(final int length, final int q) {
        final byte[] out = new byte[length];
        Arrays.fill(out, (byte) q);
        return out;
    }

    @DataProvider(name = "qualityCheck")
    public Object[][] qualityCheck() {
        return new Object[][]{
                // {bases, quals, seed length, min quality, MAX_READ_LENGTH, passes}
                {"ACGTACGTAC", flat(10, 20), 5, 20, 0, true},
                {"ACGTA", flat(5, 20), 5, 20, 0, true},
                {"ACGT", flat(4, 40), 5, 20, 0, false},
                {"ACGTNACGTA", flat(10, 40), 5, 20, 0, false},
                {"NCGTAACGTA", flat(10, 40), 5, 20, 0, false},
                {"ACnTAACGTA", flat(10, 40), 5, 20, 0, false},
                {".CGTAACGTA", flat(10, 40), 5, 20, 0, false},
                {"ACGTANCGTA", flat(10, 40), 5, 20, 0, true},
                {"ACGTANCGTA", flat(10, 40), 6, 20, 0, false},
                // the mean is truncated to an integer: 199 / 10 = 19
                {"ACGTACGTAC", quals(20, 20, 20, 20, 20, 20, 20, 20, 20, 19), 5, 20, 0, false},
                {"ACGTACGTAC", quals(20, 20, 20, 20, 20, 20, 20, 20, 20, 19), 5, 19, 0, true},
                {"ACGTACGTAC", flat(10, 0), 5, 0, 0, true},
                // MAX_READ_LENGTH limits the bases averaged
                {"ACGTACGTAC", quals(30, 30, 30, 30, 30, 0, 0, 0, 0, 0), 5, 20, 0, false},
                {"ACGTACGTAC", quals(30, 30, 30, 30, 30, 0, 0, 0, 0, 0), 5, 20, 5, true},
                {"ACGTACGTAC", quals(30, 30, 30, 30, 30, 0, 0, 0, 0, 0), 5, 20, 6, true},
                {"ACGTACGTAC", quals(30, 30, 30, 30, 30, 0, 0, 0, 0, 0), 5, 20, 7, true},
                {"ACGTACGTAC", quals(30, 30, 30, 30, 30, 0, 0, 0, 0, 0), 5, 20, 8, false},
                {"ACGTACGTAC", quals(30, 30, 30, 30, 30, 0, 0, 0, 0, 0), 5, 20, 100, false},
                {"ACGTACGTAC", quals(30, 30, 30, 30, 30, 0, 0, 0, 0, 0), 5, 20, -1, false},
        };
    }

    /**
     * MAX_READ_LENGTH limits the bases whose mean quality is checked, and a negative value means no maximum (the rows
     * with a MAX_READ_LENGTH). The documentation does not specify this; the expected value follows the current
     * implementation.
     */
    @Test(dataProvider = "qualityCheck")
    public void testPassesQualityCheck(final String bases, final byte[] quals, final int seedLength, final int minQuality,
                                       final int maxReadLength, final boolean passes) {
        final EstimateLibraryComplexity elc = new EstimateLibraryComplexity();
        elc.MAX_READ_LENGTH = maxReadLength;
        Assert.assertEquals(elc.passesQualityCheck(bases.getBytes(), quals, seedLength, minQuality), passes);
    }

    private static final String MAX_READ_LENGTH_ERROR = "MAX_READ_LENGTH must be greater than MIN_IDENTICAL_BASES";
    private static final String MIN_IDENTICAL_BASES_ERROR = "MIN_IDENTICAL_BASES must be greater than 0";

    @DataProvider(name = "validation")
    public Object[][] validation() {
        return new Object[][]{
                {5, 0, null},
                {5, 6, null},
                {5, -3, null},
                {1, 0, null},
                {5, 4, new String[]{MAX_READ_LENGTH_ERROR}},
                {5, 1, new String[]{MAX_READ_LENGTH_ERROR}},
                {0, 0, new String[]{MIN_IDENTICAL_BASES_ERROR}},
                {0, 5, new String[]{MIN_IDENTICAL_BASES_ERROR}},
                {-2, 0, new String[]{MIN_IDENTICAL_BASES_ERROR}},
        };
    }

    /**
     * A negative MAX_READ_LENGTH passes validation, like 0, which means no maximum (the third row). The documentation
     * does not specify this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "validation")
    public void testCustomCommandLineValidation(final int minIdenticalBases, final int maxReadLength, final String[] errors) {
        final EstimateLibraryComplexity elc = new EstimateLibraryComplexity();
        elc.MIN_IDENTICAL_BASES = minIdenticalBases;
        elc.MAX_READ_LENGTH = maxReadLength;
        Assert.assertEquals(elc.customCommandLineValidation(), errors);
    }
}
