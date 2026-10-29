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
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.util.LibraryIdGenerator;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesWithBarcodes;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static picard.sam.markduplicates.MarkDuplicatesQueryNameTestSupport.*;

/**
 * Tests the barcode options of MarkDuplicates on query-grouped input: BARCODE_TAG (including reads missing the
 * tag), READ_ONE_BARCODE_TAG and READ_TWO_BARCODE_TAG (including a group whose read two comes first), DUPLEX_UMI
 * strand normalization, and the MOLECULAR_IDENTIFIER_TAG strand suffixes; also the barcode fields that
 * buildReadEnds fills through the MarkDuplicatesHelper interface.
 */
public class MarkDuplicatesQueryGroupedBarcodeTest {

    private static final String NO_OPTICAL = "READ_NAME_REGEX=null";

    @DataProvider(name = "barcodeTag")
    public Object[][] barcodeTag() {
        final List<Object[]> rows = new ArrayList<>();
        for (final String[] spill : new String[][]{new String[0], spillArgs()}) {
            rows.add(new Object[]{args(spill, "BARCODE_TAG=RX"), "..DD....DDD.D",
                    Arrays.asList(0, 0, 0, 0, null, null, 6, 6, 6, 6), Arrays.asList(2, 2, 2, 2, null, null, 2, 2, 2, 2), 2L, 2L});
            rows.add(new Object[]{spill, "..DDDDDDDDDDD",
                    Arrays.asList(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), Arrays.asList(5, 5, 5, 5, 5, 5, 5, 5, 5, 5), 4L, 3L});
        }
        return rows.toArray(new Object[0][]);
    }

    /**
     * With BARCODE_TAG=RX, pairs at 100F/309R split into sets by RX: {b1, b2} with AAAA, b3 alone with CCCC, and
     * {b4, b5}, which both lack RX and are compared as having the same missing barcode. A fragment at 100F is a
     * duplicate only when a pair with its barcode ends there: b6 (AAAA) and b8 (no RX) are, b7 (GGGG) is not.
     * Without BARCODE_TAG all five pairs form one set and every fragment is a duplicate. The DI/DS tags of the
     * pair records follow the same sets; the fragments' tags are not checked.
     */
    @Test(dataProvider = "barcodeTag")
    public void testBarcodeTagWithMissingBarcodes(final String[] extraArgs, final String expectedFlags,
                                                  final List<Object> expectedDi, final List<Object> expectedDs,
                                                  final long pairDuplicates, final long unpairedDuplicates)
            throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("b1", "rgA", 100, false, 300, true, 'I', "RX:Z:AAAA")
                .pair("b2", "rgA", 100, false, 300, true, '5', "RX:Z:AAAA")
                .pair("b3", "rgA", 100, false, 300, true, '5', "RX:Z:CCCC")
                .pair("b4", "rgA", 100, false, 300, true, '?')
                .pair("b5", "rgA", 100, false, 300, true, '5')
                .fragment("b6", "rgA", 100, false, 'I', "RX:Z:AAAA")
                .fragment("b7", "rgA", 100, false, 'I', "RX:Z:GGGG")
                .fragment("b8", "rgA", 100, false, '5');

        final Result result = run(input.write(), args(extraArgs, NO_OPTICAL, "TAG_DUPLICATE_SET_MEMBERS=true"));

        assertDuplicateFlags(result, expectedFlags);
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_SET_INDEX_TAG).subList(0, 10), expectedDi);
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_SET_SIZE_TAG).subList(0, 10), expectedDs);
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
        assertMetrics(result.metricsFor("libA"), 3, 5, 0, 0, unpairedDuplicates, pairDuplicates, 0);
    }

    @DataProvider(name = "readBarcodeTags")
    public Object[][] readBarcodeTags() {
        final List<Object[]> rows = new ArrayList<>();
        for (final String[] spill : new String[][]{new String[0], spillArgs()}) {
            rows.add(new Object[]{args(spill, "READ_ONE_BARCODE_TAG=BX"), "..DD..DD.D.."});
            rows.add(new Object[]{args(spill, "READ_TWO_BARCODE_TAG=BX"), "....DDDD.DDD"});
            rows.add(new Object[]{args(spill, "READ_ONE_BARCODE_TAG=BX", "READ_TWO_BARCODE_TAG=BX"), "......DD.D.."});
            rows.add(new Object[]{spill, "..DDDDDD.DDD"});
        }
        return rows.toArray(new Object[0][]);
    }

    /**
     * READ_ONE_BARCODE_TAG is read from read one of a pair and from unpaired reads; READ_TWO_BARCODE_TAG only from
     * read two. The pairs (all at 100F/309R) carry BX on read one / read two of AAA/CCC (r1), AAA/GGG (r2),
     * TTT/CCC (r3) and AAA/CCC (r4, read two written first). The fragments at 2000F carry AAA (r5, r6), CCC (r7)
     * and none (r8). r1 and r5 score highest.
     */
    @Test(dataProvider = "readBarcodeTags")
    public void testReadOneAndReadTwoBarcodeTags(final String[] extraArgs, final String expectedFlags)
            throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("r1", "rgA", 100, false, 300, true, 'I', 'I', false, new String[]{"BX:Z:AAA"}, new String[]{"BX:Z:CCC"})
                .pair("r2", "rgA", 100, false, 300, true, '5', '5', false, new String[]{"BX:Z:AAA"}, new String[]{"BX:Z:GGG"})
                .pair("r3", "rgA", 100, false, 300, true, '5', '5', false, new String[]{"BX:Z:TTT"}, new String[]{"BX:Z:CCC"})
                .pair("r4", "rgA", 100, false, 300, true, '?', '?', true, new String[]{"BX:Z:AAA"}, new String[]{"BX:Z:CCC"})
                .fragment("r5", "rgA", 2000, false, 'I', "BX:Z:AAA")
                .fragment("r6", "rgA", 2000, false, '5', "BX:Z:AAA")
                .fragment("r7", "rgA", 2000, false, '5', "BX:Z:CCC")
                .fragment("r8", "rgA", 2000, false, '5');

        final Result result = run(input.write(), args(extraArgs, NO_OPTICAL, "ASSUME_SORT_ORDER=queryname"));

        assertDuplicateFlags(result, expectedFlags);
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
    }

    @DataProvider(name = "duplexUmi")
    public Object[][] duplexUmi() {
        return new Object[][]{
                {new String[]{"DUPLEX_UMI=true"}, "DD....", new String[]{"/A", "/B", "/A"}},
                {spillArgs(), "....DD", new String[]{"", "", ""}},
                {args(spillArgs(), "DUPLEX_UMI=true"), "DD....", new String[]{"/A", "/B", "/A"}},
                {new String[0], "....DD", new String[]{"", "", ""}},
        };
    }

    /**
     * x1 (read one at 100F, AAA-CCC) is a top strand pair; x2 (read one at 300R, CCC-AAA) is a bottom strand pair
     * whose UMI normalizes to AAA-CCC; x3 (read one at 100F, CCC-AAA) is a top strand pair. All three end at
     * 100F/309R and x2 scores highest. With DUPLEX_UMI x1 is a duplicate of x2; without it x3 is, since x2 and x3
     * carry the same raw UMI. MOLECULAR_IDENTIFIER_TAG gives both reads of a pair the same value, suffixed /A for
     * top and /B for bottom strand pairs under DUPLEX_UMI and unsuffixed otherwise; without DUPLEX_UMI the
     * duplicate x3 shares x2's value.
     */
    @Test(dataProvider = "duplexUmi")
    public void testDuplexUmiAndMolecularIdentifier(final String[] extraArgs, final String expectedFlags,
                                                    final String[] expectedSuffixes) throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("x1", "rgA", 100, false, 300, true, '5', "RX:Z:AAA-CCC")
                .pair("x2", "rgA", 300, true, 100, false, 'I', "RX:Z:CCC-AAA")
                .pair("x3", "rgA", 100, false, 300, true, '?', "RX:Z:CCC-AAA");

        final Result result = run(input.write(),
                args(extraArgs, NO_OPTICAL, "BARCODE_TAG=RX", "MOLECULAR_IDENTIFIER_TAG=MI"));

        assertDuplicateFlags(result, expectedFlags);
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
        final List<Object> mi = result.tagValues("MI");
        final List<String> bases = new ArrayList<>();
        for (int pair = 0; pair < 3; pair++) {
            final String first = (String) mi.get(2 * pair);
            Assert.assertNotNull(first);
            Assert.assertEquals(mi.get(2 * pair + 1), first, "both reads of a pair share the molecular identifier");
            if (expectedSuffixes[pair].isEmpty()) {
                Assert.assertFalse(first.endsWith("/A") || first.endsWith("/B"), first);
                bases.add(first);
            } else {
                Assert.assertTrue(first.endsWith(expectedSuffixes[pair]), first);
                bases.add(first.substring(0, first.length() - 2));
            }
        }
        if (expectedSuffixes[0].isEmpty()) {
            Assert.assertEquals(bases.get(2), bases.get(1), "duplicates share the molecular identifier");
        } else {
            Assert.assertEquals(bases.get(0), bases.get(1), "the strands of one molecule share a base identifier");
        }
    }

    /** Without MOLECULAR_IDENTIFIER_TAG, BARCODE_TAG adds no molecular identifier tag. */
    @Test
    public void testNoMolecularIdentifierByDefault() throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("x1", "rgA", 100, false, 300, true, '5', "RX:Z:AAA-CCC")
                .pair("x2", "rgA", 300, true, 100, false, 'I', "RX:Z:CCC-AAA");

        final Result result = run(input.write(), NO_OPTICAL, "BARCODE_TAG=RX", "DUPLEX_UMI=true");

        assertDuplicateFlags(result, "DD..");
        Assert.assertEquals(result.tagValues("MI"), Arrays.asList(null, null, null, null));
    }

    /**
     * buildReadEnds, called through the MarkDuplicatesHelper interface, records the index it is given (the first
     * index of the query group), takes the barcode from BARCODE_TAG on every read, READ_ONE_BARCODE_TAG only from
     * read one and unpaired reads, and READ_TWO_BARCODE_TAG only from read two. Fields are compared between reads
     * rather than to fixed hash values.
     */
    @Test
    public void testBuildReadEndsBarcodeFields() throws IOException {
        final File file = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("p", "rgA", 100, false, 300, true, 'I', "RX:Z:AAA", "B1:Z:GGG", "B2:Z:TTT")
                .fragment("f", "rgA", 500, false, 'I', "RX:Z:AAA", "B1:Z:GGG", "B2:Z:TTT")
                .fragment("g", "rgA", 500, false, 'I')
                .write();
        final SAMFileHeader header;
        final List<SAMRecord> records = new ArrayList<>();
        try (SamReader reader = SamReaderFactory.makeDefault().open(file)) {
            header = reader.getFileHeader();
            reader.iterator().forEachRemaining(records::add);
        }

        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.READ_NAME_REGEX = null;
        markDuplicates.setupOpticalDuplicateFinder();
        markDuplicates.BARCODE_TAG = "RX";
        markDuplicates.READ_ONE_BARCODE_TAG = "B1";
        markDuplicates.READ_TWO_BARCODE_TAG = "B2";
        markDuplicates.libraryIdGenerator = new LibraryIdGenerator(header);
        final MarkDuplicatesHelper helper = markDuplicates;

        final List<ReadEndsForMarkDuplicatesWithBarcodes> ends = new ArrayList<>();
        for (final SAMRecord rec : records) {
            ends.add((ReadEndsForMarkDuplicatesWithBarcodes) helper.buildReadEnds(header, 7, rec, true));
        }
        final ReadEndsForMarkDuplicatesWithBarcodes readOne = ends.get(0);
        final ReadEndsForMarkDuplicatesWithBarcodes readTwo = ends.get(1);
        final ReadEndsForMarkDuplicatesWithBarcodes tagged = ends.get(2);
        final ReadEndsForMarkDuplicatesWithBarcodes untagged = ends.get(3);

        ends.forEach(e -> Assert.assertEquals(e.read1IndexInFile, 7L));
        Assert.assertEquals(readOne.barcode, tagged.barcode);
        Assert.assertEquals(readTwo.barcode, tagged.barcode);
        Assert.assertNotEquals(untagged.barcode, tagged.barcode);

        Assert.assertEquals(readOne.readOneBarcode, tagged.readOneBarcode);
        Assert.assertNotEquals(untagged.readOneBarcode, tagged.readOneBarcode);
        Assert.assertEquals(readTwo.readOneBarcode, untagged.readOneBarcode);

        Assert.assertNotEquals(readTwo.readTwoBarcode, untagged.readTwoBarcode);
        Assert.assertEquals(readOne.readTwoBarcode, untagged.readTwoBarcode);
        Assert.assertEquals(tagged.readTwoBarcode, untagged.readTwoBarcode);
    }
}
