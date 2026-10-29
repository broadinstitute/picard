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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import static picard.sam.markduplicates.MarkDuplicatesQueryNameTestSupport.*;

/**
 * Tests the errors MarkDuplicates raises for query-grouped input whose header or content does not match the order
 * it claims, for conflicting sort-order options, and for UMIs that BARCODE_TAG or DUPLEX_UMI cannot use.
 */
public class MarkDuplicatesQueryGroupedInputErrorsTest {

    private static final String NO_OPTICAL = "READ_NAME_REGEX=null";

    private static SamText groupedPairs(final String hdOrder) {
        return new SamText(hdOrder).readGroup("rgA", "libA")
                .pair("q1", "rgA", 100, false, 300, true, 'I')
                .pair("q2", "rgA", 100, false, 300, true, '5');
    }

    /**
     * A header claiming queryname order is written back with that order, so names out of order make the writer
     * fail; names compare as plain strings, so "q10" may precede "q9" but not follow it.
     */
    @DataProvider(name = "namesOutOfOrder")
    public Object[][] namesOutOfOrder() {
        return new Object[][]{{"q2", "q1"}, {"q9", "q10"}};
    }

    @Test(dataProvider = "namesOutOfOrder")
    public void testQuerynameHeaderWithNamesOutOfOrder(final String firstName, final String secondName)
            throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair(firstName, "rgA", 100, false, 300, true, 'I')
                .pair(secondName, "rgA", 100, false, 300, true, '5');

        final Throwable thrown = runExpectingFailure(input.write(), NO_OPTICAL);

        Assert.assertTrue(thrown instanceof IllegalArgumentException, thrown.toString());
        Assert.assertTrue(thrown.getMessage().contains("out of order"), thrown.getMessage());
    }

    /** The same out-of-order names pass once ASSUME_SORT_ORDER=queryname treats the input as only query grouped. */
    @Test(dataProvider = "namesOutOfOrder")
    public void testAssumeQuerynameAcceptsNamesOutOfOrder(final String firstName, final String secondName)
            throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair(firstName, "rgA", 100, false, 300, true, 'I')
                .pair(secondName, "rgA", 100, false, 300, true, '5');

        final Result result = run(input.write(), NO_OPTICAL, "ASSUME_SORT_ORDER=queryname");

        assertDuplicateFlags(result, "..DD");
        Assert.assertEquals(result.records.get(0).getReadName(), firstName);
    }

    @DataProvider(name = "unusableOrders")
    public Object[][] unusableOrders() {
        return new Object[][]{
                {"SO:unsorted\tGO:query", new String[0], "unsorted"},
                {"SO:unknown\tGO:query", new String[0], "unknown"},
                {"SO:unsorted", new String[0], "unsorted"},
                {"SO:queryname", new String[]{"ASSUME_SORT_ORDER=unsorted"}, "unsorted"},
                {"SO:unsorted\tGO:query", new String[]{"ASSUME_SORT_ORDER=unknown"}, "unknown"},
                {"SO:queryname", new String[]{"ASSUME_SORT_ORDER=duplicate"}, "duplicate"},
        };
    }

    /**
     * A query-grouped header (GO:query) is not enough without ASSUME_SORT_ORDER=queryname, and no sort order other
     * than queryname or coordinate is accepted, from the header or from ASSUME_SORT_ORDER.
     */
    @Test(dataProvider = "unusableOrders")
    public void testUnusableSortOrders(final String hdOrder, final String[] extraArgs, final String expectedMessage)
            throws IOException {
        final Throwable thrown = runExpectingFailure(groupedPairs(hdOrder).write(), args(extraArgs, NO_OPTICAL));

        Assert.assertTrue(thrown instanceof PicardException, thrown.toString());
        Assert.assertTrue(thrown.getMessage().contains(expectedMessage), thrown.getMessage());
    }

    /** ASSUME_SORTED and ASSUME_SORT_ORDER are mutually exclusive, so the command line is rejected. */
    @Test
    public void testAssumeSortedConflictsWithAssumeSortOrder() throws IOException {
        final File input = groupedPairs("SO:unsorted\tGO:query").write();
        final File dir = Files.createTempDirectory("mdQueryName").toFile();
        dir.deleteOnExit();
        final File output = new File(dir, "out.sam");
        final File metrics = new File(dir, "metrics.txt");

        final int status = new MarkDuplicates().instanceMain(new String[]{
                "INPUT=" + input.getAbsolutePath(), "OUTPUT=" + output.getAbsolutePath(),
                "METRICS_FILE=" + metrics.getAbsolutePath(), "ASSUME_SORTED=true", "ASSUME_SORT_ORDER=queryname"});

        Assert.assertEquals(status, 1);
        Assert.assertFalse(output.exists());
        Assert.assertFalse(metrics.exists());
        Assert.assertTrue(dir.delete());
    }

    @DataProvider(name = "badUmis")
    public Object[][] badUmis() {
        return new Object[][]{
                {"AAACCC", true, "X-Y"},
                {"AA-CC-GG", true, "X-Y"},
                {"AAX-CCC", true, "illegal characters"},
                {"AAX", false, "illegal characters"},
                {"AA CC", false, "illegal characters"},
        };
    }

    /** A UMI with characters other than ACGTN and hyphen, or a duplex UMI not made of two hyphen-separated parts. */
    @Test(dataProvider = "badUmis")
    public void testUnusableUmis(final String umi, final boolean duplex, final String expectedMessage)
            throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("q1", "rgA", 100, false, 300, true, 'I', "RX:Z:AAA-CCC")
                .pair("q2", "rgA", 100, false, 300, true, '5', "RX:Z:" + umi);

        final Throwable thrown = runExpectingFailure(input.write(), NO_OPTICAL, "BARCODE_TAG=RX", "DUPLEX_UMI=" + duplex);

        Assert.assertTrue(thrown instanceof PicardException, thrown.toString());
        Assert.assertTrue(thrown.getMessage().contains(expectedMessage), thrown.getMessage());
        if (duplex && !expectedMessage.contains("illegal")) {
            Assert.assertTrue(thrown.getMessage().contains(umi), thrown.getMessage());
        }
    }
}
