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

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;

/**
 * Tests {@link MarkDuplicates} on reads clipped at the start of a contig, whose unclipped 5' coordinate (the
 * alignment start minus the leading soft and hard clips) is zero or negative. Reads are duplicates when their
 * unclipped 5' coordinates are equal, wherever they align, and are not when those coordinates differ by one, so the
 * coordinates are neither clamped to the contig nor shifted. An unclipped read cannot start before position 1, so the
 * comparison with an unclipped read is made at position 1.
 *
 * Forward fragments on chr1 (unclipped 5' coordinate in brackets): m4a, m4b, m4c [-4]; m3a [-3]; z0a, z0b [0];
 * p1a (unclipped), p1b [1]. Pairs with the forward read on chr2 and the reverse mate on chr3 at 5' end 300, which wait
 * for their mates in the map of unmatched read ends: pa, pb [-4]; pc [-3]. In each set the read or pair with the
 * higher sum of base qualities is kept.
 */
public class MarkDuplicatesContigStartClippingTest {
    private static final int CHR1 = 0, CHR2 = 1, CHR3 = 2;

    private static MarkDuplicatesEndToEndRunner input() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        input.fragment("m4a", "rg", CHR1, 1, false, "5S10M", 40);
        input.fragment("m4b", "rg", CHR1, 3, false, "7S10M", 15);
        input.fragment("m4c", "rg", CHR1, 1, false, "2H3S10M", 15);
        input.fragment("m3a", "rg", CHR1, 1, false, "4S10M", 15);
        input.fragment("z0a", "rg", CHR1, 1, false, "1S10M", 40);
        input.fragment("z0b", "rg", CHR1, 2, false, "2S10M", 15);
        input.fragment("p1a", "rg", CHR1, 1, false, "10M", 15);
        input.fragment("p1b", "rg", CHR1, 4, false, "3S10M", 40);
        input.pair("pa", "rg", CHR2, 1, false, "5S10M", CHR3, 291, true, "10M", 40);
        input.pair("pb", "rg", CHR2, 2, false, "6S10M", CHR3, 291, true, "10M", 15);
        input.pair("pc", "rg", CHR2, 1, false, "4S10M", CHR3, 291, true, "10M", 15);
        return input;
    }

    /** A SORTING_COLLECTION_SIZE_RATIO that leaves room for two read ends in memory in this JVM. */
    private static String tinySortingCollectionRatio() {
        return Double.toString(2.5 * ReadEndsForMarkDuplicates.getSizeOf() / Runtime.getRuntime().maxMemory());
    }

    @DataProvider(name = "spilling")
    public Object[][] spilling() {
        return new Object[][]{
                {new String[0]},
                {new String[]{"SORTING_COLLECTION_SIZE_RATIO=" + tinySortingCollectionRatio(), "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1"}},
        };
    }

    @Test(dataProvider = "spilling")
    public void testFivePrimeCoordinatesAtOrBeforeTheContigStart(final String[] args) {
        final MarkDuplicatesEndToEndRunner input = input();
        final MarkDuplicatesEndToEndRunner.Result result = input.run(args);
        MarkDuplicatesEndToEndRunner.assertSameRecords(input.inputRecords(), result.records, "PG");
        result.assertDuplicates("m4b", "m4c", "z0b", "p1a", "pb/1", "pb/2");
        MarkDuplicatesEndToEndRunner.assertMetrics(result.onlyMetrics(), "lib", 8, 3, 0, 0, 4, 1, 0);
    }
}
