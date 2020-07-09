/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

import org.testng.annotations.Test;

/**
 * The purpose of this class is to show that MarkDuplicates gives (almost) the same results when run on files that
 * are queryname sorted.
 *
 */
public class MarkDuplicatesTestQueryNameSorted extends MarkDuplicatesTest {
    final static String TEST_DATA_DIR = "testdata/picard/sam/MarkDuplicates";

    @Override
    protected boolean markUnmappedRecordsLikeTheirMates() {
        return true;
    }

    @Override
    protected boolean markSecondaryAndSupplementaryRecordsLikeTheCanonical() {
        return true;
    }

    @Override
    protected AbstractMarkDuplicatesCommandLineProgramTester getTester() {
        return new QuerySortedMarkDuplicatesTester();
    }

    @Test
    public void testOpticalDuplicateFinding() {
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();

        // explicitly creating 1 expected optical duplicate pair
        tester.setExpectedOpticalDuplicate(2);

        // pass in the read names manually, in order to control duplicates vs optical duplicates
        tester.addMatePair("READ0:1:1:1:1", 1, 1, 100, false, false, false, false, "50M", "50M", false, true, false,
                false, false, DEFAULT_BASE_QUALITY); // non-duplicate mapped pair to start
        tester.addMatePair("READ1:1:1:1:300", 1, 1, 100, false, false, true, true, "50M", "50M", false, true, false,
                false, false, DEFAULT_BASE_QUALITY); // duplicate pair, NOT optical duplicate (delta-Y > 100)
        tester.addMatePair("READ1:1:1:1:50", 1, 1, 100, false, false, true, true, "50M", "50M", false, true, false,
                false, false, DEFAULT_BASE_QUALITY); // two marked optical duplicates in sequence, because the sort order is queryname this might fail
        tester.addMatePair("READ2:1:1:2:50", 1, 1, 100, false, false, true, true, "50M", "50M", false, true, false,
                false, false, DEFAULT_BASE_QUALITY); // duplicate pair, expected optical duplicate (delta-X and delta-Y < 100)
        tester.runTest();
    }
}
