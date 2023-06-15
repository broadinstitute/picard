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

import java.util.Arrays;

/**
 * This class defines the individual test cases to check the duplicate set histograms.
 *
 * @author hogstrom@broadinstitute.org
 */
public class MarkDuplicatesSetSizeHistogramTest extends AbstractMarkDuplicatesCommandLineProgramTest {
    @Override
    protected MarkDuplicatesTester getTester() {
        return new MarkDuplicatesTester();
    }

    @Test
    public void TestSingleSet() {
        // This tests checks that if I add two read pairs with the same alignment start, a duplicate
        // set size entry of 2 is put into the histogram. The pair is an optical dup, so a 2 entry for
        // the 'optical_sets' should also be marked
        final MarkDuplicatesTester tester = getTester();
        tester.setExpectedOpticalDuplicate(1);
        String representativeReadName = "RUNID:1:1:16020:13352";
        tester.addMatePair(representativeReadName, 1, 485253, 485253, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:1:1:15993:13361", 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        // set expected counts in hashmap that takes the form: key=(Histogram Label, histogram bin), value=histogram entry
        tester.expectedSetSizeMap.put(Arrays.asList("all_sets", "2.0"), 1.0);
        tester.expectedSetSizeMap.put(Arrays.asList("optical_sets", "2.0"), 1.0);
        tester.runTest();
    }

    @Test
    public void testOpticalAndNonOpticalSet() {
        // This tests checks that if I have two optical dups and one non-optical in the same dup set that the correct histogram entries
        // are specified
        final MarkDuplicatesTester tester = getTester();
        tester.setExpectedOpticalDuplicate(2);
        String representativeReadName = "RUNID:1:1:16020:13352";
        tester.addMatePair(representativeReadName, 1, 485253, 485253, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:1:1:15993:13361", 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:1:1:15994:13364", 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        // add one non-optical duplicate
        tester.addMatePair("RUNID:1:1:25993:23361", 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedSetSizeMap.put(Arrays.asList("all_sets", "4.0"), 1.0);
        tester.expectedSetSizeMap.put(Arrays.asList("optical_sets", "3.0"), 1.0);
        tester.runTest();
    }

    @Test
    public void testSingleton() {
        // This tests checks if having a read pair that is not duplicated, the correct histogram entry is updated
        final MarkDuplicatesTester tester = getTester();
        tester.setExpectedOpticalDuplicate(0);
        // Add non-duplicate read
        tester.addMatePair("RUNID:1:1:15993:13375", 1, 485255, 485255, false, false, false, false, "43M2S", "43M2S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedSetSizeMap.put(Arrays.asList("all_sets", "1.0"), 1.0);
        tester.runTest();
    }

    @Test
    public void testMultiRepresentativeReadTags() {
        // This tests checks multiple different duplicate sets of varying sizes
        final MarkDuplicatesTester tester = getTester();
        tester.setExpectedOpticalDuplicate(3);
        // Duplicate set: size 2 - all optical
        String representativeReadName1 = "RUNID:1:1:16020:13352";
        tester.addMatePair(representativeReadName1, 1, 485253, 485253, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:1:1:15993:13361", 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedSetSizeMap.put(Arrays.asList("all_sets", "3.0"), 1.0);
        tester.expectedSetSizeMap.put(Arrays.asList("optical_sets", "3.0"), 1.0);

        // Duplicate set: size 3 - all optical
        String representativeReadName2 = "RUNID:1:1:15993:13360";
        tester.addMatePair(representativeReadName2, 1, 485299, 485299, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:1:1:15993:13365", 1, 485299, 485299, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:1:1:15993:13370", 1, 485299, 485299, false, false, true, true, "43M2S", "43M2S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedSetSizeMap.put(Arrays.asList("all_sets", "4.0"), 1.0);
        tester.expectedSetSizeMap.put(Arrays.asList("optical_sets", "4.0"), 1.0);

        // Add non-duplicate read
        tester.addMatePair("RUNID:1:1:15993:13375", 1, 485255, 485255, false, false, false, false, "43M2S", "43M2S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedSetSizeMap.put(Arrays.asList("all_sets", "1.0"), 1.0);

        tester.runTest();
    }

}
