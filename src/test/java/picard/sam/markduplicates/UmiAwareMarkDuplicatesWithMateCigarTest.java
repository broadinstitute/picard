/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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
import java.util.*;
import picard.PicardException;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by UmiAwareMarkDuplicatesWithMateCigarTester (see getTester).
 * @author fleharty
 */
public class UmiAwareMarkDuplicatesWithMateCigarTest extends SimpleMarkDuplicatesWithMateCigarTest {

    protected UmiAwareMarkDuplicatesWithMateCigarTester getTester() {
        return new UmiAwareMarkDuplicatesWithMateCigarTester();
    }

    @DataProvider(name = "testUmiSetsDataProvider")
    private Object[][] testUmiSetsDataProvider() {
        return new Object[][]{
                {
                        // Test basic error correction using edit distance of 1
                        Arrays.asList(new String[] {"AAAA", "AAAA", "ATTA", "AAAA", "AAAT"}), // Observed UMI
                        Arrays.asList(new String[] {"AAAA", "AAAA", "ATTA", "AAAA", "AAAA"}), // Expected inferred UMI
                        Arrays.asList(new Boolean[] {false, true, false, true, true}), // Should it be marked as duplicate?
                        1 // Edit Distance to Join
                },
                {
                        // Test basic error correction using edit distance of 2
                        Arrays.asList(new String[] {"AAAA", "AAAA", "ATTA", "AAAA", "AAAT"}),
                        Arrays.asList(new String[] {"AAAA", "AAAA", "AAAA", "AAAA", "AAAA"}),
                        Arrays.asList(new Boolean[] {false, true, true, true, true}),
                        2
                },
                {
                        // Test basic error correction using edit distance of 1 where UMIs
                        // form a chain in edit distance space so that a UMI with large
                        // edit distance will get error corrected to a distant but linked (in edit space) UMI
                        Arrays.asList(new String[] {"AAAA", "AAAA", "AAAT", "AAGT", "ACGT", "TCGT", "CCCC"}),
                        Arrays.asList(new String[] {"AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "CCCC"}),
                        Arrays.asList(new Boolean[] {false, true, true, true, true, true, false}),
                        1
                },
                {
                        // Test short UMIs
                        Arrays.asList(new String[] {"A", "A", "T", "G", "G", "C", "C", "A"}),
                        Arrays.asList(new String[] {"A", "A", "A", "A", "A", "A", "A", "A"}), // All UMIs should get corrected to A
                        Arrays.asList(new Boolean[] {false, true, true, true, true, true, true, true}), // All mate pairs should be duplicates except the first
                        1
                },
                {
                        // Test short UMIs with no allowance for errors
                        Arrays.asList(new String[] {"A", "A", "T", "G", "G", "C", "C", "A"}),
                        Arrays.asList(new String[] {"A", "A", "T", "G", "G", "C", "C", "A"}), // No UMIs should get corrected
                        Arrays.asList(new Boolean[] {false, true, false, false, true, false, true, true}), // Only exactly duplicated UMIs will give rise to a new duplicate set
                        0
                },
                {
                        // Test longish UMIs with relatively large allowance for error
                        // UMIs "TTGACATCCA", "ATGCCATCGA", "AAGTCACCGT" should belong to the same duplicate set since
                        // they are within edit distance of 4 of each other.  TTGACATCCA should be chosen as the inferred
                        // UMI even though it only occurs once.  Since all UMIs only occur once, we choose the UMI that
                        // is not marked as duplicate to be the inferred UMI.
                        Arrays.asList(new String[] {"TTGACATCCA", "ATGCCATCGA", "AAGTCACCGT"}),
                        Arrays.asList(new String[] {"TTGACATCCA", "TTGACATCCA", "TTGACATCCA"}), // All UMIs should get corrected to TTGACATCCA
                        Arrays.asList(new Boolean[] {false, true, true}), // All mate pairs should be duplicates except the first
                        4
                },
                {
                        // Test to make sure that if any reads don't have a UMI, we treat things as if there were no
                        // UMIs at all
                        Arrays.asList(new String[] {"TTGACATCCA", null, null}),
                        Arrays.asList(new String[] {null, null, null}), // Since we had missing UMIs, no UMIs should be inferred
                        Arrays.asList(new Boolean[] {false, true, true}), // All mate pairs should be duplicates except the first
                        4
                }

        };
    }

    @DataProvider(name = "testBadUmiSetsDataProvider")
    private Object[][] testBadUmiSetsDataProvider() {
        return new Object[][]{
                {
                        // The code should not support variable length UMIs, if we observe variable length UMIs
                        // ensure that an exception is thrown.
                        Arrays.asList(new String[] {"AAAA", "A"}),
                        Arrays.asList(new String[] {"AAAA", "A"}),
                        Arrays.asList(new Boolean[] {false, false}),
                        4
                },
                {
                        // The code should not support variable length UMIs, if we observe variable length UMIs
                        // ensure that an exception is thrown.
                        Arrays.asList(new String[] {"T", "GG"}),
                        Arrays.asList(new String[] {"T", "GG"}),
                        Arrays.asList(new Boolean[] {false, false}),
                        1
                }
        };
    }

    @Test(dataProvider = "testUmiSetsDataProvider")
    public void testUmi(List<String> umis, List<String> inferredUmi, final List<Boolean> isDuplicate, final int editDistanceToJoin) {
        UmiAwareMarkDuplicatesWithMateCigarTester tester = getTester();
        tester.addArg("EDIT_DISTANCE_TO_JOIN=" + editDistanceToJoin);
        tester.addArg("ADD_INFERRED_UMI=TRUE");

        for(int i = 0;i < umis.size();i++) {
            tester.addMatePairWithUmi(umis.get(i), inferredUmi.get(i), isDuplicate.get(i), isDuplicate.get(i));
        }
        tester.setExpectedInferredUmis(inferredUmi);
        tester.runTest();
    }


    @Test(dataProvider = "testBadUmiSetsDataProvider", expectedExceptions = PicardException.class)
    public void testBadUmi(List<String> umis, List<String> inferredUmi, final List<Boolean> isDuplicate, final int editDistanceToJoin) {
        UmiAwareMarkDuplicatesWithMateCigarTester tester = getTester();
        tester.addArg("EDIT_DISTANCE_TO_JOIN=" + editDistanceToJoin);
        tester.addArg("ADD_INFERRED_UMI=TRUE");

        for(int i = 0;i < umis.size();i++) {
            tester.addMatePairWithUmi(umis.get(i), inferredUmi.get(i), isDuplicate.get(i), isDuplicate.get(i));
        }
        tester.setExpectedInferredUmis(inferredUmi);
        tester.runTest();
    }

}
