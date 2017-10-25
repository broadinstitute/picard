/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import htsjdk.samtools.util.QualityUtil;
import picard.PicardException;

import java.util.*;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by UmiAwareMarkDuplicatesWithMateCigarTester (see getTester).
 * @author fleharty
 */
public class UmiAwareMarkDuplicatesWithMateCigarTest extends SimpleMarkDuplicatesWithMateCigarTest {

    @Override
    protected UmiAwareMarkDuplicatesWithMateCigarTester getTester() {
        return new UmiAwareMarkDuplicatesWithMateCigarTester();
    }

    protected UmiAwareMarkDuplicatesWithMateCigarTester getTester(final boolean allowMissingUmis) {
        return new UmiAwareMarkDuplicatesWithMateCigarTester(allowMissingUmis);
    }

    @DataProvider(name = "testUmiSetsDataProvider")
    private Object[][] testUmiSetsDataProvider() {
        return new Object[][] {{
                // Test basic error correction using edit distance of 1
                Arrays.asList(new String[] {"AAAA", "AAAA", "ATTA", "AAAA", "AAAT"}), // Observed UMI
                Arrays.asList(new String[] {"AAAA", "AAAA", "ATTA", "AAAA", "AAAA"}), // Expected inferred UMI
                Arrays.asList(new Boolean[] {false, true, false, true, true}), // Should it be marked as duplicate?
                1 // Edit Distance to Join
        }, {
                // Test basic error correction using edit distance of 1 including dashes
                Arrays.asList(new String[] {"AA-AA", "--AAAA", "A-T-TA", "-AAAA----", "A-AA-T"}), // Observed UMI
                Arrays.asList(new String[] {"AAAA", "AAAA", "ATTA", "AAAA", "AAAA"}), // Expected inferred UMI
                Arrays.asList(new Boolean[] {false, true, false, true, true}), // Should it be marked as duplicate?
                1 // Edit Distance to Join
        }, {
                // Test basic error correction using edit distance of 2
                Arrays.asList(new String[] {"AAAA", "AAAA", "ATTA", "AAAA", "AAAT"}),
                Arrays.asList(new String[] {"AAAA", "AAAA", "AAAA", "AAAA", "AAAA"}),
                Arrays.asList(new Boolean[] {false, true, true, true, true}),
                2
        }, {
                // Test basic error correction using edit distance of 2 including dashes
                Arrays.asList(new String[] {"AAA-A", "A--AAA", "A---TT-A", "----AAAA", "A-AAT"}),
                Arrays.asList(new String[] {"AAAA", "AAAA", "AAAA", "AAAA", "AAAA"}),
                Arrays.asList(new Boolean[] {false, true, true, true, true}),
                2
        }, {
                // Test basic error correction using edit distance of 1 where UMIs
                // form a chain in edit distance space so that a UMI with large
                // edit distance will get error corrected to a distant but linked (in edit space) UMI
                Arrays.asList(new String[] {"AAAA", "AAAA", "AAAT", "AAGT", "ACGT", "TCGT", "CCCC"}),
                Arrays.asList(new String[] {"AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "CCCC"}),
                Arrays.asList(new Boolean[] {false, true, true, true, true, true, false}),
                1
        }, {
                // Test short UMIs
                Arrays.asList(new String[] {"A", "A", "T", "G", "G", "C", "C", "A"}),
                Arrays.asList(new String[] {"A", "A", "A", "A", "A", "A", "A", "A"}), // All UMIs should get corrected to A
                Arrays.asList(new Boolean[] {false, true, true, true, true, true, true, true}), // All mate pairs should be duplicates except the first
                1
        }, {
                // Test short UMIs with no allowance for errors
                Arrays.asList(new String[] {"A", "A", "T", "G", "G", "C", "C", "A"}),
                Arrays.asList(new String[] {"A", "A", "T", "G", "G", "C", "C", "A"}), // No UMIs should get corrected
                Arrays.asList(new Boolean[] {false, true, false, false, true, false, true, true}), // Only exactly duplicated UMIs will give rise to a new duplicate set
                0
        }, {
                // Test longish UMIs with relatively large allowance for error
                // UMIs "TTGACATCCA", "ATGCCATCGA", "AAGTCACCGT" should belong to the same duplicate set since
                // they are within edit distance of 4 of each other.  TTGACATCCA should be chosen as the inferred
                // UMI even though it only occurs once.  Since all UMIs only occur once, we choose the UMI that
                // is not marked as duplicate to be the inferred UMI.
                Arrays.asList(new String[] {"TTGACATCCA", "ATGCCATCGA", "AAGTCACCGT"}),
                Arrays.asList(new String[] {"TTGACATCCA", "TTGACATCCA", "TTGACATCCA"}), // All UMIs should get corrected to TTGACATCCA
                Arrays.asList(new Boolean[] {false, true, true}), // All mate pairs should be duplicates except the first
                4
        }, {
                // Test that the inferred UMI is correct with N
                Arrays.asList(new String[] {"AAAN", "AANN"}),
                Arrays.asList(new String[] {"AAAN", "AAAN"}), // Only N containing UMI
                Arrays.asList(new Boolean[] {false, true}), // All mate pairs should be duplicates except the first
                1
        }, {
                // Test that the majority with no Ns wins
                Arrays.asList(new String[] {"AAAN", "AAAN", "AAAA", "AAAA", "AAAN" }),
                Arrays.asList(new String[] {"AAAA", "AAAA", "AAAA", "AAAA", "AAAA"}), // Even though AAAN is majority, AAAA is represented
                Arrays.asList(new Boolean[] {false, true, true, true, true}), // All mate pairs should be duplicates except the first
                1
        }, {
                // Test that the majority with the fewest N wins when both have Ns
                Arrays.asList(new String[] {"AAAN", "AAAN", "AANN", "AANN", "AANN" }),
                Arrays.asList(new String[] {"AAAN", "AAAN", "AAAN", "AAAN", "AAAN"}), // Even though AANN is majority, AAAN is represented
                Arrays.asList(new Boolean[] {false, true, true, true, true, true}), // All mate pairs should be duplicates except the first
                1
        } };
    }

    @DataProvider(name = "testBadUmiSetsDataProvider")
    private Object[][] testBadUmiSetsDataProvider() {
        return new Object[][] {{
                // The code should not support variable length UMIs, if we observe variable length UMIs
                // ensure that an exception is thrown.
                Arrays.asList(new String[] {"AAAA", "A"}),
                Arrays.asList(new String[] {"AAAA", "A"}),
                Arrays.asList(new Boolean[] {false, false}),
                4
        }, {
                // The code should not support variable length UMIs, if we observe variable length UMIs
                // ensure that an exception is thrown.
                // Arrays.asList(new String[] {"T", "GG"}),
                Arrays.asList(new String[] {"T", "GG"}),
                Arrays.asList(new String[] {"T", "GG"}),
                Arrays.asList(new Boolean[] {false, false}),
                1
        }, {
                // Test to make sure that we throw an exception with missing UMIs when allowMissingUmis is false
                // This throws an exception because the UMIs have differing lengths.
                Arrays.asList(new String[] {"TTGA", "TTAT", null}),
                Arrays.asList(new String[] {"TTGA", "TTAT", null}),
                Arrays.asList(new Boolean[] {false, false, false}),
                4
        }};
    }

    @DataProvider(name = "testEmptyUmiDataProvider")
    private Object[][] testEmptyUmiDataProvider() {
        return new Object[][] {{
                // Test to make sure we treat empty UMIs correctly when they are allowed
                Arrays.asList(new String[] {null, null, null}),
                Arrays.asList(new String[] {null, null, null}),
                Arrays.asList(new Boolean[] {false, true, true}),
                4
        }};
    }

    @Test(dataProvider = "testUmiSetsDataProvider")
    public void testUmi(List<String> umis, List<String> assignedUmi, final List<Boolean> isDuplicate, final int editDistanceToJoin) {
        UmiAwareMarkDuplicatesWithMateCigarTester tester = getTester(false);
        tester.addArg("MAX_EDIT_DISTANCE_TO_JOIN=" + editDistanceToJoin);

        for(int i = 0;i < umis.size();i++) {
            tester.addMatePairWithUmi(umis.get(i), assignedUmi.get(i), isDuplicate.get(i), isDuplicate.get(i));
        }
        tester.setExpectedAssignedUmis(assignedUmi).runTest();
    }

    @Test(dataProvider = "testEmptyUmiDataProvider")
    public void testEmptyUmis(List<String> umis, List<String> assignedUmi, final List<Boolean> isDuplicate, final int editDistanceToJoin) {
        UmiAwareMarkDuplicatesWithMateCigarTester tester = getTester(true);
        tester.addArg("MAX_EDIT_DISTANCE_TO_JOIN=" + editDistanceToJoin);

        for(int i = 0;i < umis.size();i++) {
            tester.addMatePairWithUmi(umis.get(i), assignedUmi.get(i), isDuplicate.get(i), isDuplicate.get(i));
        }
        tester.setExpectedAssignedUmis(assignedUmi).runTest();
    }

    @Test(dataProvider = "testBadUmiSetsDataProvider", expectedExceptions = {IllegalArgumentException.class, PicardException.class})
    public void testBadUmis(List<String> umis, List<String> assignedUmi, final List<Boolean> isDuplicate, final int editDistanceToJoin) {
        UmiAwareMarkDuplicatesWithMateCigarTester tester = getTester(false);
        tester.addArg("MAX_EDIT_DISTANCE_TO_JOIN=" + editDistanceToJoin);

        for(int i = 0;i < umis.size();i++) {
            tester.addMatePairWithUmi(umis.get(i), assignedUmi.get(i), isDuplicate.get(i), isDuplicate.get(i));
        }
        tester.setExpectedAssignedUmis(assignedUmi).runTest();
    }

    @DataProvider(name = "testUmiMetricsDataProvider")
    private Object[][] testUmiMetricsDataProvider() {

        // Calculate values of metrics by hand to ensure they are right
        // effectiveLength4_1 is the effective UMI length observing 5 UMIs where 4 are the same
        double effectiveLength4_1 = -(4./5.)*Math.log(4./5.)/Math.log(4.) -(1./5.)*Math.log(1./5.)/Math.log(4.);
        // effectiveLength4_1 is the effective UMI length observing 5 UMIs where 3 are the same and the other two are
        // unique
        double effectiveLength3_1_1 = -(3./5.)*Math.log(3./5.)/Math.log(4.) -2*(1./5.)*Math.log(1./5.)/Math.log(4.);

        double effectiveLength_N = -(3./4.)*Math.log(3./4.)/Math.log(4.) -(1./4.)*Math.log(1./4.)/Math.log(4.);

        // estimatedBaseQualityk_n is the phred scaled base quality score where k of n bases are incorrect
        double estimatedBaseQuality1_20 = QualityUtil.getPhredScoreFromErrorProbability(1./20.);
        double estimatedBaseQuality3_20 = QualityUtil.getPhredScoreFromErrorProbability(3./20.);
        double estimatedBaseQuality_N = QualityUtil.getPhredScoreFromErrorProbability(2./16.);

        double estimatedPercentWithN3_7 = 3./7.;

        return new Object[][]{{
                // Test basic error correction using edit distance of 1
                Arrays.asList(new String[]{"AAAA", "AAAA", "ATTA", "AAAA", "AAAT"}), // Observed UMI
                Arrays.asList(new String[]{"AAAA", "AAAA", "ATTA", "AAAA", "AAAA"}), // Expected inferred UMI
                Arrays.asList(new Boolean[]{false, true, false, true, true}), // Should it be marked as duplicate?
                1, // Edit Distance to Join
                new UmiMetrics(4.0,                         // MEAN_UMI_LENGTH
                               3,                           // OBSERVED_UNIQUE_UMIS
                               2,                           // INFERRED_UNIQUE_UMIS
                               2,                           // OBSERVED_BASE_ERRORS (Note: This is 2 rather than 1 because we are using paired end reads)
                               2,                           // DUPLICATE_SETS_WITHOUT_UMI
                               4,                           // DUPLICATE_SETS_WITH_UMI
                                effectiveLength4_1,         // EFFECTIVE_LENGTH_OF_INFERRED_UMIS
                                effectiveLength3_1_1,       // EFFECTIVE_LENGTH_OF_OBSERVED_UMIS
                                estimatedBaseQuality1_20,   // ESTIMATED_BASE_QUALITY_OF_UMIS
                               0)                           // UMI_WITH_N
        }, {
                // Test basic error correction using edit distance of 2
                Arrays.asList(new String[]{"AAAA", "AAAA", "ATTA", "AAAA", "AAAT"}),
                Arrays.asList(new String[]{"AAAA", "AAAA", "AAAA", "AAAA", "AAAA"}),
                Arrays.asList(new Boolean[]{false, true, true, true, true}),
                2,
                new UmiMetrics(4.0,                         // MEAN_UMI_LENGTH
                               3,                           // OBSERVED_UNIQUE_UMIS
                               1,                           // INFERRED_UNIQUE_UMIS
                               6,                           // OBSERVED_BASE_ERRORS
                               2,                           // DUPLICATE_SETS_WITHOUT_UMI
                               2,                           // DUPLICATE_SETS_WITH_UMI
                               0.0,                         // EFFECTIVE_LENGTH_OF_INFERRED_UMIS
                                effectiveLength3_1_1,       // EFFECTIVE_LENGTH_OF_OBSERVED_UMIS
                                estimatedBaseQuality3_20,   // ESTIMATED_BASE_QUALITY_OF_UMIS
                               0)                           // UMI_WITH_N
        }, {
                // Test basic error correction using edit distance of 2 including dashes
                Arrays.asList(new String[]{"AAAA", "AAA-A", "A-TTA", "--AAA-A", "AAA-T"}),
                Arrays.asList(new String[]{"AAAA", "AAAA", "AAAA", "AAAA", "AAAA"}),
                Arrays.asList(new Boolean[]{false, true, true, true, true}),
                2,
                new UmiMetrics(4.0,                     // MEAN_UMI_LENGTH
                        3,                              // OBSERVED_UNIQUE_UMIS
                        1,                              // INFERRED_UNIQUE_UMIS
                        6,                              // OBSERVED_BASE_ERRORS
                        2,                              // DUPLICATE_SETS_WITHOUT_UMI
                        2,                              // DUPLICATE_SETS_WITH_UMI
                        0.0,                            // EFFECTIVE_LENGTH_OF_INFERRED_UMIS
                         effectiveLength3_1_1,          // EFFECTIVE_LENGTH_OF_OBSERVED_UMIS
                         estimatedBaseQuality3_20,      // ESTIMATED_BASE_QUALITY_OF_UMIS
                        0)                              // UMI_WITH_N
        },{
                // Test basic error correction using edit distance of 2 - Ns metrics should not include the umis with Ns
                Arrays.asList(new String[]{"AAAA", "AAAA","AANA","ANNA", "ATTA", "AAAA", "ANAT"}),
                Arrays.asList(new String[]{"AAAA", "AAAA","AAAA","AAAA","AAAA", "AAAA", "AAAA"}),
                Arrays.asList(new Boolean[]{false, true, true, true, true, true, true}),
                2,
                new UmiMetrics(4.0,                 // MEAN_UMI_LENGTH
                        2,                          // OBSERVED_UNIQUE_UMIS
                        1,                          // INFERRED_UNIQUE_UMIS
                        4,                          // OBSERVED_BASE_ERRORS
                        2,                          // DUPLICATE_SETS_WITHOUT_UMI
                        2,                          // DUPLICATE_SETS_WITH_UMI
                        0.0,                        // EFFECTIVE_LENGTH_OF_INFERRED_UMIS
                         effectiveLength_N,         // EFFECTIVE_LENGTH_OF_OBSERVED_UMIS
                         estimatedBaseQuality_N,    // ESTIMATED_BASE_QUALITY_OF_UMIS
                        estimatedPercentWithN3_7)   // UMI_WITH_N
        }, {
                // Test basic error correction using edit distance of 2 including Ns and dashes
                Arrays.asList(new String[]{"AAAA-", "AA-AA","AAN-A","ANNA", "AT-TA", "AAA-A-", "A--NAT-"}),
                Arrays.asList(new String[]{"AAAA", "AAAA","AAAA","AAAA","AAAA", "AAAA", "AAAA"}),
                Arrays.asList(new Boolean[]{false, true, true, true, true, true, true}),
                2,
                new UmiMetrics(4.0,                 // MEAN_UMI_LENGTH
                        2,                          // OBSERVED_UNIQUE_UMIS
                        1,                          // INFERRED_UNIQUE_UMIS
                        4,                          // OBSERVED_BASE_ERRORS
                        2,                          // DUPLICATE_SETS_WITHOUT_UMI
                        2,                          // DUPLICATE_SETS_WITH_UMI
                        0.0,                        // EFFECTIVE_LENGTH_OF_INFERRED_UMIS
                        effectiveLength_N,         // EFFECTIVE_LENGTH_OF_OBSERVED_UMIS
                        estimatedBaseQuality_N,    // ESTIMATED_BASE_QUALITY_OF_UMIS
                        estimatedPercentWithN3_7)   // UMI_WITH_N
        }, {
                // Test maximum entropy (EFFECTIVE_LENGTH_OF_INFERRED_UMIS)
                Arrays.asList(new String[]{"AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG", "CA", "CT", "CC", "CG", "GA", "GT", "GC", "GG"}),
                Arrays.asList(new String[]{"AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG", "CA", "CT", "CC", "CG", "GA", "GT", "GC", "GG"}),
                Arrays.asList(new Boolean[]{false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false}),
                0,
                new UmiMetrics(2.0,    // MEAN_UMI_LENGTH
                               16,     // OBSERVED_UNIQUE_UMIS
                               16,     // INFERRED_UNIQUE_UMIS
                               0,      // OBSERVED_BASE_ERRORS
                               2,      // DUPLICATE_SETS_WITHOUT_UMI
                               32,     // DUPLICATE_SETS_WITH_UMI
                               2.0,    // EFFECTIVE_LENGTH_OF_INFERRED_UMIS
                               2,      // EFFECTIVE_LENGTH_OF_OBSERVED_UMIS
                               -1,     // ESTIMATED_BASE_QUALITY_OF_UMIS
                               0)      // UMI_WITH_N
        }};
    }

    @Test(dataProvider = "testUmiMetricsDataProvider")
    public void testUmiMetrics(List<String> umis, List<String> assignedUmi, final List<Boolean> isDuplicate,
                               final int editDistanceToJoin, final UmiMetrics expectedMetrics) {
        UmiAwareMarkDuplicatesWithMateCigarTester tester = getTester(false);
        tester.addArg("MAX_EDIT_DISTANCE_TO_JOIN=" + editDistanceToJoin);

        for( int i = 0;i < umis.size();i++ ) {
            tester.addMatePairWithUmi(umis.get(i), assignedUmi.get(i), isDuplicate.get(i), isDuplicate.get(i));
        }
        tester.setExpectedAssignedUmis(assignedUmi);
        tester.setExpectedMetrics(expectedMetrics);
        tester.runTest();
    }

    @DataProvider(name = "testUmiUtilDataProvider")
    private Object[][] testUmiUtilDataProvider() {
        return new Object[][]{{
                Arrays.asList(new String[]{"AAAA", "AA-AA", "-A-T-A", "AAAAA--", "---A", "---", ""}), // Observed UMI
                Arrays.asList(new String[]{"AAAA", "AAAA", "ATA", "AAAAA", "A", "", ""}) // Sanitized UMI
        }};
    }

    @Test(dataProvider = "testUmiUtilDataProvider")
    public void testUmiUtil(List<String> observed, List<String> expected) {
        for( int i = 0;i < observed.size();i++ ) {
            SAMRecord rec = new SAMRecord(new SAMFileHeader());
            rec.setAttribute("RX", observed.get(i));
            Assert.assertEquals(UmiUtil.getSanitizedUMI(rec, "RX"), expected.get(i));
        }
    }
}
