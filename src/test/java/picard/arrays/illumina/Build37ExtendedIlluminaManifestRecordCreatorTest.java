package picard.arrays.illumina;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class Build37ExtendedIlluminaManifestRecordCreatorTest {

    @DataProvider(name = "illuminaIndelDataProvider")
    public Object[][] illuminaIndelDataProvider() {
        return new Object[][]{
                {  // 1:100316615-CAG-C
                  "TCTTTTAATTCTTATGAAAGATTTCAAATCCTCTAGAAGCCAAAATGGGACACAGTAAAC", "AG", "ATTCGAATTTTACTTCTGAACGAAATGGAGAAACTGGAAAAGACCCTCTTCAGACTTGA", true,
                  "TCTTTTAATTCTTATGAAAGATTTCAAATCCTCTAGAAGCCAAAATGGGACACAGTAAAC",       "ATTCGAATTTTACTTCTGAACGAAATGGAGAAACTGGAAAAGACCCTCTTCAGACTTGA", 1.0, 59,
                   "CTTTTAATTCTTATGAAAGATTTCAAATCCTCTAGAAGCCAAAATGGGACACAGTAAACA",     "GATTCGAATTTTACTTCTGAACGAAATGGAGAAACTGGAAAAGACCCTCTTCAGACTTG", 0.0, 0
                },
                {  // 1:107600181-T-TC
                  "GCTTGGAGCAGGAGCTGGAGGCCGGAGTGGGCGGGCGCTTCCGCTGCAGCTGCTATGGCT", "C", "CGGCGCCCATGCATGGCTTTGCCATCTGGTTCCAGGTGACCTTCCCTGGAGGGGAGTCG", false,
                 "GGCTTGGAGCAGGAGCTGGAGGCCGGAGTGGGCGGGCGCTTCCGCTGCAGCTGCTATGG",        "CCGGCGCCCATGCATGGCTTTGCCATCTGGTTCCAGGTGACCTTCCCTGGAGGGGAGTCG", 0.0, 0,
                  "GCTTGGAGCAGGAGCTGGAGGCCGGAGTGGGCGGGCGCTTCCGCTGCAGCTGCTATGGCT",      "CGGCGCCCATGCATGGCTTTGCCATCTGGTTCCAGGTGACCTTCCCTGGAGGGGAGTCG", 1.0, 59
                }

        };
    }

    @Test(dataProvider = "illuminaIndelDataProvider")
    public void testCalculateIsDeletion(final String fivePrimeSeq, final String indelSeq, final String threePrimeSeq, final boolean indelSeqMatch,
                                        final String genomicDeletionSeqFivePrime, final String genomicDeletionSeqThreePrime,
                                        final Double expectedDeletionContextScore, final Integer expectedDeletionContextLength,
                                        final String unused1, final String unused2, final Double unused3, final Integer unused4) {
        ImmutablePair<Double, Integer> deletionContextInfo = Build37ExtendedIlluminaManifestRecordCreator.calculateIsDeletion(fivePrimeSeq, threePrimeSeq, indelSeqMatch,
                genomicDeletionSeqFivePrime, indelSeq, genomicDeletionSeqThreePrime);
        Assert.assertEquals(deletionContextInfo.left, expectedDeletionContextScore, 0.000001);
        Assert.assertEquals(deletionContextInfo.right, expectedDeletionContextLength);
    }

    @Test(dataProvider = "illuminaIndelDataProvider")
    public void testCalculateIsInsertion(final String fivePrimeSeq, final String indelSeq, final String threePrimeSeq, final boolean unused1,
                                         final String unused2, final String unused3,
                                         final Double unused4, final Integer unused5,
                                         final String genomicInsertionSeqFivePrime, final String genomicInsertionSeqThreePrime,
                                         final Double expectedInsertionContextScore, final Integer expectedInsertionContextLength
                                         ) {
        ImmutablePair<Double, Integer> insertionContextInfo = Build37ExtendedIlluminaManifestRecordCreator.calculateIsInsertion(fivePrimeSeq,
                threePrimeSeq, genomicInsertionSeqFivePrime,
                indelSeq, genomicInsertionSeqThreePrime);
        Assert.assertEquals(insertionContextInfo.left, expectedInsertionContextScore, 0.000001);
        Assert.assertEquals(insertionContextInfo.right, expectedInsertionContextLength);
    }

    @DataProvider(name = "illuminaLeftShiftDataProvider")
    public Object[][] illuminaLeftShiftDataProvider() {
        return new Object[][]{
                {"TCTTTTAATTCTTATGAAAGATTTCAAATCCTCTAGAAGCCAAAATGGGACACAGTAAAC", "AG", "ATTCGAATTTTACTTCTGAACGAAATGGAGAAACTGGAAAAGACCCTCTTCAGACTTGA",
                 "TCTTTTAATTCTTATGAAAGATTTCAAATCCTCTAGAAGCCAAAATGGGACACAGTAAAC",       "ATTCGAATTTTACTTCTGAACGAAATGGAGAAACTGGAAAAGACCCTCTTCAGACTTGA"},
                {"CTTTTAATTCTTATGAAAGATTTCAAATCCTCTAGAAGCCAAAATGGGACACAGTAAACA", "AG", "GATTCGAATTTTACTTCTGAACGAAATGGAGAAACTGGAAAAGACCCTCTTCAGACTTG",
                 "CTTTTAATTCTTATGAAAGATTTCAAATCCTCTAGAAGCCAAAATGGGACACAGTAAACA",       "GATTCGAATTTTACTTCTGAACGAAATGGAGAAACTGGAAAAGACCCTCTTCAGACTTG"},
                {"TGTTTTGAATGATTAAAACTACCATGTCTTATGTCATTTTTCAGGAAAGGCTATAAAGGT", "T", "CTCATATGATGAGTGGAACAGAAAAATACAAGACAACTTTGAAAAGCTATTTCATGTTT",
                 "TGTTTTGAATGATTAAAACTACCATGTCTTATGTCATTTTTCAGGAAAGGCTATAAAGG",      "TCTCATATGATGAGTGGAACAGAAAAATACAAGACAACTTTGAAAAGCTATTTCATGTTT"},
                {"TTCCAGATTGATGGGCCCGGAGACTACTGCAAAGACTATAGTTTTGGTTAAAAATGTTCT", "T", "TTCCCGACATTATGTTCATCTTGAGAGGTAAGTCATCAGGAGCATGTAATTTCCATAAC",
                 "TTCCAGATTGATGGGCCCGGAGACTACTGCAAAGACTATAGTTTTGGTTAAAAATGTTC",      "TTTCCCGACATTATGTTCATCTTGAGAGGTAAGTCATCAGGAGCATGTAATTTCCATAAC"},
                {"CGCCCATGCATGGCTTTGCCATCTGGTTCCAGGTGACCTTCCCTGGAGGGGAGTCGGAGA", "AA", "AACCCCTGGTGCTGTCCACCTCGCCTTTTCACCCGGCCACTCACTGGAAACAGGCGCTC",
                 "CGCCCATGCATGGCTTTGCCATCTGGTTCCAGGTGACCTTCCCTGGAGGGGAGTCGGAG",       "AAACCCCTGGTGCTGTCCACCTCGCCTTTTCACCCGGCCACTCACTGGAAACAGGCGCTC"},
                {"TACCAAACTTGGATTCTTTCCTGACCCTAGACCTTTTCCTCTGCCCTTATCATCGCTTTT", "T", "CAGTGATGGAGGAAATGTTGGTTGTGTTGATGTAATTATTCAAAGAGCATACCCTATACA",
                 "TACCAAACTTGGATTCTTTCCTGACCCTAGACCTTTTCCTCTGCCCTTATCATCGC",      "TTTTCAGTGATGGAGGAAATGTTGGTTGTGTTGATGTAATTATTCAAAGAGCATACCCTATACA"},
                {"TATCAAGGGATGTCACAACCGTGTGGAAGTTGCGTATTGTAAGCTATTCAAAAAAAGAAA", "AA", "GATTCAGGTAAGTATGTAAATGCTTTGTTTTTATCAGTTTTATTAACTTAAAAAATGACC",
                 "TATCAAGGGATGTCACAACCGTGTGGAAGTTGCGTATTGTAAGCTATTCAAAAAAAG",       "AAAGATTCAGGTAAGTATGTAAATGCTTTGTTTTTATCAGTTTTATTAACTTAAAAAATGACC"},
                {"GGACCAATAAGTCTTAATTGGTTTGAAGAACTTTCTTCAGAAGCTCCACCCTATAATTCT", "ATAATTCT", "GAACCTGCAGAAGAATCTGAACATAAAAACAACAATTACGAACCAAACCTATTTAAAACT",
                 "GGACCAATAAGTCTTAATTGGTTTGAAGAACTTTCTTCAGAAGCTCCACCCT",             "ATAATTCTGAACCTGCAGAAGAATCTGAACATAAAAACAACAATTACGAACCAAACCTATTTAAAACT"},
                {"GATAGACAGGTCGGCCACGGGAAAGTAGGTGAGGCCGCCAGGCGGGTCGGGGGCGGGGCT", "GGGCT", "GGGCTGGAAGGTGAGCTCGGGAACGTTGGTAGGGATGACGCGGTTGACAGCCAGAAAATG",
                 "GATAGACAGGTCGGCCACGGGAAAGTAGGTGAGGCCGCCAGGCGGGTCGGGGGCG",          "GGGCTGGGCTGGAAGGTGAGCTCGGGAACGTTGGTAGGGATGACGCGGTTGACAGCCAGAAAATG"},
                {"AAATTTTTAAAGGCACAAGAGGCCCTAGATTTCTATGGGGAAGTAAGGACCAGAGACAAA", "AA", "GGTAAGTTATTTTTTGATGTTTTTCCTTTCCTCTTCCTGGATCTGAGAATTTATTGGAAA",
                 "AAATTTTTAAAGGCACAAGAGGCCCTAGATTTCTATGGGGAAGTAAGGACCAGAGAC",       "AAAGGTAAGTTATTTTTTGATGTTTTTCCTTTCCTCTTCCTGGATCTGAGAATTTATTGGAAA"},
                {"TTATGGAAGATGATGAACTGACAGATTCTAAACTGCCAAGTCATGCCACACATTCNCTTT", "TT", "ACATGTCCCGAAAATGAGGAAATGGTTTTGTCAAATTCAAGAATTGGAAAAAGAAGAGGA",
                 "TTATGGAAGATGATGAACTGACAGATTCTAAACTGCCAAGTCATGCCACACATTCNC",       "TTTACATGTCCCGAAAATGAGGAAATGGTTTTGTCAAATTCAAGAATTGGAAAAAGAAGAGGA"},
                {"TTCTCTGTTACTTCTCCTAATGATTCTTCTTTTTTTTTCTTTTTTCTTTTTTTTTTTTTT", "TTT", "GCATACAGCAGAGAATGCCTGTGGGTTTAAAAATTTTTTAACTAATCTTGATGAATCAAT",
                 "TTCTCTGTTACTTCTCCTAATGATTCTTCTTTTTTTTTCTTTTTTC",        "TTTTTTTTTTTTTTGCATACAGCAGAGAATGCCTGTGGGTTTAAAAATTTTTTAACTAATCTTGATGAATCAAT"}
        };
    }

    @Test(dataProvider = "illuminaLeftShiftDataProvider")
    public void testIlluminaLeftShift(final String fivePrimeSeq, final String indel, final String threePrimeSeq,
                                      final String expectedLeftShiftedFivePrimeSeq,
                                      final String expectedRightShiftedThreePrimeSeq) {
        ImmutablePair<String, String> leftShiftedSeqs = Build37ExtendedIlluminaManifestRecordCreator.illuminaLeftShift(fivePrimeSeq, indel, threePrimeSeq);
        Assert.assertEquals(leftShiftedSeqs.left, expectedLeftShiftedFivePrimeSeq);
        Assert.assertEquals(leftShiftedSeqs.right, expectedRightShiftedThreePrimeSeq);
    }

    @DataProvider(name = "commonPrefixAndSuffixLengthDataProvider")
    public Object[][] commonPrefixAndSuffixLengthDataProvider() {
        return new Object[][] {
                {"", "A", 0, 0},
                {"A", "", 0, 0},
                {"A", "A", 1, 1},
                {"AA", "A", 1, 1},
                {"A", "AA", 1, 1},
                {"AC", "AA", 1, 0},
                {"AA", "AC", 1, 0},
                {"AAA", "ACA", 1, 1},
                {"ACA", "AAA", 1, 1},
                {"AAA", "CAA", 0, 2},
                {"CAA", "AAA", 0, 2},
                {"AAA", "AAC", 2, 0},
                {"AAC", "AAA", 2, 0},
                {"ACGGTTG", "CGGTTG", 0, 6},
                {"CGGTTG", "CGGTTGA", 6, 0}
        };
    }

    @Test(dataProvider = "commonPrefixAndSuffixLengthDataProvider")
    public void testCommonPrefixLength(final String seq1, final String seq2, final int expectedLength, final int ignoreMe) {
        final int commonLength = Build37ExtendedIlluminaManifestRecordCreator.commonPrefixLength(seq1, seq2);
        Assert.assertEquals(commonLength, expectedLength);
    }

    @Test(dataProvider = "commonPrefixAndSuffixLengthDataProvider")
    public void testCommonSuffixLength(final String seq1, final String seq2, final int ignoreMe, final int expectedLength) {
        final int commonLength = Build37ExtendedIlluminaManifestRecordCreator.commonSuffixLength(seq1, seq2);
        Assert.assertEquals(commonLength, expectedLength);
    }
}

