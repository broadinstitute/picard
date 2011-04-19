/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.illumina.parser;

import net.sf.samtools.util.StringUtil;
import net.sf.picard.util.SolexaQualityConverter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
public class IlluminaDataProviderTest {

    public static final File TEST_DATA_LOCATION = new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir");
    public static final String RUN_BARCODE = "305PJAAXX080716";
    private static final File RTA_BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/rta/BasecallsDir");


    /**
     * This overload uses the standard TEST_DATA_LOCATION Basecalls directory.
     */
    @Test(dataProvider="data")
    public void testNonRtaIlluminaDataProvider(
            final String testName, final int lane, final int size, final boolean pe, final String firstRead1Seq,
            final String firstRead1Qual, final String firstRead2Seq, final String firstRead2Qual,
            final short[] firstIntensities, final boolean firstPf, final String lastRead1Seq, final String lastRead1Qual,
            final String lastRead2Seq, final String lastRead2Qual, final short[] lastIntensities, final boolean lastPf,
            final boolean isBustard_1_1, final boolean parseProcessedIntensities, final boolean parseRawIntensities,
            final boolean parseNoise, final int barcodeCycle, final int barcodeLength,
            final String firstBarcodeSeq, final String firstBarcodeQual, final int seekAfterFirstRead)
            throws Exception {

        final IlluminaDataProviderFactory.BaseCallerVersion expectedBaseCallerVersion =
                isBustard_1_1? IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_1:
                IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_3;

        final IlluminaDataProviderFactory.ImageAnalyzerVersion expectedImageAnalyzerVersion =
                isBustard_1_1? IlluminaDataProviderFactory.ImageAnalyzerVersion.Firecrest_1_1:
                        IlluminaDataProviderFactory.ImageAnalyzerVersion.Firecrest_1_3;

        testIlluminaDataProvider(testName, lane, size, pe, firstRead1Seq, firstRead1Qual, firstRead2Seq, firstRead2Qual,
                firstIntensities, firstPf, lastRead1Seq, lastRead1Qual, lastRead2Seq, lastRead2Qual, lastIntensities,
                lastPf, expectedBaseCallerVersion, expectedImageAnalyzerVersion, parseProcessedIntensities,
                parseRawIntensities, parseNoise, barcodeCycle,
                barcodeLength, firstBarcodeSeq, firstBarcodeQual, seekAfterFirstRead, TEST_DATA_LOCATION);
    }

    @Test(dataProvider="arbitraryBasecallsDirData")
    public void testIlluminaDataProvider(
            final String testName, final int lane, final int size, final boolean pe, final String firstRead1Seq,
            final String firstRead1Qual, final String firstRead2Seq, final String firstRead2Qual,
            final short[] firstIntensities, final boolean firstPf, final String lastRead1Seq, final String lastRead1Qual,
            final String lastRead2Seq, final String lastRead2Qual, final short[] lastIntensities, final boolean lastPf,
            final IlluminaDataProviderFactory.BaseCallerVersion expectedBaseCallerVersion,
            final IlluminaDataProviderFactory.ImageAnalyzerVersion expectedImageAnalyzerVersion,
            final boolean parseProcessedIntensities,
            final boolean parseRawIntensities,
            final boolean parseNoise, final int barcodeCycle, final int barcodeLength,
            final String firstBarcodeSeq, final String firstBarcodeQual, final int seekAfterFirstRead, final File basecallsDirectory)
            throws Exception {

        final boolean isBustard_1_1 = expectedBaseCallerVersion == IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_1;
        // Convert expected qualities into Phred space
        final SolexaQualityConverter qualityConverter = SolexaQualityConverter.getSingleton();
        final byte[] firstRead1QualBytes = StringUtil.stringToBytes(firstRead1Qual);
        convertSolexaQualityCharsToPhredBinary(isBustard_1_1, qualityConverter, firstRead1QualBytes);
        byte[] firstRead2QualBytes = null;
        if (firstRead2Qual != null) {
            firstRead2QualBytes = StringUtil.stringToBytes(firstRead2Qual);
            convertSolexaQualityCharsToPhredBinary(isBustard_1_1, qualityConverter, firstRead2QualBytes);
        }
        final byte[] lastRead1QualBytes = StringUtil.stringToBytes(lastRead1Qual);
        convertSolexaQualityCharsToPhredBinary(isBustard_1_1, qualityConverter, lastRead1QualBytes);
        byte[] lastRead2QualBytes = null;
        if (lastRead2Qual != null) {
            lastRead2QualBytes = StringUtil.stringToBytes(lastRead2Qual);
            convertSolexaQualityCharsToPhredBinary(isBustard_1_1, qualityConverter, lastRead2QualBytes);
        }

        final IlluminaDataProviderFactory factory;
        if (barcodeLength == 0) {
            factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, IlluminaDataType.BaseCalls,
                IlluminaDataType.QualityScores, IlluminaDataType.PF);
        } else {
            factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, barcodeCycle, barcodeLength,
                    IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF);
        }
        if (parseProcessedIntensities){
            factory.addDataType(IlluminaDataType.ProcessedIntensities);
        }
        if (parseRawIntensities){
            factory.addDataType(IlluminaDataType.RawIntensities);
        }
        if (parseNoise){
            factory.addDataType(IlluminaDataType.Noise);
        }
        final AbstractIlluminaDataProvider parser = factory.makeDataProvider();
        Assert.assertEquals(factory.getBaseCallerVersion(), expectedBaseCallerVersion);
        Assert.assertEquals(factory.getImageAnalyzerVersion(), expectedImageAnalyzerVersion);
        int count = 0;
        while (parser.hasNext()) {
            final IlluminaReadData n = parser.next();
            if (count == 0) {
                Assert.assertEquals(StringUtil.bytesToString(n.getFirstEnd().getBases()), firstRead1Seq, testName);
                Assert.assertEquals(n.getFirstEnd().getQualities(), firstRead1QualBytes, testName);
                if (pe) {
                    Assert.assertEquals(StringUtil.bytesToString(n.getSecondEnd().getBases()), firstRead2Seq, testName);
                    Assert.assertEquals(n.getSecondEnd().getQualities(), firstRead2QualBytes, testName);
                }
                if (firstIntensities != null) {
                    final FourChannelIntensityData processedIntensities = n.getFirstEnd().getProcessedIntensities();
                    Assert.assertEquals(processedIntensities.getA()[0], firstIntensities[0], testName);
                    Assert.assertEquals(processedIntensities.getC()[0], firstIntensities[1], testName);
                    Assert.assertEquals(processedIntensities.getG()[0], firstIntensities[2], testName);
                    Assert.assertEquals(processedIntensities.getT()[0], firstIntensities[3], testName);
                }
                Assert.assertEquals(firstPf, n.isPf().booleanValue(), testName);
                if (firstBarcodeSeq != null) {
                    Assert.assertEquals(StringUtil.bytesToString(n.getBarcodeRead().getBases()), firstBarcodeSeq);
                }
                if (firstBarcodeQual != null) {
                    final byte[] firstBarcodeQualBytes = StringUtil.stringToBytes(firstBarcodeQual);
                    convertSolexaQualityCharsToPhredBinary(isBustard_1_1, qualityConverter, firstBarcodeQualBytes);
                    Assert.assertEquals(n.getBarcodeRead().getQualities(), firstBarcodeQualBytes);
                }
                if (seekAfterFirstRead != 0) {
                    parser.seekToTile(seekAfterFirstRead);
                }
            }
            else if (count == (size-1)) {
                Assert.assertEquals(StringUtil.bytesToString(n.getFirstEnd().getBases()), lastRead1Seq, testName);
                Assert.assertEquals(n.getFirstEnd().getQualities(), lastRead1QualBytes, testName);
                if (pe) {
                    Assert.assertEquals(StringUtil.bytesToString(n.getSecondEnd().getBases()), lastRead2Seq, testName);
                    Assert.assertEquals(n.getSecondEnd().getQualities(), lastRead2QualBytes, testName);
                }
                if (lastIntensities != null && pe) {
                    final FourChannelIntensityData processedIntensities = n.getSecondEnd().getProcessedIntensities();
                    final int lastIntensityIndex = processedIntensities.getA().length - 1;
                    Assert.assertEquals(processedIntensities.getA()[lastIntensityIndex], lastIntensities[0], testName);
                    Assert.assertEquals(processedIntensities.getC()[lastIntensityIndex], lastIntensities[1], testName);
                    Assert.assertEquals(processedIntensities.getG()[lastIntensityIndex], lastIntensities[2], testName);
                    Assert.assertEquals(processedIntensities.getT()[lastIntensityIndex], lastIntensities[3], testName);
                }
                Assert.assertEquals(lastPf, n.isPf().booleanValue(), testName);
            }
            count++;
        }
        Assert.assertEquals(count, size, testName);
    }

    private void convertSolexaQualityCharsToPhredBinary(final boolean isBustard_1_1,
                                                        final SolexaQualityConverter qualityConverter,
                                                        final byte[] qualBytes) {
        if (isBustard_1_1) {
            qualityConverter.convertSolexaQualityCharsToPhredBinary(qualBytes);
        } else {
            qualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(qualBytes);
        }
    }

    /**
     * Data for one PE lane and one non-PE lane.  The data are truncated versions of the qseq and sig2 files
     * for actual lanes.
     */
    @DataProvider(name = "data")
    private Object[][] getBasecallsTestData()
    {
        return new Object[][]{
                {"PE Bustard Parsing Test", 1, 60, true,
                        "G....................C.....................T.....................T..........",
                        "\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD",
                        "C...A................A.....................C.....................T..........",
                        "^DDDIDDDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD",
                        new short[] { -2, 71, 370, 7}, false,
                        "C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT",
                        "^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW",
                        "TCCATCCACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC",
                        "aaa[`aa_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG",
                        new short[] { 179, 211, 55, -35}, false, false, true, true, true, 0, 0, null, null, 0},

                {"PE with Barcode Bustard Parsing Test", 1, 60, true,
                        "G....................C.....................T.....................T..........",
                        "\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD",
                        "...............A.....................C.....................T..........",
                        "DDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD",
                        new short[] { -2, 71, 370, 7}, false,
                        "C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT",
                        "^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW",
                        "CACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC",
                        "a_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG",
                        new short[] { 179, 211, 55, -35}, false, false, true, false, false, 77, 6,
                        "C...A.", "^DDDID", 0},

                {"PE, Barcode, seek Bustard Parsing Test", 1, 21, true,
                        "G....................C.....................T.....................T..........",
                        "\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD",
                        "...............A.....................C.....................T..........",
                        "DDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD",
                        new short[] { -2, 71, 370, 7}, false,
                        "C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT",
                        "^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW",
                        "CACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC",
                        "a_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG",
                        new short[] { 179, 211, 55, -35}, false, false, true, false, false, 77, 6,
                        "C...A.", "^DDDID", 3},

                {"PE Bustard Parsing No Sig2 Files Test", 4, 60, true,
                        "G....................C.....................T.....................T..........",
                        "\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD",
                        "C...A................A.....................C.....................T..........",
                        "^DDDIDDDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD",
                        null, false,
                        "C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT",
                        "^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW",
                        "TCCATCCACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC",
                        "aaa[`aa_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG",
                        null, false, false, false, false, false, 0, 0, null, null, 0},

                {"Non-PE Bustard Parsing Test", 2, 20, false,
                        "GACTTTGGGAAGGGTCATTACTGCCCTTGTAGAAAGAACACCTCATGTTCCTTATCGAGAGCGGCCGCTGCTGATC",
                        "W[`bbbb_baS\\`_\\bbabbaWR`bba``ab_bbbbbabbbaabb^^^\\aa_ab`a_`[`VaST[^SWTXWNYEHM",
                        null, null, new short[] { 52, -44, 422, 175}, true,
                        "CACACACACACACACACACACACCACCTTTTGGCTTATCTGCACGCGGCCGCGTGCCCTACCCTACCCCATGGGAT",
                        "a_aa^\\Ra\\`aaXa_aa_aaaaaaa[^X^``V[`_`a^aaO``^_SJUELTVMKVTPFOKJJMNKTRRJEEPFKTR",
                        null, null, new short[] { 5, -18, 13, 73}, true, false, true, false, false, 0, 0, null, null, 0},


                {"Non-PE Bustard Parsing No Sig2 Files Test", 5, 20, false,
                        "GACTTTGGGAAGGGTCATTACTGCCCTTGTAGAAAGAACACCTCATGTTCCTTATCGAGAGCGGCCGCTGCTGATC",
                        "W[`bbbb_baS\\`_\\bbabbaWR`bba``ab_bbbbbabbbaabb^^^\\aa_ab`a_`[`VaST[^SWTXWNYEHM",
                        null, null, null, true,
                        "CACACACACACACACACACACACCACCTTTTGGCTTATCTGCACGCGGCCGCGTGCCCTACCCTACCCCATGGGAT",
                        "a_aa^\\Ra\\`aaXa_aa_aaaaaaa[^X^``V[`_`a^aaO``^_SJUELTVMKVTPFOKJJMNKTRRJEEPFKTR",
                        null, null, null, true, false, false, false, false, 0, 0, null, null, 0},


                {"Non-PE Barcode Bustard Parsing No Sig2 Files Test", 5, 20, false,
                        "GACTTTGGGAAGGGTCATTACTGCCCTTGTAGAAAGAACACCTCATGTTCCTTATCGAGAGCGGCCGCTG",
                        "W[`bbbb_baS\\`_\\bbabbaWR`bba``ab_bbbbbabbbaabb^^^\\aa_ab`a_`[`VaST[^SWTX",
                        null, null, null, true,
                        "CACACACACACACACACACACACCACCTTTTGGCTTATCTGCACGCGGCCGCGTGCCCTACCCTACCCCA",
                        "a_aa^\\Ra\\`aaXa_aa_aaaaaaa[^X^``V[`_`a^aaO``^_SJUELTVMKVTPFOKJJMNKTRRJE",
                        null, null, null, true, false, false, false, false, 71, 6, "CTGATC", "WNYEHM", 0},


                {"PE Bustard 1.1 Parsing Test", 3, 60, true,
                        "GGTGGGGTCGTATTTCTCGATGAAGGTGCCGGTCACGAACTGCACGGTCAGGGCGGAGTGGCACACACAGACAGAG",
                        "hhLhhhhNhhghhh`chchhUhhhhhNhTShhWOWSh_WMTgLWKgQELYWQ_FdPN???GB?@PLBE?BDB@?@D",
                        "GAGTACACAGTGGTGGGGCTGGGCTCGGGCGGGGGAGGCAAATCCCGCCTGAACGCCGAGCTCGTGACAGGCACCT",
                        "hhhghhhhhhLhhHhh?hXOhPhJFJYhhMhhhhDBWhDOPRBDB?@EK@B???DB???G?BDJ?E?B?DKB?B@@",
                        null, false,
                        "GTTTAAGTGACAGCTGGCTTACCCTTACCCCCAGGCCACCAGATGAAATTGTGGAGCAGAGCCTGCCTGGCACGTA",
                        "hhhghehhhhhhhh_hh`e`h]LRhXheYhZSZSeWPWQOJhPRhPJUOJYH`DVhJXZY\\DDHOHDFZEBLGOBD",
                        "TCCTTTGTCCCCTCTTGTCGGAGCACTCTACCTGCTCAGTGCCACCAGGACGAGCTCACGCATTCTGCACCGCCCT",
                        "`hchhhhh_ghhhhhXhVhhhBghXUhhhGSdf]fhYLeJhOQJ`PDNNFRB@JKHGDEMEDFJGGKJBGDQBDHE",
                        null, true, true, false, false, false, 0, 0, null, null, 0},

                {"PE Barcode Bustard 1.1 Parsing Test", 3, 60, true,
                        "GGTGGGGTCGTATTTCTCGATGAAGGTGCCGGTCACGAACTGCACGGTCAGGGCGGAGTGGCACACACAG",
                        "hhLhhhhNhhghhh`chchhUhhhhhNhTShhWOWSh_WMTgLWKgQELYWQ_FdPN???GB?@PLBE?B",
                        "GAGTACACAGTGGTGGGGCTGGGCTCGGGCGGGGGAGGCAAATCCCGCCTGAACGCCGAGCTCGTGACAGGCACCT",
                        "hhhghhhhhhLhhHhh?hXOhPhJFJYhhMhhhhDBWhDOPRBDB?@EK@B???DB???G?BDJ?E?B?DKB?B@@",
                        null, false,
                        "GTTTAAGTGACAGCTGGCTTACCCTTACCCCCAGGCCACCAGATGAAATTGTGGAGCAGAGCCTGCCTGG",
                        "hhhghehhhhhhhh_hh`e`h]LRhXheYhZSZSeWPWQOJhPRhPJUOJYH`DVhJXZY\\DDHOHDFZE",
                        "TCCTTTGTCCCCTCTTGTCGGAGCACTCTACCTGCTCAGTGCCACCAGGACGAGCTCACGCATTCTGCACCGCCCT",
                        "`hchhhhh_ghhhhhXhVhhhBghXUhhhGSdf]fhYLeJhOQJ`PDNNFRB@JKHGDEMEDFJGGKJBGDQBDHE",
                        null, true, true, false, false, false, 71, 6, "ACAGAG", "DB@?@D", 0},

                {"PE Barcode, seek Bustard 1.1 Parsing Test", 3, 21, true,
                        "GGTGGGGTCGTATTTCTCGATGAAGGTGCCGGTCACGAACTGCACGGTCAGGGCGGAGTGGCACACACAG",
                        "hhLhhhhNhhghhh`chchhUhhhhhNhTShhWOWSh_WMTgLWKgQELYWQ_FdPN???GB?@PLBE?B",
                        "GAGTACACAGTGGTGGGGCTGGGCTCGGGCGGGGGAGGCAAATCCCGCCTGAACGCCGAGCTCGTGACAGGCACCT",
                        "hhhghhhhhhLhhHhh?hXOhPhJFJYhhMhhhhDBWhDOPRBDB?@EK@B???DB???G?BDJ?E?B?DKB?B@@",
                        null, false,
                        "GTTTAAGTGACAGCTGGCTTACCCTTACCCCCAGGCCACCAGATGAAATTGTGGAGCAGAGCCTGCCTGG",
                        "hhhghehhhhhhhh_hh`e`h]LRhXheYhZSZSeWPWQOJhPRhPJUOJYH`DVhJXZY\\DDHOHDFZE",
                        "TCCTTTGTCCCCTCTTGTCGGAGCACTCTACCTGCTCAGTGCCACCAGGACGAGCTCACGCATTCTGCACCGCCCT",
                        "`hchhhhh_ghhhhhXhVhhhBghXUhhhGSdf]fhYLeJhOQJ`PDNNFRB@JKHGDEMEDFJGGKJBGDQBDHE",
                        null, true, true, false, false, false, 71, 6, "ACAGAG", "DB@?@D", 3},

                {"Non-PE Bustard 1.1 Parsing Test", 6, 60, false,
                        "GGTGGGGTCGTATTTCTCGATGAAGGTGCCGGTCACGAACTGCACGGTCAGGGCGGAGTGGCACACACAGACAGAGGAGTACACAGTGGTGGGGCTGGGCTCGGGCGGGGGAGGCAAATCCCGCCTGAACGCCGAGCTCGTGACAGGCACCT",
                        "hhLhhhhNhhghhh`chchhUhhhhhNhTShhWOWSh_WMTgLWKgQELYWQ_FdPN???GB?@PLBE?BDB@?@DhhhghhhhhhLhhHhh?hXOhPhJFJYhhMhhhhDBWhDOPRBDB?@EK@B???DB???G?BDJ?E?B?DKB?B@@",
                        null,
                        null,
                        null, false,
                        "GTTTAAGTGACAGCTGGCTTACCCTTACCCCCAGGCCACCAGATGAAATTGTGGAGCAGAGCCTGCCTGGCACGTATCCTTTGTCCCCTCTTGTCGGAGCACTCTACCTGCTCAGTGCCACCAGGACGAGCTCACGCATTCTGCACCGCCCT",
                        "hhhghehhhhhhhh_hh`e`h]LRhXheYhZSZSeWPWQOJhPRhPJUOJYH`DVhJXZY\\DDHOHDFZEBLGOBD`hchhhhh_ghhhhhXhVhhhBghXUhhhGSdf]fhYLeJhOQJ`PDNNFRB@JKHGDEMEDFJGGKJBGDQBDHE",
                        null,
                        null,
                        null, true, true, false, false, false, 0, 0, null, null, 0},

                {"Non-PE Barcode Bustard 1.1 Parsing Test", 6, 60, false,
                        "GTCGTATTTCTCGATGAAGGTGCCGGTCACGAACTGCACGGTCAGGGCGGAGTGGCACACACAGACAGAGGAGTACACAGTGGTGGGGCTGGGCTCGGGCGGGGGAGGCAAATCCCGCCTGAACGCCGAGCTCGTGACAGGCACCT",
                        "hNhhghhh`chchhUhhhhhNhTShhWOWSh_WMTgLWKgQELYWQ_FdPN???GB?@PLBE?BDB@?@DhhhghhhhhhLhhHhh?hXOhPhJFJYhhMhhhhDBWhDOPRBDB?@EK@B???DB???G?BDJ?E?B?DKB?B@@",
                        null,
                        null,
                        null, false,
                        "GTGACAGCTGGCTTACCCTTACCCCCAGGCCACCAGATGAAATTGTGGAGCAGAGCCTGCCTGGCACGTATCCTTTGTCCCCTCTTGTCGGAGCACTCTACCTGCTCAGTGCCACCAGGACGAGCTCACGCATTCTGCACCGCCCT",
                        "hhhhhhhh_hh`e`h]LRhXheYhZSZSeWPWQOJhPRhPJUOJYH`DVhJXZY\\DDHOHDFZEBLGOBD`hchhhhh_ghhhhhXhVhhhBghXUhhhGSdf]fhYLeJhOQJ`PDNNFRB@JKHGDEMEDFJGGKJBGDQBDHE",
                        null,
                        null,
                        null, true, true, false, false, false, 1, 6, "GGTGGG", "hhLhhh", 0},

        };
    }

    @DataProvider(name = "arbitraryBasecallsDirData")
    private Object[][] getRtaBasecallsTestData()
    {
        return new Object[][]{
                {"Barcode-aware PE Bustard Parsing Test", 7, 10, true,
                        "TAGAGATGGC.CT.........T......C........G...TCCAGACCGCCCATTCTCTGCCTGCC",
                        "^`^BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                        "CCTCTAATCCCAGCACTATCCGAGACCAAATCAGGCAAATCACTTGAAGTCAGGAGTTCGAGACCAGC",
                        "]]VISQK_M\\`MHX\\ZFMaPWHXYUa]ZHJGaULGPXTRS\\W[`ZH_GMUa_[]M]PTUZX]VaZaU\\",
                        null, false,
                        "CAACTCTTGT.GT........GT......A........G...AATATATTCTGAAACTCAGCAATGTT",
                        "aaabaaBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                        "TAACTTTCAGAGGCCCTTCAGGAGGCCCTGGCCTGTCAAGTACCTTTACAGTGATGGGTATAGACTTT",
                        "abbbabba_bbabb_]S]ab_ab__^abSORYX^RF[aa`_a_[XVa`[WUN`a^__a\\ZL\\BBBBBB",
                        null, false,
                        IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_5,
                        IlluminaDataProviderFactory.ImageAnalyzerVersion.rta,
                        false, false, false, 69, 8,
                        "CCACCCAC", "_ZFZ^]BB", 0, RTA_BASECALLS_DIR},
        };
    }


    @Test
    public void testBarcodeParsing() {
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(TEST_DATA_LOCATION, 6, 1, 6,
                    IlluminaDataType.Barcodes);
        factory.computeReadConfiguration();
        Assert.assertEquals(factory.getBaseCallerVersion(), IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_1);
        Assert.assertEquals(factory.getImageAnalyzerVersion(), IlluminaDataProviderFactory.ImageAnalyzerVersion.Firecrest_1_1);
        final AbstractIlluminaDataProvider parser = factory.makeDataProvider();
        for (int i = 0; parser.hasNext(); ++i) {
            final IlluminaReadData read = parser.next();
            final String matchedBarcode = read.getMatchedBarcode();
            if (i % 2 == 0) {
                // The barcode are not actual sequence in the test data, just 0-padded numbers
                Assert.assertNotNull(matchedBarcode);
                final int barcodeAsInt = Integer.parseInt(matchedBarcode);
                Assert.assertEquals(barcodeAsInt, i+1);
            } else {
                Assert.assertNull(matchedBarcode);
            }
        }
    }

    @Test
    public void testIntensitiesAvailable() {
        // Check for intensities that are present
        IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(RTA_BASECALLS_DIR, 1,
                    IlluminaDataType.RawIntensities);
        factory.setImageAnalyzerVersion(IlluminaDataProviderFactory.ImageAnalyzerVersion.rta);
        factory.setBaseCallerVersion(IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_4);
        factory.setPairedEnd(false);
        Assert.assertTrue(factory.intensitiesAvailable(1));

        // Check for intensities that are not present
        factory = new IlluminaDataProviderFactory(RTA_BASECALLS_DIR, 2,
                    IlluminaDataType.RawIntensities);
        factory.setImageAnalyzerVersion(IlluminaDataProviderFactory.ImageAnalyzerVersion.rta);
        factory.setBaseCallerVersion(IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_4);
        factory.setPairedEnd(false);
        Assert.assertFalse(factory.intensitiesAvailable(1));


    }
}
