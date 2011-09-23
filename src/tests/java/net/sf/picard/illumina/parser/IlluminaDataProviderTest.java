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
import java.util.Arrays;

/**
 * @author alecw@broadinstitute.org
 */
public class IlluminaDataProviderTest {

    public static final File PARSING_TEST_BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaBarcodeParsingTest/Basecalls");
    public static final File TEST_DATA_LOCATION = new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir");
    public static final String RUN_BARCODE = "305PJAAXX080716";
    private static final File RTA_BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/rta/BasecallsDir");
    private static final IlluminaDataType [] DEFAULT_DATA_TYPES = new IlluminaDataType[]{IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF};


    /**
     * This overload uses the standard TEST_DATA_LOCATION Basecalls directory.
     */
    @Test(dataProvider="data")
    public void testNonRtaIlluminaDataProvider(
            final String testName, final int lane, final int size, final boolean pe,
            final IlluminaReadData firstRead,
            final IlluminaReadData lastRead,
            final IlluminaDataType [] extraDataTypes,
            final int barcodeCycle, final int barcodeLength,
            final int seekAfterFirstRead)
            throws Exception {

        final IlluminaDataProviderFactory.BaseCallerVersion expectedBaseCallerVersion = IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_3;

        testIlluminaDataProvider(testName, lane, size, pe,
                firstRead, lastRead,
                expectedBaseCallerVersion,
                extraDataTypes,
                barcodeCycle, barcodeLength,
                seekAfterFirstRead,
                TEST_DATA_LOCATION);
    }

    @Test(dataProvider="arbitraryBasecallsDirData")
    public void testIlluminaDataProvider(
            final String testName, final int lane, final int size, final boolean pe,
            final IlluminaReadData firstRead,
            final IlluminaReadData lastRead,
            final IlluminaDataProviderFactory.BaseCallerVersion expectedBaseCallerVersion,
            final IlluminaDataType [] extraDataTypes,
            final int barcodeCycle, final int barcodeLength,
            final int seekAfterFirstRead,
            final File basecallsDirectory)
            throws Exception {

        final IlluminaDataProviderFactory factory;
        if (barcodeLength == 0) {
            factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, DEFAULT_DATA_TYPES);
        } else {
            factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, barcodeCycle, barcodeLength, DEFAULT_DATA_TYPES);
        }

        for(final IlluminaDataType dt : extraDataTypes) {
            factory.addDataType(dt);
        }

        final AbstractIlluminaDataProvider parser = factory.makeDataProvider();
        Assert.assertEquals(factory.getBaseCallerVersion(), expectedBaseCallerVersion);

        int count = 0;
        while (parser.hasNext()) {
            final IlluminaReadData n = parser.next();
            if (count == 0) {
                compareReadData(pe, n, firstRead, testName);
                if (seekAfterFirstRead != 0) {
                    parser.seekToTile(seekAfterFirstRead);
                }
            }
            else if (count == (size-1)) {
                compareReadData(pe, n, lastRead, testName);
            }
            count++;
        }
        Assert.assertEquals(count, size, testName);
    }

    private void compareBasesAndQuals(final IlluminaEndData ied1, final IlluminaEndData ied2, final String testName) {
        Assert.assertEquals(ied1.getBases(),     ied2.getBases(),     testName);
        Assert.assertEquals(ied1.getQualities(), ied2.getQualities(), testName);
    }

    private void compareReadData(final boolean pe, final IlluminaReadData ird1, final IlluminaReadData ird2, final String testName) {
        compareBasesAndQuals(ird1.getFirstEnd(), ird2.getFirstEnd(), testName);

        if(pe) {
            Assert.assertTrue(ird1.isPairedEnd());
            Assert.assertTrue(ird2.isPairedEnd());
            compareBasesAndQuals(ird1.getSecondEnd(), ird2.getSecondEnd(), testName);
        }

        Assert.assertEquals(ird1.getBarcodeRead() == null, ird2.getBarcodeRead() == null, testName);
        if(ird1.getBarcodeRead() != null) {
            compareBasesAndQuals(ird1.getBarcodeRead(), ird2.getBarcodeRead(), testName);
        }

        Assert.assertEquals(ird1.isPf().booleanValue(),    ird2.isPf().booleanValue(), testName);
    }

    private byte[] solexaQualityCharsToPhredBinary(final byte[] qualBytes) {
        final byte[] outQualBytes = Arrays.copyOf(qualBytes, qualBytes.length);
        SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(outQualBytes);
        return outQualBytes;
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
                makeReadData(false, //first read
                    makeEndData(toBytes("G....................C.....................T.....................T.........."), //first end of first read
                                toBytes("\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD")),
                    makeEndData(toBytes("C...A................A.....................C.....................T.........."),
                                toBytes("^DDDIDDDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD")),
                    null),
                makeReadData(false,
                    makeEndData(toBytes("C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT"),
                                toBytes("^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW")),
                    makeEndData(toBytes("TCCATCCACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC"),
                                toBytes("aaa[`aa_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG")),
                    null),
                new IlluminaDataType[]{},
                0, 0, 0
            },

            {"PE with Barcode Bustard Parsing Test", 1, 60, true,
                makeReadData(false, //first read
                    makeEndData(toBytes("G....................C.....................T.....................T.........."),
                                toBytes("\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD")),
                    makeEndData(toBytes("...............A.....................C.....................T.........."),
                                toBytes("DDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD")),
                    makeEndData(toBytes("C...A."),
                                toBytes("^DDDID"))),
                makeReadData(false,
                    makeEndData(toBytes("C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT"),
                                toBytes("^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW")),
                    makeEndData(toBytes("CACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC"),
                                toBytes("a_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG")),
                    makeEndData(toBytes("TCCATC"),
                                toBytes("aaa[`a"))),
                new IlluminaDataType[]{},
                77, 6, 0
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20, true,
                makeReadData(false, //first read
                    makeEndData(toBytes("G..."), toBytes("\\DDD")),
                    makeEndData(toBytes("C..."), toBytes("^DDD")),
                    null),
                makeReadData(false,
                    makeEndData(toBytes("A..."), toBytes("GDDD")),
                    makeEndData(toBytes("A..."), toBytes("RDDD")),
                    null),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                0, 0, 0
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20, true,
                makeReadData(false, //first read
                    makeEndData(toBytes("G..."), toBytes("\\DDD")),
                    makeEndData(toBytes(".."),   toBytes("DD")),
                    makeEndData(toBytes("C."),   toBytes("^D"))),
                makeReadData(false,
                    makeEndData(toBytes("A..."), toBytes("GDDD")),
                    makeEndData(toBytes(".."),   toBytes("DD")),
                    makeEndData(toBytes("A."),   toBytes("RD"))),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                5, 2, 0
            },

            {"PE, Barcode, seek Bustard Parsing Test", 1, 21, true,
                makeReadData(false, //first read
                    makeEndData(toBytes("G....................C.....................T.....................T.........."),
                                toBytes("\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD")),
                    makeEndData(toBytes("...............A.....................C.....................T.........."),
                                toBytes("DDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD")),
                    makeEndData(toBytes("C...A."),
                                toBytes("^DDDID"))),
                makeReadData(false,
                    makeEndData(toBytes("C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT"),
                                toBytes("^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW")),
                    makeEndData(toBytes("CACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC"),
                                toBytes("a_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG")),
                    makeEndData(toBytes("TCCATC"),
                                toBytes("aaa[`a"))),
                new IlluminaDataType[]{},
                77, 6, 3
            },

            {"Non-PE Bustard Parsing Test", 5, 20, false,
                makeReadData(true,
                    makeEndData(toBytes("GACTTTGGGAAGGGTCATTACTGCCCTTGTAGAAAGAACACCTCATGTTCCTTATCGAGAGCGGCCGCTGCTGATC"),
                                toBytes("W[`bbbb_baS\\`_\\bbabbaWR`bba``ab_bbbbbabbbaabb^^^\\aa_ab`a_`[`VaST[^SWTXWNYEHM")),
                    null, null),
                makeReadData(true,
                    makeEndData(toBytes("CACACACACACACACACACACACCACCTTTTGGCTTATCTGCACGCGGCCGCGTGCCCTACCCTACCCCATGGGAT"),
                                toBytes("a_aa^\\Ra\\`aaXa_aa_aaaaaaa[^X^``V[`_`a^aaO``^_SJUELTVMKVTPFOKJJMNKTRRJEEPFKTR")),
                    null, null),
                new IlluminaDataType[]{},
                0, 0, 0
            },


            {"Non-PE Barcode Bustard Parsing Test", 5, 20, false,
                makeReadData(true,
                    makeEndData(toBytes("GACTTTGGGAAGGGTCATTACTGCCCTTGTAGAAAGAACACCTCATGTTCCTTATCGAGAGCGGCCGCTG"),
                                toBytes("W[`bbbb_baS\\`_\\bbabbaWR`bba``ab_bbbbbabbbaabb^^^\\aa_ab`a_`[`VaST[^SWTX")),
                    null,
                    makeEndData(toBytes("CTGATC"),
                                toBytes("WNYEHM"))),
                makeReadData(true,
                    makeEndData(toBytes("CACACACACACACACACACACACCACCTTTTGGCTTATCTGCACGCGGCCGCGTGCCCTACCCTACCCCA"),
                                toBytes("a_aa^\\Ra\\`aaXa_aa_aaaaaaa[^X^``V[`_`a^aaO``^_SJUELTVMKVTPFOKJJMNKTRRJE")),
                    null,
                    makeEndData(toBytes("TGGGAT"),
                                toBytes("EPFKTR"))),
                new IlluminaDataType[]{},
                71, 6, 0
            }
        };
    }

    @DataProvider(name = "arbitraryBasecallsDirData")
    private Object[][] getRtaBasecallsTestData()
    {
        return new Object[][]{
            {"Barcode-aware PE Bustard Parsing Test", 7, 10, true,
                makeReadData(false,
                    makeEndData(toBytes("TAGAGATGGC.CT.........T......C........G...TCCAGACCGCCCATTCTCTGCCTGCC"),
                                toBytes("^`^BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")),
                    makeEndData(toBytes("CCTCTAATCCCAGCACTATCCGAGACCAAATCAGGCAAATCACTTGAAGTCAGGAGTTCGAGACCAGC"),
                                toBytes("]]VISQK_M\\`MHX\\ZFMaPWHXYUa]ZHJGaULGPXTRS\\W[`ZH_GMUa_[]M]PTUZX]VaZaU\\")),
                    makeEndData(toBytes("CCACCCAC"),
                                toBytes("_ZFZ^]BB"))),
                makeReadData(false,
                    makeEndData(toBytes("CAACTCTTGT.GT........GT......A........G...AATATATTCTGAAACTCAGCAATGTT"),
                                toBytes("aaabaaBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")),
                    makeEndData(toBytes("TAACTTTCAGAGGCCCTTCAGGAGGCCCTGGCCTGTCAAGTACCTTTACAGTGATGGGTATAGACTTT"),
                                toBytes("abbbabba_bbabb_]S]ab_ab__^abSORYX^RF[aa`_a_[XVa`[WUN`a^__a\\ZL\\BBBBBB")),
                    makeEndData(toBytes("TCGGAATG"),
                                toBytes("abY`bb_^"))),
                    IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_5,
                    new IlluminaDataType[]{},
                    69, 8, 0, RTA_BASECALLS_DIR
            },
        };
    }

    @Test
    public void testBarcodeParsing() {
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(PARSING_TEST_BASECALLS_DIR, 6, 153, 6, IlluminaDataType.Barcodes);
        factory.computeReadConfiguration();

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
        factory.setBaseCallerVersion(IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_4);
        factory.setPairedEnd(false);
        Assert.assertTrue(factory.intensitiesAvailable(1));

        // Check for intensities that are not present
        factory = new IlluminaDataProviderFactory(RTA_BASECALLS_DIR, 2,
                    IlluminaDataType.RawIntensities);
        factory.setBaseCallerVersion(IlluminaDataProviderFactory.BaseCallerVersion.Bustard_1_4);
        factory.setPairedEnd(false);
        Assert.assertFalse(factory.intensitiesAvailable(1));
    }

    private byte [] toBytes(final String str) {
        return StringUtil.stringToBytes(str);
    }

    private IlluminaEndData makeEndData(final byte [] bases, final byte [] qualities) {
        final IlluminaEndData ied = new IlluminaEndData();
        ied.setBases(bases);
        ied.setQualities(solexaQualityCharsToPhredBinary(qualities));
        return ied;
    }

    private IlluminaEndData makeEndData(final byte [] bases, final byte [] qualities, final FourChannelIntensityData rawIntensities, final FourChannelIntensityData noise) {
        final IlluminaEndData ied = makeEndData(bases, qualities);
        ied.setRawIntensities(rawIntensities);
        ied.setNoise(noise);
        return ied;
    }

    private IlluminaReadData makeReadData(final boolean pf, final IlluminaEndData firstEnd, final IlluminaEndData secondEnd, final IlluminaEndData barcode) {
        final IlluminaReadData ird = new IlluminaReadData();
        ird.setPf(pf);
        ird.setFirstEnd(firstEnd);
        ird.setSecondEnd(secondEnd);
        ird.setBarcodeRead(barcode);

        return ird;
    }
}
