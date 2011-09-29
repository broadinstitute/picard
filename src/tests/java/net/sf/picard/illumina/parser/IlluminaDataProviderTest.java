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

    public static final File PARSING_TEST_BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaBarcodeParsingTest/BaseCalls");
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
            final ClusterData firstRead,
            final ClusterData lastRead,
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
            final ClusterData firstCluster,
            final ClusterData lastCluster,
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

        final IlluminaDataProvider parser = factory.makeDataProvider();
        Assert.assertEquals(factory.getBaseCallerVersion(), expectedBaseCallerVersion);

        int count = 0;
        while (parser.hasNext()) {
            final ClusterData cluster = parser.next();
            if (count == 0) {
                compareReadData(pe, cluster, firstCluster, testName);
                if (seekAfterFirstRead != 0) {
                    parser.seekToTile(seekAfterFirstRead);
                }
            }
            else if (count == (size-1)) {
                compareReadData(pe, cluster, lastCluster, testName);
            }
            count++;
        }
        Assert.assertEquals(count, size, testName);
    }

    private void compareBasesAndQuals(final ReadData rd1, final ReadData rd2, final String testName) {
        Assert.assertEquals(rd1.getBases(),     rd2.getBases(),     testName);
        Assert.assertEquals(rd1.getQualities(), rd2.getQualities(), testName);
    }

    private void compareReadData(final boolean pe, final ClusterData cd1, final ClusterData cd2, final String testName) {
        compareBasesAndQuals(cd1.getFirstEnd(), cd2.getFirstEnd(), testName);

        if(pe) {
            Assert.assertTrue(cd1.isPairedEnd());
            Assert.assertTrue(cd2.isPairedEnd());
            compareBasesAndQuals(cd1.getSecondEnd(), cd2.getSecondEnd(), testName);
        }

        Assert.assertEquals(cd1.getBarcodeRead() == null, cd2.getBarcodeRead() == null, testName);
        if(cd1.getBarcodeRead() != null) {
            compareBasesAndQuals(cd1.getBarcodeRead(), cd2.getBarcodeRead(), testName);
        }

        Assert.assertEquals(cd1.isPf().booleanValue(),    cd2.isPf().booleanValue(), testName);
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
                makeClusterData(false, //first read
                        makeReadData(toBytes("G....................C.....................T.....................T.........."), //first end of first read
                                toBytes("\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD")),
                        makeReadData(toBytes("C...A................A.....................C.....................T.........."),
                                toBytes("^DDDIDDDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD")),
                        null),
                makeClusterData(false,
                        makeReadData(toBytes("C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT"),
                                toBytes("^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW")),
                        makeReadData(toBytes("TCCATCCACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC"),
                                toBytes("aaa[`aa_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG")),
                        null),
                new IlluminaDataType[]{},
                0, 0, 0
            },

            {"PE with Barcode Bustard Parsing Test", 1, 60, true,
                makeClusterData(false, //first read
                        makeReadData(toBytes("G....................C.....................T.....................T.........."),
                                toBytes("\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD")),
                        makeReadData(toBytes("...............A.....................C.....................T.........."),
                                toBytes("DDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD")),
                        makeReadData(toBytes("C...A."),
                                toBytes("^DDDID"))),
                makeClusterData(false,
                        makeReadData(toBytes("C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT"),
                                toBytes("^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW")),
                        makeReadData(toBytes("CACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC"),
                                toBytes("a_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG")),
                        makeReadData(toBytes("TCCATC"),
                                toBytes("aaa[`a"))),
                new IlluminaDataType[]{},
                77, 6, 0
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20, true,
                makeClusterData(false, //first read
                        makeReadData(toBytes("G..."), toBytes("\\DDD")),
                        makeReadData(toBytes("C..."), toBytes("^DDD")),
                        null),
                makeClusterData(false,
                        makeReadData(toBytes("A..."), toBytes("GDDD")),
                        makeReadData(toBytes("A..."), toBytes("RDDD")),
                        null),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                0, 0, 0
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20, true,
                makeClusterData(false, //first read
                        makeReadData(toBytes("G..."), toBytes("\\DDD")),
                        makeReadData(toBytes(".."), toBytes("DD")),
                        makeReadData(toBytes("C."), toBytes("^D"))),
                makeClusterData(false,
                        makeReadData(toBytes("A..."), toBytes("GDDD")),
                        makeReadData(toBytes(".."), toBytes("DD")),
                        makeReadData(toBytes("A."), toBytes("RD"))),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                5, 2, 0
            },

            {"PE, Barcode, seek Bustard Parsing Test", 1, 21, true,
                makeClusterData(false, //first read
                        makeReadData(toBytes("G....................C.....................T.....................T.........."),
                                toBytes("\\DDDDDDDDDDDDDDDDDDDD\\DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDD")),
                        makeReadData(toBytes("...............A.....................C.....................T.........."),
                                toBytes("DDDDDDDDDDDDDDDKDDDDDDDDDDDDDDDDDDDDD[DDDDDDDDDDDDDDDDDDDDDMDDDDDDDDDD")),
                        makeReadData(toBytes("C...A."),
                                toBytes("^DDDID"))),
                makeClusterData(false,
                        makeReadData(toBytes("C.TAT.......C..CGGGTACCACAGTTGAGGACTGACATTCTGAACCCTGATGTTTCTAAAGAAACGACAGTAT"),
                                toBytes("^DU_WDDDDDDD^DD^_^_U```^^U[]]_UNV^`^^U][[W\\_QTQZ]_WS[X]TW_^VMLUZVWZ[SFXL[YUW")),
                        makeReadData(toBytes("CACTTCCCTGAGCCTCAGAAAAGGGCAAGGCATGGCTCACATACTCTCAGCCACGGCCTGGCCTGCTGCC"),
                                toBytes("a_aaaaaa\\__`aa^aT_VVV\\ZZ`a`X`Za^\\][aa_U_``^a]aX]I]``X`TR^]GDWXGMX]Z[YG")),
                        makeReadData(toBytes("TCCATC"),
                                toBytes("aaa[`a"))),
                new IlluminaDataType[]{},
                77, 6, 3
            },

            {"Non-PE Bustard Parsing Test", 5, 20, false,
                makeClusterData(true,
                        makeReadData(toBytes("GACTTTGGGAAGGGTCATTACTGCCCTTGTAGAAAGAACACCTCATGTTCCTTATCGAGAGCGGCCGCTGCTGATC"),
                                toBytes("W[`bbbb_baS\\`_\\bbabbaWR`bba``ab_bbbbbabbbaabb^^^\\aa_ab`a_`[`VaST[^SWTXWNYEHM")),
                        null, null),
                makeClusterData(true,
                        makeReadData(toBytes("CACACACACACACACACACACACCACCTTTTGGCTTATCTGCACGCGGCCGCGTGCCCTACCCTACCCCATGGGAT"),
                                toBytes("a_aa^\\Ra\\`aaXa_aa_aaaaaaa[^X^``V[`_`a^aaO``^_SJUELTVMKVTPFOKJJMNKTRRJEEPFKTR")),
                        null, null),
                new IlluminaDataType[]{},
                0, 0, 0
            },


            {"Non-PE Barcode Bustard Parsing Test", 5, 20, false,
                makeClusterData(true,
                        makeReadData(toBytes("GACTTTGGGAAGGGTCATTACTGCCCTTGTAGAAAGAACACCTCATGTTCCTTATCGAGAGCGGCCGCTG"),
                                toBytes("W[`bbbb_baS\\`_\\bbabbaWR`bba``ab_bbbbbabbbaabb^^^\\aa_ab`a_`[`VaST[^SWTX")),
                        null,
                        makeReadData(toBytes("CTGATC"),
                                toBytes("WNYEHM"))),
                makeClusterData(true,
                        makeReadData(toBytes("CACACACACACACACACACACACCACCTTTTGGCTTATCTGCACGCGGCCGCGTGCCCTACCCTACCCCA"),
                                toBytes("a_aa^\\Ra\\`aaXa_aa_aaaaaaa[^X^``V[`_`a^aaO``^_SJUELTVMKVTPFOKJJMNKTRRJE")),
                        null,
                        makeReadData(toBytes("TGGGAT"),
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
                makeClusterData(false,
                        makeReadData(toBytes("TAGAGATGGC.CT.........T......C........G...TCCAGACCGCCCATTCTCTGCCTGCC"),
                                toBytes("^`^BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")),
                        makeReadData(toBytes("CCTCTAATCCCAGCACTATCCGAGACCAAATCAGGCAAATCACTTGAAGTCAGGAGTTCGAGACCAGC"),
                                toBytes("]]VISQK_M\\`MHX\\ZFMaPWHXYUa]ZHJGaULGPXTRS\\W[`ZH_GMUa_[]M]PTUZX]VaZaU\\")),
                        makeReadData(toBytes("CCACCCAC"),
                                toBytes("_ZFZ^]BB"))),
                makeClusterData(false,
                        makeReadData(toBytes("CAACTCTTGT.GT........GT......A........G...AATATATTCTGAAACTCAGCAATGTT"),
                                toBytes("aaabaaBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")),
                        makeReadData(toBytes("TAACTTTCAGAGGCCCTTCAGGAGGCCCTGGCCTGTCAAGTACCTTTACAGTGATGGGTATAGACTTT"),
                                toBytes("abbbabba_bbabb_]S]ab_ab__^abSORYX^RF[aa`_a_[XVa`[WUN`a^__a\\ZL\\BBBBBB")),
                        makeReadData(toBytes("TCGGAATG"),
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

        final IlluminaDataProvider parser = factory.makeDataProvider();
        for (int i = 0; parser.hasNext(); ++i) {
            final ClusterData cluster = parser.next();
            final String matchedBarcode = cluster.getMatchedBarcode();
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

    private ReadData makeReadData(final byte[] bases, final byte[] qualities) {
        final ReadData readData = new ReadData();
        readData.setBases(bases);
        readData.setQualities(solexaQualityCharsToPhredBinary(qualities));
        return readData;
    }

    private ReadData makeReadData(final byte[] bases, final byte[] qualities, final FourChannelIntensityData rawIntensities, final FourChannelIntensityData noise) {
        final ReadData readData = makeReadData(bases, qualities);
        readData.setRawIntensities(rawIntensities);
        readData.setNoise(noise);
        return readData;
    }

    private ClusterData makeClusterData(final boolean pf, final ReadData firstEnd, final ReadData secondEnd, final ReadData barcode) {
        final ClusterData cluster = new ClusterData();
        cluster.setPf(pf);
        cluster.setFirstEnd(firstEnd);
        cluster.setSecondEnd(secondEnd);
        cluster.setBarcodeRead(barcode);

        return cluster;
    }
}
