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

import net.sf.picard.util.SolexaQualityConverter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Map;
import static net.sf.picard.illumina.parser.TestDataUtil.*;
import static net.sf.picard.illumina.parser.TestDataUtil.getTiledReadData;
import static net.sf.picard.util.CollectionUtil.*;

/**
* @author alecw@broadinstitute.org
*/

public class IlluminaDataProviderTest {

    public static final File PARSING_TEST_BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaBarcodeParsingTest/BaseCalls");
    public static final File TEST_DATA_LOCATION = new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir");
    public static final String RUN_BARCODE = "305PJAAXX080716";
    private static final File RTA_BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/rta/BasecallsDir");
    private static final IlluminaDataType [] DEFAULT_DATA_TYPES = new IlluminaDataType[]{
            IlluminaDataType.Position, IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF
    };

    @Test(dataProvider="data")
    public void testIlluminaDataProviderOldMethod(
            final String testName, final int lane, final int size, final boolean pe,
            final Map<Integer, ClusterData> readNoToClusterData,
            final IlluminaDataType [] extraDataTypes,
            final String illuminaConfigStr,
            final int barcodeCycle, final int barcodeLength,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final File basecallsDirectory)
            throws Exception {
        final IlluminaDataType [] dts = getDataTypes(extraDataTypes);
        final IlluminaDataProviderFactory factory;
        if (barcodeLength == 0) {
            factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, null, dts);
        } else {
            factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, barcodeCycle, barcodeLength, dts);
        }
        final IlluminaDataProvider dataProvider   = factory.makeDataProvider();

        runTest(testName, size, pe, illuminaConfigStr, readNoToClusterData, seekAfterFirstRead, seekTestDataReadOffset, dataProvider);
    }


    @Test(dataProvider="data")
    public void testIlluminaDataProviderNewMethod(
            final String testName, final int lane, final int size, final boolean pe,
            final Map<Integer, ClusterData> readNoToClusterData,
            final IlluminaDataType [] extraDataTypes,
            final String illuminaConfigStr,
            final int barcodeCycle, final int barcodeLength,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final File basecallsDirectory)
            throws Exception {

        final IlluminaDataType [] dts = getDataTypes(extraDataTypes);
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, new ReadStructure(illuminaConfigStr), dts);
        final IlluminaDataProvider dataProvider   = factory.makeDataProvider();

        runTest(testName, size, pe, illuminaConfigStr, readNoToClusterData, seekAfterFirstRead, seekTestDataReadOffset, dataProvider);
    }

    private void runTest(
            final String testName,  final int size, final boolean pe,
            final String illuminaConfigStr,
            final Map<Integer, ClusterData> readNoToClusterData,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final IlluminaDataProvider dataProvider)
            throws Exception {

        Assert.assertEquals(new ReadStructure(illuminaConfigStr), dataProvider.getReadStructure());

        int count = 0;
        int readNum = 0;
        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();
            if (readNoToClusterData.containsKey(readNum)) {
                compareReadData(pe, cluster, readNoToClusterData.get(readNum), testName);
            }

            if (seekAfterFirstRead != 0 && count == 0) {
                dataProvider.seekToTile(seekAfterFirstRead);
                readNum += seekTestDataReadOffset;
            }

            readNum++;
            count++;
        }
        Assert.assertEquals(count, size, testName);
    }

    private IlluminaDataType [] getDataTypes(IlluminaDataType [] extraDataTypes) {
        final IlluminaDataType [] dts;

        if(extraDataTypes == null) {
            dts = DEFAULT_DATA_TYPES;
        } else {
            dts = Arrays.copyOf(DEFAULT_DATA_TYPES, DEFAULT_DATA_TYPES.length + extraDataTypes.length);
            System.arraycopy(extraDataTypes, 0, dts, DEFAULT_DATA_TYPES.length, extraDataTypes.length);
        }
        return dts;
    }


    private void compareBasesAndQuals(final ReadData rd1, final ReadData rd2, final String testName) {
        Assert.assertEquals(rd1.getBases(),     rd2.getBases(),     testName);
        Assert.assertEquals(rd1.getQualities(), rd2.getQualities(), testName);
    }

    private void comparePositionalData(final ClusterData cd1, final ClusterData cd2) {
        Assert.assertEquals(cd1.getLane(), cd2.getLane());
        Assert.assertEquals(cd1.getTile(), cd2.getTile());
        Assert.assertEquals(cd1.getX(),    cd2.getX());
        Assert.assertEquals(cd1.getY(),    cd2.getY());
    }

    //Doesn't compare intensities right now -- Do we want too?
    private void compareReadData(final boolean pe, final ClusterData cd1, final ClusterData cd2, final String testName) {
        comparePositionalData(cd1, cd2);
        Assert.assertEquals(cd1.getNumReads(), cd2.getNumReads());
        for(int i = 0; i < cd1.getNumReads(); i++) {
            compareBasesAndQuals(cd1.getRead(i), cd2.getRead(i), testName);
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
            {
              "PE Bustard Parsing Test", 1, 60, true,
                qseqDataToClusterMap(getTiledReadData(s_1_1, makeList(1, 2, 3)), getTiledReadData(s_1_2, makeList(1,2,3)), null, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "76T76T",
                0, 0,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE with Barcode Bustard Parsing Test", 1, 60, true,
                qseqDataToClusterMap(getTiledReadData(s_1_1, makeList(1, 2, 3)), getTiledReadData(s_1_2, makeList(1,2,3)), 77, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "76T6B70T",
                77, 6,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20, true,
                qseqDataToClusterMap(getReadData(s_8_1_0001), getReadData(s_8_2_0001), null, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                "4T4T",
                0, 0,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20, true,
                qseqDataToClusterMap(getReadData(s_8_1_0001), getReadData(s_8_2_0001), 5, 2, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                "4T2B2T",
                5, 2,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE, Barcode, seek Bustard Parsing Test", 1, 21, true,
                qseqDataToClusterMap(getTiledReadData(s_1_1, makeList(1, 2, 3)), getTiledReadData(s_1_2, makeList(1,2,3)), 77, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "76T6B70T",
                77, 6,
                3, getFileSize(s_1_1_0001) + getFileSize(s_1_2_0001) - 1,
                TEST_DATA_LOCATION
            },

            {"Non-PE Bustard Parsing Test", 5, 20, false,
                qseqDataToClusterMap(getReadData(s_5_1_0001), null, null, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "76T",
                0, 0,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"Non-PE Bustard Parsing Test", 5, 20, false,
                qseqDataToClusterMap(getReadData(s_5_1_0001), null, 71, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "70T6B",
                71, 6,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"Barcode-aware PE Bustard Parsing Test", 6, 10, true,
                qseqDataToClusterMap(getReadData(s_6_1_0001), getReadData(s_6_3_0001),  getReadData(s_6_2_0001), DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "68T8B68T",
                69, 8,
                0, 0,
                TEST_DATA_LOCATION
            }
        };
    }

    public void runBarcodeParsingTest(final IlluminaDataProviderFactory factory) {

        final IlluminaDataProvider dateProvider = factory.makeDataProvider();
        for (int i = 0; dateProvider.hasNext(); ++i) {
            final ClusterData cluster = dateProvider.next();
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
    public void barcodeParsingTestNoConfig() {
        runBarcodeParsingTest(new IlluminaDataProviderFactory(PARSING_TEST_BASECALLS_DIR, 6, 153, 6, IlluminaDataType.Barcodes));
    }

    @Test
    public void barcodeParsingTestWConfig() {
        runBarcodeParsingTest(new IlluminaDataProviderFactory(PARSING_TEST_BASECALLS_DIR, 6, new ReadStructure("76T76T6B"), IlluminaDataType.Barcodes));
    }

}
