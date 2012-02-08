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

import net.sf.picard.PicardException;
import net.sf.picard.util.SolexaQualityConverter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import static net.sf.picard.illumina.parser.QSeqTdUtil.*;
import static net.sf.picard.illumina.parser.QSeqTdUtil.getTiledReadData;
import static net.sf.picard.util.CollectionUtil.*;

/**
* @author jburke@broadinstitute.org
*/

public class IlluminaDataProviderTest {

    public static final File PARSING_TEST_BASECALLS_DIR = new File("testdata/net/sf/picard/illumina/IlluminaBarcodeParsingTest/BaseCalls");
    public static final File TEST_DATA_LOCATION = new File("testdata/net/sf/picard/illumina/IlluminaTests/BasecallsDir");
    public static final File BINARY_TD_LOCATION = new File("testdata/net/sf/picard/illumina/CompleteIlluminaDir/Intensities/BaseCalls");
    private static final IlluminaDataType [] DEFAULT_DATA_TYPES = new IlluminaDataType[]{
            IlluminaDataType.Position, IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF
    };

    @Test(dataProvider="data")
    public void testIlluminaDataProviderQSeqMethod(
            final String testName, final int lane, final int size,
            final Map<Integer, ClusterData> readNoToClusterData,
            final IlluminaDataType[] extraDataTypes,
            final String illuminaConfigStr,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final File basecallsDirectory)
            throws Exception {

        final IlluminaDataType [] dts = getDataTypes(extraDataTypes);
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, new ReadStructure(illuminaConfigStr), dts);
        final IlluminaDataProvider dataProvider   = factory.makeDataProvider();

        runTest(testName, size, readNoToClusterData, seekAfterFirstRead, seekTestDataReadOffset, dataProvider);
    }

    private void runTest(
            final String testName,  final int size,
            final Map<Integer, ClusterData> readNoToClusterData,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final IlluminaDataProvider dataProvider)
            throws Exception {

        int count = 0;
        int readNum = 0;
        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();
            if (readNoToClusterData.containsKey(readNum)) {
                compareReadData(cluster, readNoToClusterData.get(readNum), testName + " cluster num " + readNum);
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
        Assert.assertEquals(rd1.getReadType(),  rd2.getReadType());
    }

    private void comparePositionalData(final ClusterData cd1, final ClusterData cd2, final String testName) {
        Assert.assertEquals(cd1.getLane(), cd2.getLane(), testName);
        Assert.assertEquals(cd1.getTile(), cd2.getTile(), testName);
        Assert.assertEquals(cd1.getX(),    cd2.getX(), testName);
        Assert.assertEquals(cd1.getY(),    cd2.getY(), testName);
    }

    //Doesn't compare intensities right now -- Do we want too?
    private void compareReadData(final ClusterData cd1, final ClusterData cd2, final String testName) {
       comparePositionalData(cd1, cd2, testName);
        Assert.assertEquals(cd1.getNumReads(), cd2.getNumReads());
        for(int i = 0; i < cd1.getNumReads(); i++) {
            compareBasesAndQuals(cd1.getRead(i), cd2.getRead(i), testName);
        }

        Assert.assertEquals(cd1.getMatchedBarcode(), cd2.getMatchedBarcode(), testName);
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
              "PE Bustard Parsing Test", 1, 60,
                qseqDataToClusterMap(getTiledReadData(s_1_1, makeList(1, 2, 3)), getTiledReadData(s_1_2, makeList(1,2,3)), null, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "76T76T",
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE with Barcode Bustard Parsing Test", 1, 60,
                qseqDataToClusterMap(getTiledReadData(s_1_1, makeList(1, 2, 3)), getTiledReadData(s_1_2, makeList(1,2,3)), 77, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "76T6B70T",
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20,
                qseqDataToClusterMap(getReadData(s_8_1_0001), getReadData(s_8_2_0001), null, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                "4T4T",
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20,
                qseqDataToClusterMap(getReadData(s_8_1_0001), getReadData(s_8_2_0001), 5, 2, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                "4T2B2T",
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE, Barcode, seek Bustard Parsing Test", 1, 21,
                qseqDataToClusterMap(getTiledReadData(s_1_1, makeList(1, 2, 3)), getTiledReadData(s_1_2, makeList(1,2,3)), 77, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "76T6B70T",
                3, getFileSize(s_1_1_0001) + getFileSize(s_1_2_0001) - 1,
                TEST_DATA_LOCATION
            },

            {"Non-PE Bustard Parsing Test", 5, 20,
                qseqDataToClusterMap(getReadData(s_5_1_0001), null, null, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "76T",
                0, 0,
                TEST_DATA_LOCATION
            },

            {"Non-PE Bustard Parsing Test", 5, 20,
                qseqDataToClusterMap(getReadData(s_5_1_0001), null, 71, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "70T6B",
                0, 0,
                TEST_DATA_LOCATION
            },

            {"Barcode-aware PE Bustard Parsing Test", 6, 10,
                qseqDataToClusterMap(getReadData(s_6_1_0001), getReadData(s_6_3_0001),  getReadData(s_6_2_0001), DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                "68T8B68T",
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
    public void barcodeParsingTest() {
        runBarcodeParsingTest(new IlluminaDataProviderFactory(PARSING_TEST_BASECALLS_DIR, 6, new ReadStructure("76T76T6B"), IlluminaDataType.Barcodes));
    }

    @DataProvider(name="binaryData")
    public Object[][] binaryData() {
        return new Object[][] {
            {
                "Bustard Parsing Test(25T8B25T) w/Clocs", 1, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{IlluminaDataType.Barcodes}, //ADD TESTING FOR BARCODES
                "25T8B25T",
                0, 0,
                BINARY_TD_LOCATION
            },
            {
                "Bustard Parsing Test(25T8S25T) w/Clocs", 1, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{IlluminaDataType.Barcodes}, //ADD TESTING FOR BARCODES
                "25T8S25T",
                0, 0,
                BINARY_TD_LOCATION
            },
            {
                "Bustard Parsing Test(25S8S25T) w/Clocs", 1, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{IlluminaDataType.Barcodes}, //ADD TESTING FOR BARCODES
                "25S8S25T",
                0, 0,
                BINARY_TD_LOCATION
            },
            {
                "Bustard Parsing Test(25T8B25T) w/Clocs And Seeking", 1, 21,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{IlluminaDataType.Barcodes},
                "25T8B25T",
                2101, 39,
                BINARY_TD_LOCATION
            },
            {
                "Bustard Parsing Test(25T8B25T) w/Pos", 2, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{},
                "25T8B25T",
                0, 0,
                BINARY_TD_LOCATION
            },
            {
                "Bustard Parsing Test(25T8B25T) w/Pos and Seeking", 2, 41,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{},
                "25T8B25T",
                1201, 19,
                BINARY_TD_LOCATION
            },
            {
                "Bustard Parsing Test(19T10B10B19T) w/Pos.gz", 4, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{},
                "19T10B10B19T",
                0, 0,
                BINARY_TD_LOCATION
            },
            {
                "Bustard Parsing Test(19T10B10B19T) w/Pos.gz and Seeking", 4, 41,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{},
                "19T10B10B19T",
                1201, 19,
                BINARY_TD_LOCATION
            },
            {
                "Bustard Parsing Test(17T2S10B10B17T2S) w/Pos.gz and Seeking", 4, 41,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{},
                "17T2S10B10B17T2S",
                1201, 19,
                BINARY_TD_LOCATION
            }
        };
    }

    @Test(dataProvider="binaryData")
    public void testIlluminaDataProviderBclMethod(
            final String testName, final int lane, final int size,
            final List<Integer> tiles,
            final IlluminaDataType[] extraDataTypes,
            final String illuminaConfigStr,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final File basecallsDirectory)
            throws Exception {

        final IlluminaDataType [] dts = getDataTypes(extraDataTypes);

        Map<Integer, ClusterData> readNoToClusterData = BinTdUtil.clusterData(lane, tiles, illuminaConfigStr, dts);
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, new ReadStructure(illuminaConfigStr), dts);
        final IlluminaDataProvider dataProvider   = factory.makeDataProvider();

        runTest(testName, size, readNoToClusterData, seekAfterFirstRead, seekTestDataReadOffset, dataProvider);
    }

    //Unlike above, the data types here do not have DEFAULT_DATA_TYPES added before createing the dataProvider
    @DataProvider(name="badData")
    public Object[][] badData() {
        return new Object[][] {
            {
                "Bad Lane(5)", 5, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{IlluminaDataType.Barcodes},
                "25T8B25T",
                BINARY_TD_LOCATION
            },
            {
                "Bad Read Structure(25TB25T)", 4, 60,
                makeList(1101, 1201, 2101),
                DEFAULT_DATA_TYPES,
                "25TB25T",
                BINARY_TD_LOCATION
            },
            {
                "Bad Read Structure(25T0B25T)", 4, 60,
                makeList(1101, 1201, 2101),
                DEFAULT_DATA_TYPES,
                "25T0B25T",
                BINARY_TD_LOCATION
            },
            {
                "Bad Read Structure(-225T0B25T)", 4, 60,
                makeList(1101, 1201, 2101),
                DEFAULT_DATA_TYPES,
                "-25T0B25T",
                BINARY_TD_LOCATION
            },
            {
                "Missing Barcodes File", 9, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{IlluminaDataType.Position, IlluminaDataType.Barcodes},
                "25T8B25T",
                BINARY_TD_LOCATION
            },
            {
                "Missing Cycle File", 9, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{IlluminaDataType.BaseCalls},
                "25T8B25T",
                BINARY_TD_LOCATION
            },
            {
                "Missing Filter File", 9, 60,
                makeList(1101, 1201, 2101),
                new IlluminaDataType[]{IlluminaDataType.PF, IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores},
                "25T8B24T",
                BINARY_TD_LOCATION
            }
        };
    }

    @Test(dataProvider="badData", expectedExceptions = {PicardException.class, IllegalArgumentException.class})
    public void testIlluminaDataProviderMissingDatas(
            final String testName, final int lane, final int size,
            final List<Integer> tiles,
            final IlluminaDataType[] actualDts,
            final String illuminaConfigStr,
            final File basecallsDirectory)
            throws Exception {
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, new ReadStructure(illuminaConfigStr), actualDts);
        final IlluminaDataProvider dataProvider   = factory.makeDataProvider();
    }
}
