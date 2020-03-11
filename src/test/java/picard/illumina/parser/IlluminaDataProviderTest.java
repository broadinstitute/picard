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
package picard.illumina.parser;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static htsjdk.samtools.util.CollectionUtil.makeList;
//import static htsjdk.samtools.util.CollectionUtil.*;

/**
 * @author jburke@broadinstitute.org
 */

public class IlluminaDataProviderTest {

    public static final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY);
    public static final File BINARY_TD_LOCATION = new File("testdata/picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final IlluminaDataType[] DEFAULT_DATA_TYPES = new IlluminaDataType[]{
            IlluminaDataType.Position, IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF
    };

    private void runTest(
            final String testName, final int size,
            final Map<Integer, ClusterData> readNoToClusterData,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final BaseIlluminaDataProvider dataProvider)
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
        dataProvider.close();
    }

    private IlluminaDataType[] getDataTypes(final IlluminaDataType[] extraDataTypes) {
        final IlluminaDataType[] dts;

        if (extraDataTypes == null) {
            dts = DEFAULT_DATA_TYPES;
        } else {
            dts = Arrays.copyOf(DEFAULT_DATA_TYPES, DEFAULT_DATA_TYPES.length + extraDataTypes.length);
            System.arraycopy(extraDataTypes, 0, dts, DEFAULT_DATA_TYPES.length, extraDataTypes.length);
        }
        return dts;
    }


    private void compareBasesAndQuals(final ReadData rd1, final ReadData rd2, final String testName) {
        Assert.assertEquals(rd1.getBases(), rd2.getBases(), testName);
        Assert.assertEquals(rd1.getQualities(), rd2.getQualities(), testName);
        Assert.assertEquals(rd1.getReadType(), rd2.getReadType());
    }

    private void comparePositionalData(final ClusterData cd1, final ClusterData cd2, final String testName) {
        Assert.assertEquals(cd1.getLane(), cd2.getLane(), testName);
        Assert.assertEquals(cd1.getTile(), cd2.getTile(), testName);
        Assert.assertEquals(cd1.getX(), cd2.getX(), testName);
        Assert.assertEquals(cd1.getY(), cd2.getY(), testName);
    }

    //Doesn't compare intensities right now -- Do we want too?
    private void compareReadData(final ClusterData cd1, final ClusterData cd2, final String testName) {
        comparePositionalData(cd1, cd2, testName);
        Assert.assertEquals(cd1.getNumReads(), cd2.getNumReads());
        for (int i = 0; i < cd1.getNumReads(); i++) {
            compareBasesAndQuals(cd1.getRead(i), cd2.getRead(i), testName);
        }

        Assert.assertEquals(cd1.getMatchedBarcode(), cd2.getMatchedBarcode(), testName);
        Assert.assertEquals(cd1.isPf().booleanValue(), cd2.isPf().booleanValue(), testName);
    }

    public void runBarcodeParsingTest(final IlluminaDataProviderFactory factory) {
        int total = 0;
        final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider();
        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();
            final String matchedBarcode = cluster.getMatchedBarcode();
            if (matchedBarcode != null) {
                Assert.assertEquals(matchedBarcode, new String(cluster.getRead(1).getBases()));
            }
            if(total > 10){
                break;
            }
            total++;
        }
        dataProvider.close();
    }

    @Test
    public void barcodeParsingTest() {
        runBarcodeParsingTest(new IlluminaDataProviderFactory(BINARY_TD_LOCATION, 1, new ReadStructure("25T8B25T"), bclQualityEvaluationStrategy, IlluminaDataType.BaseCalls,
                IlluminaDataType.Barcodes));
    }

    @DataProvider(name = "binaryData")
    public Object[][] binaryData() {
        return new Object[][]{
                {
                        "Bustard Parsing Test(25T8B25T) w/Clocs", 1, 180,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8B25T",
                        0, 0,
                        BINARY_TD_LOCATION
                },
                {
                        "Bustard Parsing Test(25T8S25T) w/Clocs", 1, 180,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8S25T",
                        0, 0,
                        BINARY_TD_LOCATION
                },
                {
                        "Bustard Parsing Test(25T8S25T) w/Clocs with ending skip", 1, 180,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8B1S",
                        0, 0,
                        BINARY_TD_LOCATION
                },
                {
                        "Bustard Parsing Test(25S8S25T) w/Clocs", 1, 180,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25S8S25T",
                        0, 0,
                        BINARY_TD_LOCATION
                },
                {
                        "Bustard Parsing Test(25T8B25T) w/Clocs And Seeking", 1, 61,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8B25T",
                        2101, 4631,
                        BINARY_TD_LOCATION
                }
        };
    }

    @Test(dataProvider = "binaryData")
    public void testIlluminaDataProviderBclMethod(
            final String testName, final int lane, final int size,
            final List<Integer> tiles,
            final IlluminaDataType[] extraDataTypes,
            final String illuminaConfigStr,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final File basecallsDirectory)
            throws Exception {

        final IlluminaDataType[] dts = getDataTypes(extraDataTypes);

        final Map<Integer, ClusterData> readNoToClusterData = BinTdUtil.clusterData(lane, tiles, illuminaConfigStr, dts);
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, new ReadStructure(illuminaConfigStr), bclQualityEvaluationStrategy, dts);
        final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider();

        runTest(testName, size, readNoToClusterData, seekAfterFirstRead, seekTestDataReadOffset, dataProvider);
    }

    //Unlike above, the data types here do not have DEFAULT_DATA_TYPES added before creating the dataProvider
    @DataProvider(name = "badData")
    public Object[][] badData() {
        return new Object[][]{
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

    @Test(dataProvider = "badData", expectedExceptions = {PicardException.class, IllegalArgumentException.class})
    public void testIlluminaDataProviderMissingDatas(
            final String testName, final int lane, final int size,
            final List<Integer> tiles,
            final IlluminaDataType[] actualDts,
            final String illuminaConfigStr,
            final File basecallsDirectory)
            throws Exception {
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, new ReadStructure(illuminaConfigStr), bclQualityEvaluationStrategy, actualDts);
        factory.makeDataProvider();
    }
}
