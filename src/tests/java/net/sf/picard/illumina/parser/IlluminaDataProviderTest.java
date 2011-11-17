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

import javax.xml.bind.annotation.XmlElementRef;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
    public void testIlluminaDataProvider(
            final String testName, final int lane, final int size, final boolean pe,
            final Map<Integer, ClusterData> readNoToClusterData,
            final IlluminaDataType [] extraDataTypes,
            final int barcodeCycle, final int barcodeLength,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final File basecallsDirectory)
            throws Exception {


        final IlluminaDataType [] dts;

        if(extraDataTypes == null) {
            dts = DEFAULT_DATA_TYPES;
        } else {
            dts = Arrays.copyOf(DEFAULT_DATA_TYPES, DEFAULT_DATA_TYPES.length + extraDataTypes.length);
            System.arraycopy(extraDataTypes, 0, dts, DEFAULT_DATA_TYPES.length, extraDataTypes.length);
        }
        
        final IlluminaDataProviderFactory factory;
        if (barcodeLength == 0) {
            factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, DEFAULT_DATA_TYPES);
        } else {
            factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, barcodeCycle, barcodeLength, DEFAULT_DATA_TYPES);
        }

        final IlluminaDataProvider parser = factory.makeDataProvider();

        int count = 0;
        int readNum = 0;
        while (parser.hasNext()) {
            final ClusterData cluster = parser.next();
            if (readNoToClusterData.containsKey(readNum)) {
                compareReadData(pe, cluster, readNoToClusterData.get(readNum), testName);
            }

            if (seekAfterFirstRead != 0 && count == 0) {
                parser.seekToTile(seekAfterFirstRead);
                readNum += seekTestDataReadOffset;
            }

            readNum++;
            count++;
        }
        Assert.assertEquals(count, size, testName);
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
                0, 0,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE with Barcode Bustard Parsing Test", 1, 60, true,
                qseqDataToClusterMap(getTiledReadData(s_1_1, makeList(1, 2, 3)), getTiledReadData(s_1_2, makeList(1,2,3)), 77, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                77, 6,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20, true,
                qseqDataToClusterMap(getReadData(s_8_1_0001), getReadData(s_8_2_0001), null, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                0, 0,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE Bustard Parsing Test with Noise/Intensity", 8, 20, true,
                qseqDataToClusterMap(getReadData(s_8_1_0001), getReadData(s_8_2_0001), 5, 2, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{IlluminaDataType.RawIntensities, IlluminaDataType.Noise},
                5, 2,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"PE, Barcode, seek Bustard Parsing Test", 1, 21, true,
                qseqDataToClusterMap(getTiledReadData(s_1_1, makeList(1, 2, 3)), getTiledReadData(s_1_2, makeList(1,2,3)), 77, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                77, 6,
                3, getFileSize(s_1_1_0001) + getFileSize(s_1_2_0001) - 1,
                TEST_DATA_LOCATION
            },

            {"Non-PE Bustard Parsing Test", 5, 20, false,
                qseqDataToClusterMap(getReadData(s_5_1_0001), null, null, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                0, 0,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"Non-PE Bustard Parsing Test", 5, 20, false,
                qseqDataToClusterMap(getReadData(s_5_1_0001), null, 71, 6, DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                71, 6,
                0, 0,
                TEST_DATA_LOCATION
            },

            {"Barcode-aware PE Bustard Parsing Test", 6, 10, true,
                qseqDataToClusterMap(getReadData(s_6_1_0001), getReadData(s_6_3_0001),  getReadData(s_6_2_0001), DEFAULT_DATA_TYPES),
                new IlluminaDataType[]{},
                69, 8,
                0, 0,
                TEST_DATA_LOCATION
            },

        };
    }

    @Test
    public void testBarcodeParsing() {
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(PARSING_TEST_BASECALLS_DIR, 6, 153, 6, IlluminaDataType.Barcodes);

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

    /* //Put this back in, so as not to change Public API, or maybe provide a static method that instead you provide, tile, lane, directory
       //to keep the factory immutable and still be able to query for available data
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
    } */

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
        final ClusterData cluster = new ClusterData(firstEnd, barcode, secondEnd);
        cluster.setPf(pf);

        return cluster;
    }
}
