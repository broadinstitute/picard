/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
public class ClusterIntensityFileReaderTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/rta");
    private static final File INTENSITY_FILES_DIR = new File(TEST_DATA_DIR, "L001/C1.1");

    private static final int NUM_BASES = 2;

    @Test(dataProvider="data")
    public void testBasic(final File cifFile, final String fileType) throws Exception {
        final ClusterIntensityFileReader cifr = new ClusterIntensityFileReader(cifFile);
        System.out.println("File: " + cifr.getFile() + "; First cycle: " + cifr.getFirstCycle() +
                "; Num cycles: " + cifr.getNumCycles() + "; Num clusters: " + cifr.getNumClusters());
        System.out.println("Text file type: " + fileType);
        final ReadConfiguration configuration = new ReadConfiguration();
        configuration.setFirstStart(1);
        configuration.setFirstEnd(NUM_BASES);
        // This parser takes whatever file type you give it, but then stuff the values into the raw intensities slot.
        final HackyIntParser intParser = new HackyIntParser(configuration, TEST_DATA_DIR, 1, fileType);
        int cluster;
        for (cluster = 0; intParser.hasNext(); ++cluster) {
            final IlluminaReadData read = new IlluminaReadData();
            read.setFirstEnd(new IlluminaEndData());
            intParser.next(read);
            final FourChannelIntensityData rawIntensities = read.getFirstEnd().getRawIntensities();
            for (int cycle = cifr.getFirstCycle(); cycle < cifr.getFirstCycle() + cifr.getNumCycles(); ++cycle) {
                final String msg = "Cycle: " + cycle + "; Cluster: " + cluster;
                Assert.assertEquals(cifr.getValue(cluster, IntensityChannel.A_CHANNEL, cycle),
                        rawIntensities.getA()[cycle-1], msg);
                Assert.assertEquals(cifr.getValue(cluster, IntensityChannel.C_CHANNEL, cycle),
                        rawIntensities.getC()[cycle-1], msg);
                Assert.assertEquals(cifr.getValue(cluster, IntensityChannel.G_CHANNEL, cycle),
                        rawIntensities.getG()[cycle-1], msg);
                Assert.assertEquals(cifr.getValue(cluster, IntensityChannel.T_CHANNEL, cycle),
                        rawIntensities.getT()[cycle-1], msg);
            }
        }
        Assert.assertEquals(cluster, cifr.getNumClusters());
    }

    @DataProvider(name = "data")
    private Object[][] getTestData()
    {
        return new Object[][]{
                {new File(INTENSITY_FILES_DIR, "s_1_1.cnf"),
                 "nse"},
                {new File(INTENSITY_FILES_DIR, "s_1_1.cif"),
                 "int"},
        };
    }

    /**
     * Little kludge to open either nse or int file, but have the values stuffed into rawIntensities slot.
     */
    public class HackyIntParser extends IntensitiesOrNoiseParser{
        public HackyIntParser(final ReadConfiguration readConfiguration, final File directory, final int lane,
                              final String fileType) {
            super(readConfiguration, directory, lane, fileType, null);
        }

        @Override
        protected void setValues(final IlluminaEndData end, final FourChannelIntensityData values) {
            end.setRawIntensities(values);
        }
    }
}
