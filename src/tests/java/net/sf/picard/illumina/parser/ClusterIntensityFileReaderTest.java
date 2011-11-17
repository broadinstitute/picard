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
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/");
    private static final File INTENSITY_FILES_DIR = new File(TEST_DATA_DIR, "L001/C1.1");

    private static final int NUM_BASES = 2;
    private static final int NUM_CLUSTERS = 153010;
    private static enum DATA_FIELDS {CLUSTER_INDEX, A_CHANNEL, C_CHANNEL, G_CHANNEL, T_CHANNEL}

    //Cluster#, Expected A,C,G,T
    private static final int[][] CNF_DATA = new int[][]{
                   //cluster    A   C       G   T
        new int[]{0,         0,  9,      0,  0},
        new int[]{1,         0,  15,     0,  0},
        new int[]{2,         0,  15,     0,  11},
        new int[]{51,        9,  23,     0,  9},
        new int[]{10001,     8,  22,     7,  11},
        new int[]{10002,     10, 9,      7,  9},
        new int[]{80000,     13, 55,     12, 76},
        new int[]{80001,     10, 54,     12, 73},
        new int[]{100001,    12, 33,     9,  45},
        new int[]{100010,    11, 52,     11, 70},
        new int[]{153009,    0,  0,      0,  0}
    };

    //Cluster#, Expected A,C,G,T
    private static final int[][] CIF_DATA = new int[][]{
                   //cluster    A      C       G     T
        new int[]{0,         0,     0,      0,    0},
        new int[]{1,         0,     0,      0,    0},
        new int[]{2,         0,     0,      0,    0},
        new int[]{20001,     436,   333,    -164, -63},
        new int[]{20002,     476,   354,    -31,  -61},
        new int[]{90001,     59,    332,    -14,  -240},
        new int[]{90002,     -91,   -29,    -330, 1090},
        new int[]{153008,    0,     0,      0,    0},
        new int[]{153009,    0,     0,      0,    0}
    };

    @Test(dataProvider="data")
    public void testBasic(final File cifFile, final int [][] goldData) throws Exception {
        final ClusterIntensityFileReader cifr = new ClusterIntensityFileReader(cifFile);
        System.out.println("File: " + cifr.getFile() + "; First cycle: " + cifr.getFirstCycle() +
                "; Num cycles: " + cifr.getNumCycles() + "; Num clusters: " + cifr.getNumClusters());

        for (int [] goldRow : goldData) {
            final int cluster = goldRow[DATA_FIELDS.CLUSTER_INDEX.ordinal()];
            Assert.assertEquals(1, cifr.getNumCycles(), "Expected only 1 cycle, cluster: " + cluster);

            final int cycle = cifr.getFirstCycle();
            final String msg = "Cycle: " + cycle + "; Cluster: " + cluster;
            Assert.assertEquals(cifr.getValue(cluster, IntensityChannel.A, cycle), goldRow[DATA_FIELDS.A_CHANNEL.ordinal()], msg);
            Assert.assertEquals(cifr.getValue(cluster, IntensityChannel.C, cycle), goldRow[DATA_FIELDS.C_CHANNEL.ordinal()], msg);
            Assert.assertEquals(cifr.getValue(cluster, IntensityChannel.G, cycle), goldRow[DATA_FIELDS.G_CHANNEL.ordinal()], msg);
            Assert.assertEquals(cifr.getValue(cluster, IntensityChannel.T, cycle), goldRow[DATA_FIELDS.T_CHANNEL.ordinal()], msg);
        }

        Assert.assertEquals(NUM_CLUSTERS, cifr.getNumClusters());
    }

    @DataProvider(name = "data")
    private Object[][] getTestData()
    {
        return new Object[][]{
                {new File(INTENSITY_FILES_DIR, "s_1_1.cnf"), CNF_DATA},
                {new File(INTENSITY_FILES_DIR, "s_1_1.cif"), CIF_DATA},
        };
    }
}
