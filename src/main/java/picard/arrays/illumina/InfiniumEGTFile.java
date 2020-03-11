/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

package picard.arrays.illumina;
import org.apache.commons.io.IOUtils;
import picard.PicardException;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * A class to parse the contents of an Illumina Infinium cluster (EGT) file
 *
 * A cluster file contains information about the clustering information used in mapping red / green intensity information
 * to genotype calls
 */
public class InfiniumEGTFile extends InfiniumDataFile {
    public static final String EXTENSION = "egt";

    private static final int VALID_GENTRAIN_DATA_TYPE = 9;
    private static final int INVALID_FILE_VERSION = 2;

    public int[][] n;
    public float[][] meanR;
    public float[][] meanTheta;
    public float[][] devR;
    public float[][] devTheta;

    public float[] totalScore;

    public String[] rsNames;

    public String manifestName;

    public int numCodes;

    public Map<String, Integer> rsNameToIndex;


    public InfiniumEGTFile(final File clusterFile) throws IOException {
        super(new DataInputStream(new FileInputStream(clusterFile)), false);
        parse();
    }

    private void parse() throws IOException {
        try {
            readHeaderData();
            readFileData();
        } finally {
            IOUtils.closeQuietly(stream);
        }
    }

    private void readFileData() throws IOException {
        final int version = parseInt();
        if (version > VALID_GENTRAIN_DATA_TYPE) {
            throw new IOException("Error. Cannot read file - unknown version " + version + " for gentrain data type");
        }

        manifestName = parseString();
        numCodes = parseInt();

        initializeArrays(numCodes);

        for (int i = 0; i < numCodes; i++) {
            parseInts(n[i]);
            parseFloats(devR[i]);
            parseFloats(meanR[i]);
            parseFloats(devTheta[i]);
            parseFloats(meanTheta[i]);

            // 15 unused floats
            skipFloats(15);
        }
        for (int i = 0; i < numCodes; i++) {
            // skip cSepScore
            skipFloat();

            totalScore[i] = parseFloat();

            // skip original score
            skipFloat();

            skipBoolean();
        }
        for (int i = 0; i < numCodes; i++) {
            skipString();
        }
        for (int i = 0; i < numCodes; i++) {
            rsNames[i] = parseString();
            if (rsNameToIndex.put(rsNames[i], i) != null) {
                throw new PicardException("Non-unique rsName '" + rsNames[i] + "' found in cluster file");
            }
        }
    }

    private void initializeArrays(int numCodes) {

        n = new int[numCodes][InfiniumVcfFields.NUM_GENOTYPE_VALUES];
        devR = new float[numCodes][InfiniumVcfFields.NUM_GENOTYPE_VALUES];
        meanR = new float[numCodes][InfiniumVcfFields.NUM_GENOTYPE_VALUES];
        devTheta = new float[numCodes][InfiniumVcfFields.NUM_GENOTYPE_VALUES];
        meanTheta = new float[numCodes][InfiniumVcfFields.NUM_GENOTYPE_VALUES];

        totalScore = new float[numCodes];
        rsNames = new String[numCodes];
        rsNameToIndex = new HashMap<>();
    }

    private void parseInts(int[] array) throws IOException {
        for (int i = 0; i < array.length; i++) {
            array[i] = parseInt();
        }
    }

    private void parseFloats(float[] array) throws IOException {
        for (int i = 0; i < array.length; i++) {
            array[i] = parseFloat();
        }
    }

    private void readHeaderData() throws IOException {
        setFileVersion(parseInt());
        // skip gcVersion
        skipString();
        // skip clusterVersion
        skipString();
        // skip callVersion
        skipString();
        // skip normalizationVersion
        skipString();
        // skip dataCreated
        skipString();
        // skip isWGT
        skipBoolean();

        if (getFileVersion() == INVALID_FILE_VERSION) {
            throw new IOException("Version '" + INVALID_FILE_VERSION + "' unsupported");
        }
        manifestName = parseString();
    }
}
