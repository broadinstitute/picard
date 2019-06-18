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


import picard.PicardException;

import java.io.DataInputStream;
import java.io.IOException;

/**
 * A class to parse the contents of an Illumina Infinium genotype (GTC) file
 *
 * A GTC file is the output of Illumina's genotype calling software (either Autocall or Autoconvert) and
 * contains genotype calls, confidence scores, basecalls and raw intensities for all calls made on the chip.
 *
 * This class will parse the binary GTC file format and allow access to the genotype, scores, basecalls and raw
 * intensities.
 */
public class InfiniumGTCFile extends InfiniumDataFile {

    private static final int NUM_SNPS = 1;
    private static final int PLOIDY = 2;
    private static final int PLOIDY_TYPE = 3;
    private static final int SAMPLE_NAME = 10;
    private static final int SAMPLE_PLATE = 11;
    private static final int SAMPLE_WELL = 12;
    private static final int CLUSTER_FILE = 100;
    private static final int SNP_MANIFEST = 101;
    private static final int IMAGING_DATE = 200;
    private static final int AUTOCALL_DATE = 201;
    private static final int AUTOCALL_VERSION = 300;
    private static final int TRANSFORMATIONS = 400;
    private static final int RAW_CONTROL_X_INTENSITIES = 500;
    private static final int RAW_CONTROL_Y_INTENSITIES = 501;
    private static final int RAW_X_INTENSITIES = 1000;
    private static final int RAW_Y_INTESITIES = 1001;
    private static final int GENOTYPES = 1002;
    private static final int BASE_CALLS = 1003;
    private static final int GENOTYPE_SCORES = 1004;
    private static final int SCANNER_INFO = 1005;
    private static final int CALL_RATE = 1006;
    private static final int GENDER = 1007;
    private static final int LOG_R_DEV = 1008;
    private static final int P_10_GC = 1009;
    private static final int DX = 1010;
    private static final int EXTENDED_SAMPLE_DATA = 1011;
    private static final int B_ALLELE_FREQS = 1012;
    private static final int LOG_R_RATIOS = 1013;
    private static final int INTENSITY_X_PERCENTILES = 1014;
    private static final int INTENSITY_Y_PERCENTILES = 1015;
    private static final int SENTRIX_ID = 1016;

    private static final int NO_CALL_CHAR = (int) '-';      // Used for string representation of no-call "--"
    private static final int IDENTIFIER_LENGTH = 3;
    private static final String GTC_IDENTIFIER = "gtc";

    public static final byte NO_CALL = 0;
    public static final byte AA_CALL = 1;
    public static final byte AB_CALL = 2;
    public static final byte BB_CALL = 3;

    private final InfiniumNormalizationManifest normalizationManifest;

    private int numberOfSnps;
    private int ploidy;
    // 1 = Diploid, 2 = Autopolyploid, 3 = Allopolyploid
    private int ploidyType;
    private String sampleName;
    private String samplePlate;
    private String sampleWell;
    private String clusterFile;
    private String snpManifest;
    private String imagingDate;
    private String autoCallDate;
    private String autoCallVersion;

    private InfiniumTransformation[] normalizationTransformations;
    private int[] rawControlXIntensities;
    private int[] rawControlYIntensities;
    private int[] rawXIntensities;
    private int[] rawYIntensities;
    private byte[] genotypeBytes;
    private float[] genotypeScores;
    private float[] bAlleleFreqs;
    private float[] logRRatios;
    private byte[][] baseCalls;

    //scanner data - 1005
    private String scannerName;
    private int pmtGreen;
    private int pmtRed;
    private String scannerVersion;
    private String imagingUser;

    private double callRate;
    private String gender;
    private float logRDev;
    private float p10GC;
    private int dx;

    //extended sample information - 1011
    private float p50GC;
    private int numCalls;
    private int numNoCalls;
    private int numIntensityOnly;

    //intensity X percentiles - 1014
    private IntensityPercentiles redIntensityPercentiles;

    //intensity Y percentiles - 1015
    private IntensityPercentiles greenIntensityPercentiles;

    private String sentrixBarcode;

    //derived values
    private float[] normalizedXIntensities;
    private float[] normalizedYIntensities;
    private float[] RIlmn;
    private float[] thetaIlmn;

    private int aaCalls = 0;
    private int abCalls = 0;
    private int bbCalls = 0;

    /**
     * Creates an InfiniumGTCFile object and parses the given input stream.
     *
     * @param gtcStream The gtc file input stream.
     * @throws IOException is thrown when there is a problem reading the stream.
     */
    public InfiniumGTCFile(final DataInputStream gtcStream, final InfiniumNormalizationManifest normalizationManifest) throws IOException {
        super(gtcStream, true);
        this.normalizationManifest = normalizationManifest;
        parse();
        normalizeAndCalculateStatistics();
    }

    InfiniumGTCFile(final DataInputStream gtcStream) throws IOException {
        super(gtcStream, true);
        this.normalizationManifest = null;
        parse();
    }

    private void calculateStatistics() {
        calculateRandTheta();
    }

    private void calculateRandTheta() {
        thetaIlmn = new float[normalizedXIntensities.length];
        RIlmn = new float[normalizedXIntensities.length];
        for (int i = 0; i < normalizedXIntensities.length; i++) {
            double x = normalizedXIntensities[i];
            double y = normalizedYIntensities[i];
            thetaIlmn[i] = (float) (2 * (Math.atan(y / x) / Math.PI));
            RIlmn[i] = (float) (x + y);
        }
    }

    /**
     * Main parsing method.
     *
     * @throws IOException thrown when there is a problem reading the stream.
     */
    private void parse() throws IOException {

        stream.mark(0);

        try {
            final byte[] curIdentifier = new byte[IDENTIFIER_LENGTH];
            for (int i = 0; i < curIdentifier.length; i++) {
                curIdentifier[i] = stream.readByte();
            }

            final String identifier = new String(curIdentifier);
            setIdentifier(identifier);
            if (!identifier.equals(GTC_IDENTIFIER)) {
                throw new PicardException("Invalid identifier '" + identifier + "' for GTC file");
            }
            setFileVersion(stream.readByte());
            setNumberOfEntries(Integer.reverseBytes(stream.readInt()));

            //parse the tables of contents
            for (InfiniumFileTOC toc : getTableOfContents()) {
                stream.reset();
                readData(stream, toc);
            }
            //older versions don't have extended sample info so we need to infer this
            if (numCalls == 0) {
                numCalls = aaCalls + abCalls + bbCalls;
            }

        } finally {
            stream.close();
        }

        if ((normalizationManifest != null) && (normalizationManifest.getNormIds() != null)) {
            normalizeIntensities();
        }
    }

    private void normalizeIntensities() {
        normalizedXIntensities = new float[numberOfSnps];
        normalizedYIntensities = new float[numberOfSnps];

        final int[] normIds = normalizationManifest.getNormIds();
        for (int i = 0; i < rawXIntensities.length; i++) {
            final int rawX = rawXIntensities[i];
            final int rawY = rawYIntensities[i];

            final int normId;
            int normIndex = -1;
            if ((normIds != null) && (normIds.length > i)) {
                normId = normIds[i];
                normIndex = getAllNormIndex(normId);
            }

            if (normIndex != -1) {
                final InfiniumTransformation xform = normalizationTransformations[normIndex];
                final float tempX = rawX - xform.getOffsetX();
                final float tempY = rawY - xform.getOffsetY();
                final float theta = xform.getTheta();
                double tempX2 = Math.cos(theta) * tempX + Math.sin(theta) * tempY;
                double tempY2 = -Math.sin(theta) * tempX + Math.cos(theta) * tempY;
                double tempX3 = tempX2 - xform.getShear() * tempY2;

                if (tempX3 < 0) {
                    tempX3 = 0;
                }
                if (tempY2 < 0) {
                    tempY2 = 0;
                }

                normalizedXIntensities[i] = (float) (tempX3 / xform.getScaleX());
                normalizedYIntensities[i] = (float) (tempY2 / xform.getScaleY());
            } else {
                normalizedXIntensities[i] = rawXIntensities[i];
                normalizedYIntensities[i] = rawYIntensities[i];
            }
        }
    }

    private int getAllNormIndex(final int normId) {
        int index = 0;

        for (int currentNormId : normalizationManifest.getAllNormIds()) {
            if (currentNormId == normId) {
                return index;
            }
            index++;
        }
        return -1;
    }

    /**
     * Reads the table of contents data from the input stream.
     *
     * @param stream The stream to be parsed.
     * @param toc    The table of contents record to be parsed.
     * @throws IOException thrown when there is a problem reading the stream.
     */
    private void readData(final DataInputStream stream, final InfiniumFileTOC toc) throws IOException {
        switch (toc.getTableOfContentsId()) {
            case NUM_SNPS:
                numberOfSnps = toc.getOffset();
                break;
            case PLOIDY:
                ploidy = toc.getOffset();
                break;
            case PLOIDY_TYPE:
                ploidyType = toc.getOffset();
                break;
            case SAMPLE_NAME:
                sampleName = parseString(toc);
                break;
            case SAMPLE_PLATE:
                samplePlate = parseString(toc);
                break;
            case SAMPLE_WELL:
                sampleWell = parseString(toc);
                break;
            case CLUSTER_FILE:
                clusterFile = parseString(toc);
                break;
            case SNP_MANIFEST:
                snpManifest = parseString(toc);
                break;
            case IMAGING_DATE:
                imagingDate = parseString(toc);
                break;
            case AUTOCALL_DATE:
                autoCallDate = parseString(toc);
                break;
            case AUTOCALL_VERSION:
                autoCallVersion = parseString(toc);
                break;
            case TRANSFORMATIONS:
                parseTransformations(toc);
                break;
            case RAW_CONTROL_X_INTENSITIES:
                rawControlXIntensities = parseUnsignedShortArray(toc);
                break;
            case RAW_CONTROL_Y_INTENSITIES:
                rawControlYIntensities = parseUnsignedShortArray(toc);
                break;
            case RAW_X_INTENSITIES:
                rawXIntensities = parseUnsignedShortArray(toc);
                break;
            case RAW_Y_INTESITIES:
                rawYIntensities = parseUnsignedShortArray(toc);
                break;
            case GENOTYPES:
                genotypeBytes = parseGenotypes(toc);
                break;
            case BASE_CALLS:
                baseCalls = parseBaseCalls(toc);
                break;
            case GENOTYPE_SCORES:
                genotypeScores = parseFloatArray(toc);
                break;
            case SCANNER_INFO:
                parseScannerInfo(toc);
                break;
            case CALL_RATE:
                callRate = parseFloat(toc);
                break;
            case GENDER:
                stream.skipBytes(toc.getOffset());
                gender = String.valueOf((char) stream.read());
                break;
            case LOG_R_DEV:
                logRDev = parseFloat(toc);
                break;
            case P_10_GC:
                p10GC = parseFloat(toc);
                break;
            case DX:
                dx = parseInt(toc);
                break;
            case EXTENDED_SAMPLE_DATA:
                parseExtendedSampleData(toc);
                break;
            case B_ALLELE_FREQS:
                bAlleleFreqs = parseFloatArray(toc);
                break;
            case LOG_R_RATIOS:
                logRRatios = parseFloatArray(toc);
                break;
            case INTENSITY_X_PERCENTILES:
                redIntensityPercentiles = new IntensityPercentiles(parseShort(toc), readShort(), readShort());
                break;
            case INTENSITY_Y_PERCENTILES:
                greenIntensityPercentiles = new IntensityPercentiles(parseShort(toc), readShort(), readShort());
                break;
            case SENTRIX_ID:
                sentrixBarcode = parseString(toc);
                break;
            default:
                // throw new MPGException(new StringBuilder().append("Unknown GTC TOC id: ").append(toc.getTableOfContentsId()).toString());
        }
    }

    private void parseExtendedSampleData(final InfiniumFileTOC toc) throws IOException {
        p50GC = parseFloat(toc);
        numCalls = Integer.reverseBytes(stream.readInt());
        numNoCalls = Integer.reverseBytes(stream.readInt());
        numIntensityOnly = Integer.reverseBytes(stream.readInt());
    }

    /**
     * Utility method for parsing the scanner information.
     *
     * @param toc The table of contents record to read the scanner information from.
     * @throws IOException is thrown when there is a problem reading the stream.
     */
    private void parseScannerInfo(final InfiniumFileTOC toc) throws IOException {
        scannerName = parseString(toc);
        pmtGreen = Integer.reverseBytes(stream.readInt());
        pmtRed = Integer.reverseBytes(stream.readInt());
        scannerVersion = parseString();
        imagingUser = parseString();
    }

    /**
     * Utility method for parsing out the genotypes (NC for no call)
     *
     * @param toc The table of contents record for parsing the genotypes.
     * @return A byte array containing all of the genotype values.
     * @throws IOException is thrown when there is a problem reading the stream.
     */
    private byte[] parseGenotypes(final InfiniumFileTOC toc) throws IOException {
        final byte[] genotypeBytes = parseByteArray(toc);

        for (byte genotypeByte : genotypeBytes) {
            switch (genotypeByte) {
                case NO_CALL: {
                    break;
                }
                case AA_CALL: {
                    aaCalls++;
                    break;
                }
                case AB_CALL: {
                    abCalls++;
                    break;
                }
                case BB_CALL: {
                    bbCalls++;
                    break;
                }
            }
        }
        return genotypeBytes;
    }

    /**
     * Utility method for parsing out the basecalls. (- for no call)
     *
     * @param toc The table of contents record for parsing the basecalls.
     * @return A string array containing all of the basecall values.
     * @throws IOException is thrown when there is a problem reading the stream.
     */
    private byte[][] parseBaseCalls(final InfiniumFileTOC toc) throws IOException {
        stream.skipBytes(toc.getOffset());
        final int arrayLen = Integer.reverseBytes(stream.readInt());
        final byte[][] curBaseCalls = new byte[arrayLen][2];
        for (int i = 0; i < arrayLen; i++) {
            byte[] baseCallBytes = curBaseCalls[i];
            for (int j = 0; j < baseCallBytes.length; j++) {
                baseCallBytes[j] = stream.readByte();
                if (baseCallBytes[j] == 0) {
                    baseCallBytes[j] = NO_CALL_CHAR;
                }
            }
        }

        return curBaseCalls;
    }

    /**
     * Utility method for parsing the normalization transformation.
     *
     * @param toc The table of contents record to parse the transformation from.
     * @throws IOException is thrown when there is a problem reading the stream.
     */
    private void parseTransformations(final InfiniumFileTOC toc) throws IOException {
        stream.skipBytes(toc.getOffset());
        final int arrayLen = Integer.reverseBytes(stream.readInt());
        final InfiniumTransformation[] transformations = new InfiniumTransformation[arrayLen];
        for (int i = 0; i < transformations.length; i++) {
            InfiniumTransformation curTransformation = new InfiniumTransformation();
            curTransformation.setVersion(Integer.reverseBytes(stream.readInt()));
            curTransformation.setOffsetX(parseFloat());
            curTransformation.setOffsetY(parseFloat());
            curTransformation.setScaleX(parseFloat());
            curTransformation.setScaleY(parseFloat());
            curTransformation.setShear(parseFloat());
            curTransformation.setTheta(parseFloat());
            curTransformation.setReserved1(parseFloat());
            curTransformation.setReserved2(parseFloat());
            curTransformation.setReserved3(parseFloat());
            curTransformation.setReserved4(parseFloat());
            curTransformation.setReserved5(parseFloat());
            curTransformation.setReserved6(parseFloat());
            transformations[i] = curTransformation;
        }
        normalizationTransformations = transformations;
    }

    private void normalizeAndCalculateStatistics() {
        if (normalizationManifest.getNormIds() != null) {
            normalizeIntensities();
        }
        calculateStatistics();
    }

    public InfiniumGTCRecord getRecord(int index) {
        return new InfiniumGTCRecord(
                rawXIntensities[index],
                rawYIntensities[index],
                genotypeBytes[index],
                genotypeScores[index],
                normalizedXIntensities[index],
                normalizedYIntensities[index],
                RIlmn[index],
                thetaIlmn[index],
                bAlleleFreqs[index],
                logRRatios[index]
        );
    }

    public double getHetPercent() {
        return (double) abCalls / (double) numCalls;
    }

    public String getSampleName() {
        return sampleName;
    }

    public String getSamplePlate() {
        return samplePlate;
    }

    public String getSampleWell() {
        return sampleWell;
    }

    public String getClusterFile() {
        return clusterFile;
    }

    public String getSnpManifest() {
        return snpManifest;
    }

    public String getImagingDate() {
        return imagingDate;
    }

    public String getAutoCallDate() {
        return autoCallDate;
    }

    public String getAutoCallVersion() {
        return autoCallVersion;
    }

    public int[] getRawControlXIntensities() {
        return rawControlXIntensities;
    }

    public int[] getRawControlYIntensities() {
        return rawControlYIntensities;
    }

    public String getScannerName() {
        return scannerName;
    }

    public int getPmtGreen() {
        return pmtGreen;
    }

    public int getPmtRed() {
        return pmtRed;
    }

    public String getScannerVersion() {
        return scannerVersion;
    }

    public String getImagingUser() {
        return imagingUser;
    }

    public double getCallRate() {
        return callRate;
    }

    public String getGender() {
        return gender;
    }

    public int getNumberOfSnps() {
        return numberOfSnps;
    }

    public int getNumCalls() {
        return numCalls;
    }

    public int getNumNoCalls() {
        return numNoCalls;
    }

    public int getRawControlXIntensity(int index) {
        return rawControlXIntensities[index];
    }

    public int getRawControlYIntensity(int index) {
        return rawControlYIntensities[index];
    }

    public int getPloidy() {
        return ploidy;
    }

    public int getPloidyType() {
        return ploidyType;
    }

    public int getP05Red() { return redIntensityPercentiles.p05; }

    public int getP50Red() { return redIntensityPercentiles.p50; }

    public int getP95Red() { return redIntensityPercentiles.p95; }

    public int getP05Green() { return greenIntensityPercentiles.p05; }

    public int getP50Green() { return greenIntensityPercentiles.p50; }

    public int getP95Green() { return greenIntensityPercentiles.p95; }

    public float getLogRDev() {
        return logRDev;
    }

    public float getP10GC() {
        return p10GC;
    }

    public float getP50GC() {
        return p50GC;
    }

    public int getNumIntensityOnly() {
        return numIntensityOnly;
    }

    public long getAaCalls() {
        return aaCalls;
    }

    public long getBbCalls() {
        return bbCalls;
    }

    public String getSentrixBarcode() {
        return sentrixBarcode;
    }

    public int getDx() {
        return dx;
    }

    public byte[][] getBaseCalls() {
        return baseCalls;
    }

    public int getAbCalls() {
        return abCalls;
    }

    /**
     * A class to store Illumina Intensity Percentiles
     */
    private class IntensityPercentiles {
        private final int p05;
        private final int p50;
        private final int p95;

        IntensityPercentiles(int p05, int p50, int p95) {
            this.p05 = p05;
            this.p50 = p50;
            this.p95 = p95;
        }
    }

}