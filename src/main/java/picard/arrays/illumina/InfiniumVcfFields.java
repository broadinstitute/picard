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

import htsjdk.samtools.util.Iso8601Date;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.text.ParseException;
import java.text.SimpleDateFormat;

/**
 * A class to store fields that are specific to a VCF generated from an Illumina GTC file.
 */
public class InfiniumVcfFields {
    // Header Fields
    public static final String ARRAY_TYPE = "arrayType";

    public static final String EXTENDED_ILLUMINA_MANIFEST_VERSION = "extendedIlluminaManifestVersion";
    public static final String CHIP_WELL_BARCODE = "chipWellBarcode";
    public static final String ANALYSIS_VERSION_NUMBER = "analysisVersionNumber";
    public static final String SAMPLE_ALIAS = "sampleAlias";

    public static final String EXPECTED_GENDER = "expectedGender";
    public static final String FINGERPRINT_GENDER = "fingerprintGender";
    public static final String AUTOCALL_GENDER = "autocallGender";
    public static final String AUTOCALL_DATE = "autocallDate";
    public static final String IMAGING_DATE = "imagingDate";
    public static final String CLUSTER_FILE = "clusterFile";
    public static final String MANIFEST_FILE = "manifestFile";
    public static final String EXTENDED_ILLUMINA_MANIFEST_FILE = "extendedManifestFile";
    public static final String AUTOCALL_VERSION = "autocallVersion";
    public static final String ZCALL_VERSION = "zcallVersion";
    public static final String ZCALL_THRESHOLDS = "zcallThresholds";
    public static final String P_95_RED = "p95Red";
    public static final String P_95_GREEN = "p95Green";
    public static final String SCANNER_NAME = "scannerName";

    //FORMAT Fields
    public static final String X = "X";
    public static final String Y = "Y";
    public static final String NORMX = "NORMX";
    public static final String NORMY = "NORMY";
    public static final String R = "R";
    public static final String THETA = "THETA";
    public static final String BAF = "BAF";
    public static final String LRR = "LRR";
    public static final String IGC = "IGC";
    public static final String GTA = "GTA";
    public static final String GTZ = "GTZ";

    //INFO Fields
    public static final String ALLELE_A = "ALLELE_A";
    public static final String ALLELE_B = "ALLELE_B";
    public static final String ILLUMINA_STRAND = "ILLUMINA_STRAND";
    public static final String PROBE_A = "PROBE_A";
    public static final String PROBE_B = "PROBE_B";
    public static final String BEADSET_ID = "BEADSET_ID";
    public static final String ILLUMINA_CHR = "ILLUMINA_CHR";
    public static final String ILLUMINA_POS = "ILLUMINA_POS";
    public static final String ILLUMINA_BUILD = "ILLUMINA_BUILD";
    public static final String SOURCE = "SOURCE";
    public static final String GC_SCORE = "GC_SCORE";
    public static final String[] N = new String[GENOTYPE_VALUES.values().length];
    public static final String[] MEAN_R = new String[GENOTYPE_VALUES.values().length];
    public static final String[] DEV_R = new String[GENOTYPE_VALUES.values().length];
    public static final String[] MEAN_THETA = new String[GENOTYPE_VALUES.values().length];
    public static final String[] DEV_THETA = new String[GENOTYPE_VALUES.values().length];
    public static final String[] MEAN_X = new String[GENOTYPE_VALUES.values().length];
    public static final String[] DEV_X = new String[GENOTYPE_VALUES.values().length];
    public static final String[] MEAN_Y = new String[GENOTYPE_VALUES.values().length];
    public static final String[] DEV_Y = new String[GENOTYPE_VALUES.values().length];
    public static final String RS_ID = "refSNP";
    public static final String ZTHRESH_X = "zthresh_X";
    public static final String ZTHRESH_Y = "zthresh_Y";

    static {
        for (GENOTYPE_VALUES gtValue : GENOTYPE_VALUES.values()) {
            N[gtValue.ordinal()] = "N_" + gtValue.name();
            MEAN_R[gtValue.ordinal()] = "meanR_" + gtValue.name();
            DEV_R[gtValue.ordinal()] = "devR_" + gtValue.name();
            MEAN_THETA[gtValue.ordinal()] = "meanTHETA_" + gtValue.name();
            DEV_THETA[gtValue.ordinal()] = "devTHETA_" + gtValue.name();
            MEAN_X[gtValue.ordinal()] = "meanX_" + gtValue.name();
            DEV_X[gtValue.ordinal()] = "devX_" + gtValue.name();
            MEAN_Y[gtValue.ordinal()] = "meanY_" + gtValue.name();
            DEV_Y[gtValue.ordinal()] = "devY_" + gtValue.name();
        }
    }

    //FILTER Fields
    public static final String TRIALLELIC = "TRIALLELIC";
    public static final String DUPE = "DUPE";
    public static final String FAIL_REF = "FAIL_REF";
    public static final String ZEROED_OUT_ASSAY = "ZEROED_OUT_ASSAY";

    public enum GENOTYPE_VALUES {
        AA(0),
        AB(1),
        BB(2);

        int v; //value is the index of the value for (AA, AB, BB) as stored in Illumina files.

        GENOTYPE_VALUES(final int v) {
            this.v = v;
        }
    }
    final static int NUM_GENOTYPE_VALUES = GENOTYPE_VALUES.values().length;

    static String getValueFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String keyName) {
        VCFHeaderLine otherHeaderLine = vcfHeader.getOtherHeaderLine(keyName);
        if (otherHeaderLine != null) {
            return otherHeaderLine.getValue();
        } else {
            throw new IllegalArgumentException("Input VCF file is missing header line of type '" + keyName + "'");
        }
    }

    static Iso8601Date getDateFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String keyName, final SimpleDateFormat dateformat) {
        String dateString = InfiniumVcfFields.getValueFromVcfOtherHeaderLine(vcfHeader, keyName);
        try {
            return new Iso8601Date(dateformat.parse(dateString));
        } catch (ParseException pe) {
            throw new IllegalArgumentException("Unrecognized date for '" + keyName + "' in VCF header (" + dateString + ")");
        }
    }

    static Integer getIntegerFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String keyName) {
        VCFHeaderLine otherHeaderLine = vcfHeader.getOtherHeaderLine(keyName);
        if (otherHeaderLine != null) {
            return Integer.valueOf(otherHeaderLine.getValue());
        } else {
            throw new IllegalArgumentException("Input VCF file is missing header line of type '" + keyName + "'");
        }
    }

    static String getOptionalValueFromVcfOtherHeaderLine(final VCFHeader vcfHeader, final String keyName) {
        VCFHeaderLine otherHeaderLine = vcfHeader.getOtherHeaderLine(keyName);
        if (otherHeaderLine != null) {
            return otherHeaderLine.getValue();
        } else {
            return null;
        }
    }
}
