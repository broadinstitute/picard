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

import org.apache.commons.lang.ArrayUtils;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * A class to represent an 'Extended' Illumina Manifest file.
 *
 * An Extended Illumina Manifest extends a 'standard' Illumina Manifest by adding seven new fields (columns) to the manifest
 * These are currently specific to the reference NCBI Build 37 / HG38.
 * The columns are:
 *   'build37Chr' - the chromosome of the manifest entry, on build 37
 *   'build37Pos' - the position of the manifest entry, on build 37
 *   'build37RefAllele' - the reference allele of the manifest entry, on build 37
 *   'build37AlleleA' - allele A of the manifest entry, on build 37
 *   'build37AlleleB' - allele B of the manifest entry, on build 37
 *   'build37Rsid' - The rsid the manifest entry, on build 37
 *   'build37Flag' - A flag describing the validation status of the manifest entry
 *
 * Like the class IlluminaManifest which this class extends, this class reads the extended manifest header, stores the contents,
 * and then provides an iterator to allow access to the ExtendedIlluminaManifestRecords
 * (currently this only supports iterating over the assay records).
 */
public class Build37ExtendedIlluminaManifest extends IlluminaManifest {

    private static final String BUILD37_CHR_HEADER_NAME = "build37Chr";
    private static final String BUILD37_POS_HEADER_NAME = "build37Pos";
    private static final String BUILD37_REF_ALLELE_HEADER_NAME = "build37RefAllele";
    private static final String BUILD37_ALLELE_A_HEADER_NAME = "build37AlleleA";
    private static final String BUILD37_ALLELE_B_HEADER_NAME = "build37AlleleB";
    private static final String BUILD37_RSID_NAME = "build37Rsid";
    private static final String BUILD37_FLAG_HEADER_NAME = "build37Flag";

    static final String[] EXTENDED_MANIFEST_HEADERS = {
            BUILD37_CHR_HEADER_NAME,
            BUILD37_POS_HEADER_NAME,
            BUILD37_REF_ALLELE_HEADER_NAME,
            BUILD37_ALLELE_A_HEADER_NAME,
            BUILD37_ALLELE_B_HEADER_NAME,
            BUILD37_RSID_NAME,
            BUILD37_FLAG_HEADER_NAME
    };

    static final String EXTENDED_MANIFEST_VERSION_HEADER_NAME = "CreateExtendedIlluminaManifest.version";
    static final String EXTENDED_MANIFEST_TARGET_BUILD_HEADER_NAME = "Target Build";
    static final String EXTENDED_MANIFEST_TARGET_REFERENCE_HEADER_NAME = "Target Reference File";
    static final String EXTENDED_MANIFEST_CLUSTER_FILE_HEADER_NAME = "Cluster File";
    static final String EXTENDED_MANIFEST_DBSNP_FILE_HEADER_NAME = "dbSNP File";
    static final String EXTENDED_MANIFEST_SUPPORTED_BUILD_HEADER_NAME = "Supported Build";
    static final String EXTENDED_MANIFEST_SUPPORTED_REFERENCE_HEADER_NAME = "Supported Reference File";
    static final String EXTENDED_MANIFEST_SUPPORTED_CHAIN_FILE_HEADER_NAME = "Supported Chain File";

    @Override
    public String[] getAllPossibleHeaderNames() {
        return (String[])ArrayUtils.addAll(super.getAllPossibleHeaderNames(), EXTENDED_MANIFEST_HEADERS);
    }

    public Build37ExtendedIlluminaManifest(final File manifestFile) throws IOException {
        super(manifestFile);
    }

    public Iterator<Build37ExtendedIlluminaManifestRecord> extendedIterator() {

        return new Iterator<Build37ExtendedIlluminaManifestRecord>() {
            private int assayCount = 0;
            private int numAssays = getNumAssays();

            public boolean hasNext() {
                return (assayCount < numAssays) && (manifestFileParser.hasNext()) ;
            }

            public Build37ExtendedIlluminaManifestRecord next() {
                assayCount++;
                return new Build37ExtendedIlluminaManifestRecord(getAssayHeaderNameToIndex(), manifestFileParser.next(), assayCount - 1);
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    public String getExtendedManifestVersion() {
        String version = "?";
        for (String[] headerLine: getHeaderContents()) {
            if (headerLine[0].equals(EXTENDED_MANIFEST_VERSION_HEADER_NAME)) {
                version = headerLine[1];
            }
        }
        return version;
    }
}
