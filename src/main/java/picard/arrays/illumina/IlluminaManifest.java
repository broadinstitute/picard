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

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import picard.PicardException;
import picard.util.CsvInputParser;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A class to represent an Illumina Manifest file.
 *
 * Reads the header, stores the contents, and then provides an iterator to allow
 * access to the IlluminaManifestRecords (currently this only supports iterating over the assay records).
 */
public class IlluminaManifest  {
    public static final String HG17 = "HG17";
    public static final String HG18 = "HG18";
    public static final String HG19 = "HG19";

    public static final String NCBI_35 = "35";
    public static final String NCBI_36 = "36";
    public static final String NCBI_37 = "37";

    public static final BidiMap HG_TO_NCBI = new DualHashBidiMap();

    static {
        HG_TO_NCBI.put(HG17, NCBI_35);
        HG_TO_NCBI.put(HG18, NCBI_36);
        HG_TO_NCBI.put(HG19, NCBI_37);
    }

    private static final String[] ALLELES_LIST = {"A", "C", "G", "T", "I", "D"};
    public static final HashSet<String> VALID_ALLELES = new HashSet<>(Arrays.asList(ALLELES_LIST));

    public static final String ILLUMINA_ID_HEADER_NAME = "IlmnID";
    public static final String NAME_HEADER_NAME = "Name";
    public static final String ILLUMINA_STRAND_HEADER_NAME = "IlmnStrand";
    public static final String SNP_HEADER_NAME = "SNP";
    public static final String ADDRESS_A_ID_HEADER_NAME = "AddressA_ID";
    public static final String ALLELE_A_PROBE_SEQ_HEADER_NAME = "AlleleA_ProbeSeq";
    public static final String ADDRESS_B_ID_HEADER_NAME = "AddressB_ID";
    public static final String ALLELE_B_PROBE_SEQ_HEADER_NAME = "AlleleB_ProbeSeq";
    public static final String GENOME_BUILD_HEADER_NAME = "GenomeBuild";
    public static final String CHROMOSOME_HEADER_NAME = "Chr";
    public static final String MAP_INFO_HEADER_NAME = "MapInfo";
    public static final String PLOIDY_HEADER_NAME = "Ploidy";
    public static final String SPECIES_HEADER_NAME = "Species";
    public static final String SOURCE_HEADER_NAME = "Source";
    public static final String SOURCE_VERSION_HEADER_NAME = "SourceVersion";
    public static final String SOURCE_STRAND_HEADER_NAME = "SourceStrand";
    public static final String SOURCE_SEQ_HEADER_NAME = "SourceSeq";
    public static final String TOP_GENOMIC_SEQ_HEADER_NAME = "TopGenomicSeq";
    public static final String BEAD_SET_ID_HEADER_NAME = "BeadSetID";
    public static final String EXP_CLUSTERS_HEADER_NAME = "Exp_Clusters";
    public static final String REF_STRAND_HEADER_NAME = "RefStrand";
    public static final String INTENSITY_ONLY_HEADER_NAME = "Intensity_Only";

    public static final String[] MANIFEST_FILE_HEADER_NAMES = {
            ILLUMINA_ID_HEADER_NAME,
            NAME_HEADER_NAME,
            ILLUMINA_STRAND_HEADER_NAME,
            SNP_HEADER_NAME,
            ADDRESS_A_ID_HEADER_NAME,
            ALLELE_A_PROBE_SEQ_HEADER_NAME,
            ADDRESS_B_ID_HEADER_NAME,
            ALLELE_B_PROBE_SEQ_HEADER_NAME,
            GENOME_BUILD_HEADER_NAME,
            CHROMOSOME_HEADER_NAME,
            MAP_INFO_HEADER_NAME,
            PLOIDY_HEADER_NAME,
            SPECIES_HEADER_NAME,
            SOURCE_HEADER_NAME,
            SOURCE_VERSION_HEADER_NAME,
            SOURCE_STRAND_HEADER_NAME,
            SOURCE_SEQ_HEADER_NAME,
            TOP_GENOMIC_SEQ_HEADER_NAME,
            BEAD_SET_ID_HEADER_NAME,
            EXP_CLUSTERS_HEADER_NAME,
            REF_STRAND_HEADER_NAME,
            INTENSITY_ONLY_HEADER_NAME
    };
    public static final HashSet<String> HEADER_NAMES = new HashSet<>(Arrays.asList(MANIFEST_FILE_HEADER_NAMES));

    private final Log log = Log.getInstance(IlluminaManifest.class);

    private final File manifestFile;
    private final List<String[]> headerContents = new ArrayList<>();
    protected CsvInputParser manifestFileParser;

    private String descriptorFileName;
    private String assayFormat;
    private String dateManufactured;
    private int lociCount;
    private int numAssays;

    private String[] assayHeaderNames;
    private Map<String, Integer> assayHeaderNameToIndex;

    public String[] getAllPossibleHeaderNames() {
        return MANIFEST_FILE_HEADER_NAMES;
    }

    public String[] getAssayHeaderNames() {
        return assayHeaderNames;
    }

    public Map<String, Integer> getAssayHeaderNameToIndex() {
        return assayHeaderNameToIndex;
    }

    public IlluminaManifest(final File manifestFile) throws IOException {
        IOUtil.assertFileIsReadable(manifestFile);
        this.manifestFile = manifestFile;
        readHeader();
        setNumAssays(lociCount);
        init();
    }

    public IlluminaManifest(final File manifestFile, final int numAssays) throws IOException {
        IOUtil.assertFileIsReadable(manifestFile);
        this.manifestFile = manifestFile;
        readHeader();
        setNumAssays(numAssays);
        init();
    }

    private void init() throws IOException {
        // Reopen the file, and advance to the [Assay] Section
        final FileInputStream fileInputStream = new FileInputStream(getManifestFile());
        advanceDataStream(fileInputStream, "[Assay]");

        // Read the next line from the file - the header for the Assay records.  Then validate it.
        final String headerLine = readLineFromStream(fileInputStream);            // Read the header.
        validateManifestRecordHeader(headerLine);

        // Begin the parser just after the header.
        this.manifestFileParser = new CsvInputParser(false, fileInputStream);
    }

    public Iterator<IlluminaManifestRecord> iterator() {

        return new Iterator<IlluminaManifestRecord>() {
            private int assayCount = 0;
            private int numAssays = getNumAssays();
            private boolean firstCall = true;
            private String[] nextLine = null;

            public boolean hasNext() {
                if (firstCall) {
                    firstCall = false;
                    nextLine = manifestFileParser.next();
                }
                return ((assayCount < numAssays) && (nextLine != null) && (IlluminaManifestRecord.isValidManifestRecordLine(assayHeaderNames.length, nextLine)));
            }

            public IlluminaManifestRecord next() {
                String[] currentLine = nextLine;
                assayCount++;
                nextLine = (manifestFileParser.hasNext()) ? manifestFileParser.next() : null;       // Pre-fetch the next line
                return new IlluminaManifestRecord(assayHeaderNameToIndex, currentLine, assayCount - 1);
            }
        };
    }

    private void readHeader() throws IOException {
        final FileInputStream fileInputStream = new FileInputStream(getManifestFile());

        // Now that the inputStream is advanced to the line just beyond the sectionName, create a CsvInputParser for that section
        this.manifestFileParser = new CsvInputParser(false, fileInputStream);

        boolean inHeader = true;
        while (manifestFileParser.hasNext()) {
            final String[] row = manifestFileParser.next();
            headerContents.add(row);
            final String tagName = row[0].trim();         // Remove trailing whitespace.
            switch (tagName) {
                case "Descriptor File Name":
                    setDescriptorFileName(row[1]);
                    break;
                case "Assay Format":
                    setAssayFormat(row[1]);
                    break;
                case "Date Manufactured":
                    setDateManufactured(row[1]);
                    break;
                case "Loci Count":
                    setLociCount(new Integer(row[1]));
                    inHeader = false;           // Need to end the header now or parser will blow up on next line
                                                // Note that this means that the header ends at this field
                    break;
            }
            if (!inHeader) {
                break;
            }
        }
        fileInputStream.close();
    }

    private void validateManifestRecordHeader(final String line) throws PicardException {
        final String[] columns = line.trim().split(",");
        Map<String, Integer> columnNameToIndex = new HashMap<>();
        int index = 0;
        Set<String> validHeaderNames = new HashSet<>(Arrays.asList(getAllPossibleHeaderNames()));
        for (String columnName : columns) {
            if (!validHeaderNames.contains(columnName)) {
                throw new PicardException("Unrecognized Column '" + columnName + "' in Manifest file: " + manifestFile);
            }
            columnNameToIndex.put(columnName, index++);
        }
        if (!columnNameToIndex.containsKey(REF_STRAND_HEADER_NAME)) {
            // Some Illumina manifests do not have ref_strand defined.  We will use the illumina strand instead.
            log.warn("Illumina Manifest does not contain '" + REF_STRAND_HEADER_NAME + "' - we will use '" + ILLUMINA_STRAND_HEADER_NAME + "'");
        }
        assayHeaderNames = columns;
        assayHeaderNameToIndex = columnNameToIndex;
    }

    public String[] getManifestFileHeaderNames() {
        return assayHeaderNames;
    }

    /**
     * Advance the stream to the line after the sectionName.
     * Not using buffers since they read too much of the stream in and confounds AbstractInputParser later
     * (i.e. the buffer has already pulled in too much from the stream)
     *
     * Need to do this because AbstractInputParser cannot handle files with sections of differing lengths of tab-delimited fields
     * i.e. the gtc.txt file with a block of one tab data, then a block of eight tab data
     *
     * @param inputStream
     * @param sectionName
     */
    private void advanceDataStream(final FileInputStream inputStream, final String sectionName) throws IOException {
        while (true) {
            final String line = readLineFromStream(inputStream).trim();     // trim to remove trailing '/r'
            if (line.contains(sectionName)) {
                break;
            }
        }
    }

    /**
     * Read a single line of text from the stream.
     *
     * @param inputStream the input stream
     * @return A line of text
     * @throws IOException on error
     */
    private String readLineFromStream(final InputStream inputStream) throws IOException {
        final StringBuilder sb = new StringBuilder();
        int ch;
        while ((ch = inputStream.read()) > -1) {
            if (ch == '\n') {
                break;
            }
            sb.append((char) ch);
        }
        if (ch == -1) {
            throw new PicardException("Unexpected end of file");
        }
        return sb.toString();
    }



    public File getManifestFile() {
        return manifestFile;
    }

    public CsvInputParser getManifestFileParser() {
        return manifestFileParser;
    }

    public void setManifestFileParser(final CsvInputParser manifestFileParser) {
        this.manifestFileParser = manifestFileParser;
    }

    public String getDescriptorFileName() {
        return descriptorFileName;
    }

    public void setDescriptorFileName(final String descriptorFileName) {
        this.descriptorFileName = descriptorFileName;
    }

    public String getAssayFormat() {
        return assayFormat;
    }

    public void setAssayFormat(final String assayFormat) {
        this.assayFormat = assayFormat;
    }

    public String getDateManufactured() {
        return dateManufactured;
    }

    public void setDateManufactured(final String dateManufactured) {
        this.dateManufactured = dateManufactured;
    }

    public void setLociCount(final int lociCount) {
        this.lociCount = lociCount;
    }

    public int getLociCount() {
        return lociCount;
    }

    public int getNumAssays() {
        return numAssays;
    }

    public void setNumAssays(final int numAssays) {
        this.numAssays = numAssays;
    }

    public List<String[]> getHeaderContents() { return headerContents; }
}
