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

import htsjdk.tribble.annotation.Strand;
import org.apache.commons.lang.StringUtils;
import picard.PicardException;

import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A class to represent a record (line) from an Illumina Manifest [Assay] entry
 */
public class IlluminaManifestRecord {
    private final String[] originalLine;        // A copy of the line from the original file
    private final int index;

    private final String ilmnId;
    private final String name;
    private IlluminaStrand ilmnStrand;

    private final String snp;
    private final boolean isIndel;              // Not part of the file...
    private final boolean isAmbiguous;          // Not part of the file

    private final String addressAId;
    private final String alleleAProbeSeq;
    private final String addressBId;
    private final String alleleBProbeSeq;

    private final String genomeBuild;
    private final String majorGenomeBuild;      // Not part of the file... (35, 36, 37)
    private String hgGenomeBuild;               // Not part of the file... (HG17, HG18, HG19)

    private final String chr;
    private final int position;                       // annotated as 'MapInfo'

    private final String ploidy;
    private final String species;
    private final String source;
    private final String sourceVersion;
    private IlluminaStrand sourceStrand;
    private final String sourceSeq;
    private final String topGenomicSeq;
    private final int beadSetId;
    private final String expClusters;
    private Strand refStrand;
    private final boolean intensityOnly;

    private static final Set<String> indels = Stream.of("[D/I]", "[I/D]").collect(Collectors.toSet());
    private static final Set<String> ambiguousSnps = Stream.of("[A/T]", "[T/A]", "[G/C]", "[C/G]").collect(Collectors.toSet());

    protected enum IlluminaStrand {
        PLUS,
        MINUS,
        TOP,
        BOT,
        NONE
    }

    public static final String ILLUMINA_FLAGGED_BAD_CHR = "0";

    IlluminaManifestRecord(final Map<String, Integer> columnNameToIndex, final String[] line, final int index) {
        this.originalLine = line;
        this.index = index;

        ilmnId = getColumnValue(columnNameToIndex, IlluminaManifest.ILLUMINA_ID_HEADER_NAME);
        name = getColumnValue(columnNameToIndex, IlluminaManifest.NAME_HEADER_NAME);
        ilmnStrand = getIlmnStrandFromManifest(columnNameToIndex);
        snp = getColumnValue(columnNameToIndex, IlluminaManifest.SNP_HEADER_NAME).toUpperCase();         // This is of the form [A/T] or [D/I].
        isIndel = indels.contains(getSnp());
        isAmbiguous = ambiguousSnps.contains(getSnp());
        addressAId = getColumnValue(columnNameToIndex, IlluminaManifest.ADDRESS_A_ID_HEADER_NAME);
        alleleAProbeSeq = getColumnValue(columnNameToIndex, IlluminaManifest.ALLELE_A_PROBE_SEQ_HEADER_NAME);
        addressBId = getColumnValue(columnNameToIndex, IlluminaManifest.ADDRESS_B_ID_HEADER_NAME);
        alleleBProbeSeq = getColumnValue(columnNameToIndex, IlluminaManifest.ALLELE_B_PROBE_SEQ_HEADER_NAME);

        genomeBuild = getColumnValue(columnNameToIndex, IlluminaManifest.GENOME_BUILD_HEADER_NAME);   // This is (usually...) the NCBI Genome Build (36, 36.1, 37, 37.2)
        chr = getColumnValue(columnNameToIndex, IlluminaManifest.CHROMOSOME_HEADER_NAME);
        position = parseIntOrNull(getColumnValue(columnNameToIndex, IlluminaManifest.MAP_INFO_HEADER_NAME));

        ploidy = getColumnValue(columnNameToIndex, IlluminaManifest.PLOIDY_HEADER_NAME);
        species = getColumnValue(columnNameToIndex, IlluminaManifest.SPECIES_HEADER_NAME);
        source = getColumnValue(columnNameToIndex, IlluminaManifest.SOURCE_HEADER_NAME);
        sourceVersion = getColumnValue(columnNameToIndex, IlluminaManifest.SOURCE_VERSION_HEADER_NAME);
        sourceStrand = getSourceStrandFromManifest(columnNameToIndex);
        sourceSeq = getColumnValue(columnNameToIndex, IlluminaManifest.SOURCE_SEQ_HEADER_NAME);
        topGenomicSeq = getColumnValue(columnNameToIndex, IlluminaManifest.TOP_GENOMIC_SEQ_HEADER_NAME);
        beadSetId = parseIntOrNull(getColumnValue(columnNameToIndex, IlluminaManifest.BEAD_SET_ID_HEADER_NAME));
        expClusters = getColumnValueIfPresentInManifest(columnNameToIndex, IlluminaManifest.EXP_CLUSTERS_HEADER_NAME);
        refStrand = getRefStrandFromManifest(columnNameToIndex);
        intensityOnly = getIntensityOnlyFromManifest(columnNameToIndex);

        majorGenomeBuild = getMajorGenomeBuild(genomeBuild);
        hgGenomeBuild = getHgGenomeBuild(majorGenomeBuild);
    }

    private IlluminaStrand getIlmnStrandFromManifest(final Map<String, Integer> columnNameToIndex) {
        return getIlluminaStrandFromManifest(columnNameToIndex, IlluminaManifest.ILLUMINA_STRAND_HEADER_NAME);
    }

    private Strand getRefStrandFromManifest(final Map<String, Integer> columnNameToIndex) {
        final String strandValue = getColumnValueIfPresentInManifest(columnNameToIndex, IlluminaManifest.REF_STRAND_HEADER_NAME);
        if (strandValue == null) {
            return Strand.NONE;
        }
        return Strand.decode(strandValue.charAt(0));
    }

    private IlluminaStrand getSourceStrandFromManifest(final Map<String, Integer> columnNameToIndex) {
        return getIlluminaStrandFromManifest(columnNameToIndex, IlluminaManifest.SOURCE_STRAND_HEADER_NAME);
    }

    private IlluminaStrand getIlluminaStrandFromManifest(final Map<String, Integer> columnNameToIndex, final String strandHeaderName) {
        final String strandValue = getColumnValueIfPresentInManifest(columnNameToIndex, strandHeaderName);
        if (strandValue == null) {
            return IlluminaStrand.NONE;
        }
        return IlluminaStrand.valueOf(strandValue);
    }

    private boolean getIntensityOnlyFromManifest(final Map<String, Integer> columnNameToIndex) {
        final String intensityOnlyValue = getColumnValueIfPresentInManifest(columnNameToIndex, IlluminaManifest.INTENSITY_ONLY_HEADER_NAME);
        if (intensityOnlyValue == null) {
            return false;
        }
        switch (intensityOnlyValue) {
            case "1":
                return true;
            case "0":
                return false;
            default:
                throw new PicardException("Unrecognized value ('" + intensityOnlyValue + "') for '" + IlluminaManifest.INTENSITY_ONLY_HEADER_NAME + "'");
        }
    }

    private String getColumnValueIfPresentInManifest(final Map<String, Integer> columnNameToIndex, String columnHeaderName) {
        Integer index = columnNameToIndex.get(columnHeaderName);
        return index != null ? originalLine[index] : null;
    }

    private String getColumnValue(final Map<String, Integer> columnNameToIndex, String columnHeaderName) {
        return originalLine[columnNameToIndex.get(columnHeaderName)];
    }


    Integer parseIntOrNull(String s) {
        if (s != null) {
            try {
                return Integer.valueOf(s);
            } catch (NumberFormatException e) {};
        }
        return null;
    }

    static boolean isValidManifestRecordLine(final int numColumns, final String[] line) {
        if ((line == null) || (line.length != numColumns)) {
            return false;
        }
        return !line[0].startsWith("[");
    }

    /**
     * Copy Constructor
     */
    IlluminaManifestRecord(final IlluminaManifestRecord record) {
        index = record.getIndex();
        originalLine = record.originalLine;

        ilmnId = record.getIlmnId();
        name = record.getName();
        ilmnStrand = record.getIlmnStrand();
        snp = record.getSnp();

        isIndel = record.isIndel();
        isAmbiguous = record.isAmbiguous();

        addressAId = record.getAddressAId();
        alleleAProbeSeq = record.getAlleleAProbeSeq();
        addressBId = record.getAddressBId();
        alleleBProbeSeq = record.getAlleleBProbeSeq();
        genomeBuild = record.getGenomeBuild();

        majorGenomeBuild = record.getMajorGenomeBuild();
        hgGenomeBuild = record.getHgGenomeBuild();

        chr = record.getChr();
        position = record.getPosition();
        ploidy = record.getPloidy();
        species = record.getSpecies();
        source = record.getSource();
        sourceVersion = record.getSourceVersion();
        sourceStrand = record.getSourceStrand();
        sourceSeq = record.getSourceSeq();
        topGenomicSeq = record.getTopGenomicSeq();
        beadSetId = record.getBeadSetId();
        expClusters = record.getExpClusters();
        refStrand = record.getRefStrand();
        intensityOnly = record.getIntensityOnly();
    }

    public String getLine() {
        return StringUtils.join(originalLine, ",");
    }

    public String getIlmnId() {
        return ilmnId;
    }

    public String getName() {
        return name;
    }

    public IlluminaStrand getIlmnStrand() {
        return ilmnStrand;
    }

    public String getSnp() {
        return snp;
    }

    public boolean isIndel() {
        return isIndel;
    }

    public boolean isAmbiguous() {
        return isAmbiguous;
    }

    public boolean isSnp() {
        return !isIndel;
    }

    public String getAddressAId() {
        return addressAId;
    }

    public String getAlleleAProbeSeq() {
        return alleleAProbeSeq;
    }

    public String getAddressBId() {
        return addressBId;
    }

    public String getAlleleBProbeSeq() {
        return alleleBProbeSeq;
    }

    public String getGenomeBuild() {
        return genomeBuild;
    }

    public String getMajorGenomeBuild() {
        return majorGenomeBuild;
    }

    public String getHgGenomeBuild() {
        return hgGenomeBuild;
    }

    public void setHgGenomeBuild(final String hgGenomeBuild) {
        this.hgGenomeBuild = hgGenomeBuild;
    }

    public String getChr() {
        return chr.trim().equals("XY") ? "X" : chr;
    }

    public int getPosition() {
        return position;
    }

    public String getPloidy() {
        return ploidy;
    }

    public String getSpecies() {
        return species;
    }

    public String getSource() {
        return source;
    }

    public String getSourceVersion() {
        return sourceVersion;
    }

    public IlluminaStrand getSourceStrand() {
        return sourceStrand;
    }

    public String getSourceSeq() {
        return sourceSeq;
    }

    public String getTopGenomicSeq() {
        return topGenomicSeq;
    }

    public int getBeadSetId() {
        return beadSetId;
    }

    public String getExpClusters() {
        return expClusters;
    }

    public Strand getRefStrand() { return refStrand; }

    public boolean getIntensityOnly() { return intensityOnly; }

    private String getMajorGenomeBuild(String genomeBuild) {

        final String majorGenomeBuild;
        // determine majorGenomeBuild and hgGenomeBuild
        final int indexOfDot = genomeBuild.indexOf('.');
        if (indexOfDot != -1) {
            majorGenomeBuild = genomeBuild.substring(0, indexOfDot);
        } else {
            majorGenomeBuild = genomeBuild;
        }

        return majorGenomeBuild;
    }

    private String getHgGenomeBuild(String majorGenomeBuild) {
        return IlluminaManifest.HG_TO_NCBI.inverseBidiMap().containsKey(majorGenomeBuild)
                ? IlluminaManifest.HG_TO_NCBI.inverseBidiMap().get(majorGenomeBuild).toString()
                : "UNKNOWN";
    }

    public int getIndex() {
        return index;
    }
}
