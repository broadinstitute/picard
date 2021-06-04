package picard.arrays.illumina;

import htsjdk.tribble.annotation.Strand;

/**
 * A simple class to represent a locus entry in an Illumina Bead Pool Manifest (BPM) file
 */
public class IlluminaBPMLocusEntry {
    int version;        // The LocusEntry version.

    // IlmnID (probe identifier) of locus
    String ilmnId;

    // Name (variant identifier) of locus
    String name;

    // Index of this entry.
    int index;

    // Illumina Strand value
    IlluminaManifestRecord.IlluminaStrand ilmnStrand;

    // SNP value for locus (e.g., [A/C])
    String snp;

    // Chromosome for the locus (e.g., XY)
    String chrom;

    String ploidy;

    String species;

    // Mapping location of locus
    int mapInfo;

    // Customer Strand
    String customerStrand;

    // AddressA ID of locus
    int addressA;

    // Only populated in CSV files or BPM files with version 4 data block
    String alleleAProbeSeq;

    // AddressB ID of locus (0 if none)
    int addressB;

    // Only populated in CSV files or BPM files with version 4 data block (empty if none)
    String alleleBProbeSeq;

    String genomeBuild;
    String source;
    String sourceVersion;
    IlluminaManifestRecord.IlluminaStrand sourceStrand;

    // Only populated in CSV files or BPM files with version 4 data block
    String sourceSeq;

    // Only populated in CSV files or BPM files with version 4 data block
    String topGenomicSeq;

    int expClusters;
    boolean intensityOnly;

   // Identifies type of assay (0 - Infinium II , 1 - Infinium I (A/T), 2 - Infinium I (G/C)
    int assayType;

    float fracA;
    float fracC;
    float fracT;
    float fracG;

    // Refstrand annotation
    Strand refStrand;

    // Not part of the locusEntry record in the BPM, added here for convenience
    int normalizationId;

    public IlluminaBPMLocusEntry() {
        version = -1;

        ilmnId = "";
        name = "";
        index = -1;
        ilmnStrand = IlluminaManifestRecord.IlluminaStrand.NONE;
        snp = "";
        chrom = "";
        ploidy = "";
        species = "";
        mapInfo = -1;
        customerStrand = "";

        addressA = -1;
        addressB = -1;

        genomeBuild = "";
        source = "";
        sourceVersion = "";
        sourceStrand = IlluminaManifestRecord.IlluminaStrand.NONE;

        expClusters = -1;
        intensityOnly = false;
        assayType = -1;

        fracA = 0.0f;
        fracC = 0.0f;
        fracT = 0.0f;
        fracG = 0.0f;

        refStrand = Strand.NONE;

        normalizationId = -1;
    }

    public int getVersion() {
        return version;
    }

    public String getIlmnId() {
        return ilmnId;
    }

    public String getName() {
        return name;
    }

    public int getIndex() {
        return index;
    }

    public IlluminaManifestRecord.IlluminaStrand getIlmnStrand() {
        return ilmnStrand;
    }

    public String getSnp() {
        return snp;
    }

    public String getChrom() {
        return chrom;
    }

    public String getPloidy() {
        return ploidy;
    }

    public String getSpecies() {
        return species;
    }

    public int getMapInfo() {
        return mapInfo;
    }

    public String getCustomerStrand() {
        return customerStrand;
    }

    public int getAddressA() {
        return addressA;
    }

    public String getAlleleAProbeSeq() {
        return alleleAProbeSeq;
    }

    public int getAddressB() {
        return addressB;
    }

    public String getAlleleBProbeSeq() {
        return alleleBProbeSeq;
    }

    public String getGenomeBuild() {
        return genomeBuild;
    }

    public String getSource() {
        return source;
    }

    public String getSourceVersion() {
        return sourceVersion;
    }

    public IlluminaManifestRecord.IlluminaStrand getSourceStrand() {
        return sourceStrand;
    }

    public String getSourceSeq() {
        return sourceSeq;
    }

    public String getTopGenomicSeq() {
        return topGenomicSeq;
    }

    public int getExpClusters() {
        return expClusters;
    }

    public boolean isIntensityOnly() {
        return intensityOnly;
    }

    public int getAssayType() {
        return assayType;
    }

    public float getFracA() {
        return fracA;
    }

    public float getFracC() {
        return fracC;
    }

    public float getFracT() {
        return fracT;
    }

    public float getFracG() {
        return fracG;
    }

    public Strand getRefStrand() {
        return refStrand;
    }

    public int getNormalizationId() {
        return normalizationId;
    }
}
