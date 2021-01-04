package picard.arrays.illumina;

/**
 * A simple class to represent a locus entry in an Illumina Bead Pool Manifest (BPM) file
 */
public class IlluminaBPMLocusEntry {
    // IlmnID (probe identifier) of locus
    String ilmnId;

    // Name (variant identifier) of locus
    String name;

    // Index of this entry.
    int index;

    // Illumina Strand value
    String ilmnStrand;

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
    String sourceStrand;

    // Only populated in CSV files or BPM files with version 4 data block
    String sourceSeq;

    // Only populated in CSV files or BPM files with version 4 data block
    String topGenomicSeq;

    int expClusters;
    int intensityOnly;

   // Identifies type of assay (0 - Infinium II , 1 - Infinium I (A/T), 2 - Infinium I (G/C)
    int assayType;

    float fracA;
    float fracC;
    float fracT;
    float fracG;

    // Refstrand annotation
    String refStrand;

    // Not part of the locusEntry record in the BPM, added here for convenience
    int normalizationId;

    public IlluminaBPMLocusEntry() {
        ilmnId = "";
        name = "";
        index = -1;
        ilmnStrand = "";
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
        sourceStrand = "";

        sourceStrand = "";

        expClusters = -1;
        intensityOnly = -1;
        assayType = -1;

        fracA = 0.0f;
        fracC = 0.0f;
        fracT = 0.0f;
        fracG = 0.0f;

        refStrand = "";

        normalizationId = -1;
    }

}
