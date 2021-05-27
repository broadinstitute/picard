package picard.arrays.illumina;

import htsjdk.tribble.annotation.Strand;
import picard.PicardException;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Set;
import java.util.TreeSet;

/**
 * A class to parse the contents of an Illumina Bead Pool Manifest (BPM) file
 *
 * A BPM file contains metadata (including the alleles, mapping and normalization information) on an Illumina Genotyping Array
 * Each type of genotyping array has a specific BPM .
 *
 */
public class IlluminaBPMFile extends InfiniumDataFile implements AutoCloseable {
    private static final String BPM_IDENTIFIER = "BPM";

    private String manifestName;
    private String controlConfig;

    private int numLoci = 0;
    private IlluminaBPMLocusEntry[] locusEntries = null;

    private int[] allNormalizationIds = null;
    private Integer[] uniqueNormalizationIds = null;

    public IlluminaBPMFile(final File bpmFile) throws IOException {
        super(new DataInputStream(new FileInputStream(bpmFile)), true);
        parse();
    }

    @Override
    public void close() throws IOException {
        stream.close();
    }

    /**
     * Main parsing method.
     *
     * @throws IOException thrown when there is a problem reading the Bead Pool Manifest
     */
    private void parse() throws IOException {

        try {
            final byte[] formatIdentifier = new byte[BPM_IDENTIFIER.length()];
            for (int i = 0; i < formatIdentifier.length; i++) {
                formatIdentifier[i] = parseByte();
            }

            final String identifier = new String(formatIdentifier);
            setIdentifier(identifier);
            if (!identifier.equals(BPM_IDENTIFIER)) {
                throw new PicardException("Invalid identifier '" + identifier + "' for BPM file");
            }
            setFileVersion(parseByte());
            if (getFileVersion() != 1) {
                throw new PicardException("Unknown BPM version (" + getFileVersion() + ")");
            }
            int version = parseInt();
            final int versionFlag = 0x1000;
            if ((version & versionFlag) == versionFlag) {
                version ^= versionFlag;
            }
            if (version > 5 || version < 3) {
                throw new PicardException("Unsupported BPM version (" + version + ")");
            }
            manifestName = parseString();
            controlConfig = parseString();
            numLoci = parseInt();

            readData();
        } finally {
            stream.close();
        }
    }

    private void readData() throws IOException {
        // Skip the index block
        stream.skipBytes(4 * numLoci);

        // Read the names
        String[] names = new String[numLoci];
        for (int i = 0; i < numLoci; i++) {
            names[i] = parseString();
        }

        // Read the normalization ids.
        allNormalizationIds = parseByteArrayAsInts(numLoci);

        // Create an ordered unique set of normalization ids.
        final Set<Integer> uniqueNormalizationIdsSet = new TreeSet<>();

        // Initialize the locus entries
        locusEntries = new IlluminaBPMLocusEntry[numLoci];

        // Read the locus entries.
        for (int i = 0; i < numLoci; i++) {
            IlluminaBPMLocusEntry locusEntry = parseLocusEntry();
            int normId = allNormalizationIds[locusEntry.index];
            if (normId > 100) {
                throw new PicardException("Invalid normalization ID: " + normId + " for name: " + locusEntry.name);
            }
            locusEntry.normalizationId = normId + 100 * locusEntry.assayType;
            allNormalizationIds[locusEntry.index] = locusEntry.normalizationId;
            uniqueNormalizationIdsSet.add(locusEntry.normalizationId);

            if (locusEntries[locusEntry.index] != null) {
                throw new PicardException("Duplicate locus entry for index: " + locusEntry.index + " '" + locusEntry.name);
            }
            locusEntries[locusEntry.index] = locusEntry;
            if (!names[locusEntry.index].equals(locusEntry.name)) {
                throw new PicardException("Mismatch in names at index: " + locusEntry.index);
            }
        }

        uniqueNormalizationIds = uniqueNormalizationIdsSet.toArray(new Integer[0]);
    }

    private IlluminaBPMLocusEntry parseLocusEntry() throws IOException {
        IlluminaBPMLocusEntry locusEntry = new IlluminaBPMLocusEntry();
        locusEntry.version = parseInt();
        if (locusEntry.version < 6  || locusEntry.version > 8) {
            throw new PicardException("Unsupported Locus version: " + locusEntry.version);
        }

        locusEntry.ilmnId = parseString();
        locusEntry.name = parseString();
        parseString();
        parseString();
        parseString();
        locusEntry.index = parseInt() - 1;
        parseString();
        final String ilmnStrandString = parseString();
        locusEntry.ilmnStrand = (!ilmnStrandString.equals("")) ? IlluminaManifestRecord.IlluminaStrand.valueOf(ilmnStrandString) : IlluminaManifestRecord.IlluminaStrand.NONE;
        locusEntry.snp = parseString();
        locusEntry.chrom = parseString();
        locusEntry.ploidy = parseString();
        locusEntry.species = parseString();
        locusEntry.mapInfo = Integer.parseInt(parseString());
        locusEntry.topGenomicSeq = parseString();       // Only in Version 4
        locusEntry.customerStrand = parseString();
        locusEntry.addressA = parseInt();
        locusEntry.addressB = parseInt();
        locusEntry.alleleAProbeSeq = parseString();     // Only in Version 4
        locusEntry.alleleBProbeSeq = parseString();     // Only in Version 4
        locusEntry.genomeBuild = parseString();
        locusEntry.source = parseString();
        locusEntry.sourceVersion = parseString();
        final String sourceStrandString = parseString();
        locusEntry.sourceStrand = (!sourceStrandString.equals("")) ? IlluminaManifestRecord.IlluminaStrand.valueOf(sourceStrandString) : IlluminaManifestRecord.IlluminaStrand.NONE;
        locusEntry.sourceSeq = parseString();           // Only in Version 4
        parseByte();
        locusEntry.expClusters = parseByte();
        int intensityOnlyValue = parseByte();
        if ((intensityOnlyValue == 0) || (intensityOnlyValue == 1)) {
            locusEntry.intensityOnly = intensityOnlyValue == 1;
        } else {
            throw new PicardException("Unexpected value ('" + intensityOnlyValue + "') for intensity_only field");
        }
        locusEntry.assayType = parseByte();

        if (locusEntry.version >= 7) {
            locusEntry.fracA = parseFloat();
            locusEntry.fracC = parseFloat();
            locusEntry.fracT = parseFloat();
            locusEntry.fracG = parseFloat();
        }
        if (locusEntry.version == 8) {
            final String refStrandString = parseString();
            locusEntry.refStrand = (!refStrandString.equals("")) ? Strand.decode(refStrandString.charAt(0)) : Strand.NONE;
        }

        if (locusEntry.assayType < 0 || locusEntry.assayType > 2) {
            throw new PicardException("Invalid assay_type '" + locusEntry.assayType + "' in BPM file");
        }
        if ((locusEntry.assayType != 0 && locusEntry.addressB == 0) || (locusEntry.assayType == 0 && locusEntry.addressB != 0)) {
            throw new PicardException("Invalid assay_type '" + locusEntry.assayType + "' for address B '" + locusEntry.addressB + "' in BPM file");
        }
        return locusEntry;
    }

    public String getManifestName() {
        return manifestName;
    }

    public String getControlConfig() {
        return controlConfig;
    }

    public IlluminaBPMLocusEntry[] getLocusEntries() {
        return locusEntries;
    }

    public int[] getAllNormalizationIds() {
        return allNormalizationIds;
    }

    public Integer[] getUniqueNormalizationIds() {
        return uniqueNormalizationIds;
    }
}
