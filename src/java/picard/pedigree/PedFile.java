package picard.pedigree;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;

/**
 * Represents a .ped file of family information as documented here:
 *    http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
 *
 * Stores the information in memory as a map of individualId -> Pedigree information for that individual
 */
public class PedFile extends TreeMap<String, PedFile.PedTrio> {
    private static final Log log = Log.getInstance(PedFile.class);
    static final Pattern WHITESPACE = Pattern.compile("\\s+");
    static final Pattern TAB = Pattern.compile("\\t");
    private final Pattern delimiterPattern;
    private final String delimiterString; // A textual representation of the delimiter, for output purposes

    // These two are really for PedTrio, but they can't be static in there and need to be accessed outside of PedFile
    public static final Number NO_PHENO = new Integer(-9);
    public static final Sex UNKNOWN_SEX = Sex.Unknown;

    public PedFile(final boolean isTabMode) {
        delimiterPattern = isTabMode ? TAB : WHITESPACE;
        delimiterString = isTabMode ? "tabs" : "whitespace";
    }

    /** Adds a trio to the PedFile keyed by the individual id. */
    public void add(final PedTrio trio) {
        put(trio.getIndividualId(), trio);
    }

    /**
     * Writes a set of pedigrees out to disk.
     */
    public void write(final File file) {
        IOUtil.assertFileIsWritable(file);
        final BufferedWriter out = IOUtil.openFileForBufferedWriting(file);

        try {
            for (final PedTrio trio : values()) {
                out.write(trio.getFamilyId());
                out.write("\t");
                out.write(trio.getIndividualId());
                out.write("\t");
                out.write(trio.getPaternalId());
                out.write("\t");
                out.write(trio.getMaternalId());
                out.write("\t");
                out.write(String.valueOf(trio.getSex().toCode()));
                out.write("\t");
                out.write(trio.getPhenotype().toString());
                out.newLine();
            }

            out.close();
        }
        catch (final IOException ioe) {
            throw new RuntimeIOException("IOException while writing to file " + file.getAbsolutePath(), ioe);
        }
    }

    /**
     * Attempts to read a pedigree file into memory.
     */
    public static PedFile fromFile(final File file, final boolean isTabMode) {
        final PedFile pedFile = new PedFile(isTabMode);

        IOUtil.assertFileIsReadable(file);
        for (final String line : IOUtil.readLines(file)) {
            final String[] fields = pedFile.delimiterPattern.split(line);
            if (fields.length != 6) {
                log.error("Ped file line contained invalid number of fields, skipping: " + line);
                continue;
            }

            final PedTrio trio = pedFile.new PedTrio(fields[0],
                    fields[1],
                    fields[2],
                    fields[3],
                    Sex.fromCode(Integer.parseInt(fields[4])),
                    fields[5].contains(".") ? Double.parseDouble(fields[5]) : Integer.parseInt(fields[5])
            );
            pedFile.add(trio);
        }

        return pedFile;
    }

    /**
     * Scans through the pedigrees and removes all entries that do not have both paternal and maternal ids set.
     */
    public PedFile removeIncompleteTrios() {
        final Iterator<Map.Entry<String,PedTrio>> iterator = entrySet().iterator();

        while (iterator.hasNext()) {
            if (!iterator.next().getValue().hasBothParents()) iterator.remove();
        }

        return this;
    }

    public class PedTrio {
        private final String familyId;
        private final String individualId;
        private final String paternalId;
        private final String maternalId;
        private final Sex sex;
        private final Number phenotype;

        /** Constructs a TRIO that cannot be modified after the fact. */
        public PedTrio(final String familyId, final String individualId, final String paternalId, final String maternalId, final Sex sex, final Number phenotype) {
            if (delimiterPattern.split(familyId).length != 1)     throw new IllegalArgumentException("FamilyID     cannot contain " + delimiterString + ": [" + familyId     + "]");
            if (delimiterPattern.split(individualId).length != 1) throw new IllegalArgumentException("IndividualID cannot contain " + delimiterString + ": [" + individualId + "]");
            if (delimiterPattern.split(paternalId).length != 1)   throw new IllegalArgumentException("PaternalID   cannot contain " + delimiterString + ": [" + paternalId   + "]");
            if (delimiterPattern.split(maternalId).length != 1)   throw new IllegalArgumentException("MaternalID   cannot contain " + delimiterString + ": [" + maternalId   + "]");

            this.familyId = familyId;
            this.individualId = individualId;
            this.paternalId = paternalId;
            this.maternalId = maternalId;
            this.sex = sex;
            this.phenotype = phenotype;
        }

        /** True if this record has paternal and maternal ids, otherwise false. */
        public boolean hasBothParents() {
            return this.paternalId != null && this.maternalId != null;
        }

        public String getFamilyId() { return familyId; }
        public String getIndividualId() { return individualId; }
        public String getPaternalId() { return paternalId; }
        public String getMaternalId() { return maternalId; }
        public Sex getSex() { return sex; }
        public Number getPhenotype() { return phenotype; }
    }

    /** Function that accepts a map from sample-name to its sex and creates a PEDFile documenting the sexes.
     * @param sampleSexes a map from sample-name to its sex
     * @return a PedFile object that contains data.
     */
    static public PedFile fromSexMap(final Map<String, Sex> sampleSexes) {

        final PedFile pedfile = new PedFile(true);
        for (final Map.Entry<String, Sex> sampleSex : sampleSexes.entrySet()) {
            final PedFile.PedTrio ped = pedfile.new PedTrio(
                    sampleSex.getKey(), sampleSex.getKey(),
                    "." ,
                    "." ,
                    sampleSex.getValue(), PedFile.NO_PHENO);

            pedfile.add(ped);
        }

        return pedfile;
    }
}
