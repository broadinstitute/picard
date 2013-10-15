package net.sf.picard.pedigree;

/**
 * Represents a single trio within a ped file.
 *
 * @author Tim Fennell
 */
public class PedTrio {
    public static final Number NO_PHENO = new Integer(-9);
    public static final Sex UNKNOWN_SEX = Sex.Unknown;

    private final String familyId;
    private final String individualId;
    private final String paternalId;
    private final String maternalId;
    private final Sex sex;
    private final Number phenotype;

    /** Constructs a TRIO that cannot be modified after the fact. */
    public PedTrio(final String familyId, final String individualId, final String paternalId, final String maternalId, final Sex sex, final Number phenotype) {
        if (PedFile.WHITESPACE.split(familyId).length != 1)     throw new IllegalArgumentException("FamilyID     cannot contain whitespace: [" + familyId     + "]");
        if (PedFile.WHITESPACE.split(individualId).length != 1) throw new IllegalArgumentException("IndividualID cannot contain whitespace: [" + individualId + "]");
        if (PedFile.WHITESPACE.split(paternalId).length != 1)   throw new IllegalArgumentException("PaternalID   cannot contain whitespace: [" + paternalId   + "]");
        if (PedFile.WHITESPACE.split(maternalId).length != 1)   throw new IllegalArgumentException("MaternalID   cannot contain whitespace: [" + maternalId   + "]");

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
