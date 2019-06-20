package picard.pedigree;

/**
 * Represents the sex of an individual.
 */
public enum Sex {
    Male(1), Female(2), Unknown(-9);

    /** The integer code used when reading/writing ped files. */
    private final int code;

    /** Private constructor that takes the pedigree code for sex. */
    Sex(final int code) {
        this.code = code;
    }

    /** Returns the code used to encode this sex in a ped/fam file. */
    public int toCode() { return this.code;}

    /** Decodes the Sex from a numeric code. Note that any value other than 1 or 2 will return Unknown. */
    public static Sex fromCode(final int code) {
        if (code == Male.code) return Male;
        else if (code == Female.code) return Female;
        else return Unknown;
    }
}
