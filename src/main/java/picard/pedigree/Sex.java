package picard.pedigree;

import picard.PicardException;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents the sex of an individual.
 */
public enum Sex {
    Male("M", 1), Female("F", 2), Unknown("U",-9), NotReported("N",-9);

    /** The integer code used when reading/writing ped files. */
    private final int code;

    /** The single-character symbol used when reading/writing VCF files */
    private final String symbol;

    /** Private constructor that takes the pedigree code for sex. */
    Sex(final String symbol, final int code) {
        this.code = code;
        this.symbol = symbol;
    }

    /** Returns the code used to encode this sex in a ped/fam file. */
    public int toCode() { return this.code;}

    /** Decodes the Sex from a numeric code. Note that any value other than 1 or 2 will return Unknown. */
    public static Sex fromCode(final int code) {
        if (code == Male.code) return Male;
        else if (code == Female.code) return Female;
        else return Unknown;
    }

    /** Returns the single-character symbol used to encode sex in a VCF file */
    public String toSymbol() {
        return this.symbol;
    }

    /** Decodes the Sex from a String. This can be the full string or a single-character symbol representation */
    public static Sex fromString(final String sexString) {
        final Predicate<Sex> match =
                g -> sexString.equalsIgnoreCase(g.symbol)
                        ||   sexString.equalsIgnoreCase(g.name());
        final List<Sex> genders = Stream.of(Sex.values()).filter(match).collect(Collectors.toList());
        if (genders.size() == 1) return genders.get(0);
        throw new PicardException("Unrecognized Sex string: " + sexString);
    }
}
