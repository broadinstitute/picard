package picard.annotation;

import java.util.Arrays;

/**
 * Enum describing possible types of Feature field of GTF file record.
 */
public enum Feature {
    GENE("gene"),
    TRANSCRIPT("transcript"),
    EXON("exon"),
    CDS("CDS"),
    UTR("UTR"),
    SELENOCYSTEINE("selenocysteine"),
    START_CODON("start_codon"),
    STOP_CODON("stop_codon"),
    VARIATION("variation"),
    SIMILARITY("similarity"),
    UNKNOWN("?");

    private final String typeName;

    Feature(final String type) {
        this.typeName = type;
    }

    public String type() {
        return typeName;
    }

    public static Feature of(final String type) {
        return Arrays.stream(values())
                     .filter(feature -> feature.typeName.equals(type))
                     .findFirst()
                     .orElse(UNKNOWN);
    }
}
