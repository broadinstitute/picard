package picard.annotation;

import htsjdk.samtools.util.Tuple;
import htsjdk.tribble.annotation.Strand;

import java.util.Arrays;
import java.util.Map;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.util.stream.Collectors.toMap;

/**
 * Class represents one GTF record
 */

public class GtfRecord {
    private static final String COLUMN_DELIMITER = "\t";
    private static final String ATTRIBUTE_DELIMITER = ";";
    private static Pattern pattern = Pattern.compile("([^\\s]++)?+\\s\"([^\"]++)?+\"");
    /**
     * This is the ID of the sequence that is used to establish the coordinate system of the annotation.
     * In the example above, the reference sequence is "Chr1".
     */
    private final String sequenceName;
    /**
     * The source of the annotation. This field describes how the annotation was derived.
     * In the example above, the source is "curated" to indicate that the feature is the result of human curation.
     * The names and versions of software programs are often used for the source field, as in "tRNAScan-SE/1.2".
     */
    private final String source;
    /**
     * The annotation method, also known as type. This field describes the type of the annotation, such as "CDS".
     * Together the method and source describe the annotation type.
     */
    private final Feature feature;
    /**
     * The start of the annotation relative to the reference sequence.
     */
    private final int start;
    /**
     * The stop of the annotation relative to the reference sequence.
     * Start is always less than or equal to stop.
     */
    private final int end;
    /**
     * For annotations that are associated with a numeric score (for example, a sequence similarity),
     * this field describes the score. The score units are completely unspecified, but for sequence similarities,
     * it is typically percent identity. Annotations that do not have a score can use "."
     */
    private final String score;
    /**
     * For those annotations which are strand-specific, this field is the strand on which the annotation resides.
     * It is "+" for the forward strand, "-" for the reverse strand, or "." for annotations that are not stranded.
     */
    private final Strand strand;
    /**
     * For annotations that are linked to proteins, this field describes the phase of the annotation on the codons.
     * It is a number from 0 to 2, or "." for features that have no phase.
     */
    private final String frame;
    /**
     * A semicolon-separated list of tag-value pairs, providing additional information about each feature.
     */
    private final Map<String, String> attributes;

    public GtfRecord(String sequenceName, String source, Feature feature, int start, int end, String score, Strand strand, String frame, Map<String, String> attributes) {
        this.sequenceName = sequenceName;
        this.source = source;
        this.feature = feature;
        this.start = start;
        this.end = end;
        this.score = score;
        this.strand = strand;
        this.frame = frame;
        this.attributes = attributes;
    }

    public static Optional<GtfRecord> fromRow(final String line) {
        final String[] fields = line.split(COLUMN_DELIMITER);
        final Optional<GtfRecord> record;
        if (fields.length == 9) {
            final String sequenceName = fields[0];
            final String source = fields[1];
            final Feature feature = Feature.of(fields[2]);
            final int start = Integer.parseInt(fields[3]);
            final int end = Integer.parseInt(fields[4]);
            final String score = fields[5];
            final Strand strand = Strand.toStrand(fields[6]);
            final String frame = fields[7];
            final Map<String, String> attributes = Arrays.stream(fields[8].split(ATTRIBUTE_DELIMITER))
                    .flatMap(attribute -> {
                        final Matcher matcher = pattern.matcher(attribute);
                        return matcher.find()
                                ? Stream.of(new Tuple<>(matcher.group(1), matcher.group(2)))
                                : Stream.empty();
                    })
                    .collect(toMap(tuple -> tuple.a, tuple -> tuple.b));
            record = Optional.of(
                    new GtfRecord(sequenceName, source, feature, start, end, score, strand, frame, attributes)
            );
        } else {
            record = Optional.empty();
        }
        return record;
    }

    public String toRow() {
        final String strand;
        switch (strand()) {
            case POSITIVE:
                strand = "+";
                break;
            case NEGATIVE:
                strand = "-";
                break;
            case NONE:
            default:
                strand = ".";
                break;
        }
        return String.join(
                COLUMN_DELIMITER,
                sequenceName(), source(), feature().type(),
                Integer.toString(start), Integer.toString(end),
                score(), strand, frame(),
                attributes.entrySet()
                        .stream()
                        .map(entry -> String.format("%s \"%s\"", entry.getKey(), entry.getValue()))
                        .collect(Collectors.joining(ATTRIBUTE_DELIMITER))
        );
    }

    public String sequenceName() {
        return sequenceName;
    }

    public String source() {
        return source;
    }

    public Feature feature() {
        return feature;
    }

    public int start() {
        return start;
    }

    public int end() {
        return end;
    }

    public String score() {
        return score;
    }

    public Strand strand() {
        return strand;
    }

    public String frame() {
        return frame;
    }

    public Map<String, String> attributes() {
        return attributes;
    }
}
