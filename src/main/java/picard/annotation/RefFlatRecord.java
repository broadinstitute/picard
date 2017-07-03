package picard.annotation;

import htsjdk.tribble.annotation.Strand;

import java.util.Arrays;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Class represents one RefFlat record
 */
public class RefFlatRecord {
    private static final String COLUMN_DELIMITER = "\t";
    private static final String COORDINATE_DELIMITER = ",";
    /**
     * Name of gene as it appears in Genome Browser.
     */
    private final String geneName;
    /**
     * Name of gene
     */
    private final String name;
    /**
     * Chromosome name
     */
    private final String chromosomeName;
    /**
     * {@link Strand}
     */
    private final Strand strand;
    /**
     * Transcription start position
     */
    private final int transcriptStart;
    /**
     * Transcription end position
     */
    private final int transcriptEnd;
    /**
     * Coding region start
     */
    private final int codingStart;
    /**
     * Coding region end
     */
    private final int codingEnd;
    /**
     * Number of exons
     */
    private final int exonCount;
    /**
     * Exon start positions
     */
    private final int[] exonStarts;
    /**
     * Exon end positions
     */
    private final int[] exonEnds;

    public RefFlatRecord(String geneName, String name, String chromosomeName, Strand strand, int transcriptStart, int transcriptEnd, int codingStart, int codingEnd, int exonCount, int[] exonStarts, int[] exonEnds) {
        this.geneName = geneName;
        this.name = name;
        this.chromosomeName = chromosomeName;
        this.strand = strand;
        this.transcriptStart = transcriptStart;
        this.transcriptEnd = transcriptEnd;
        this.codingStart = codingStart;
        this.codingEnd = codingEnd;
        this.exonCount = exonCount;
        this.exonStarts = exonStarts;
        this.exonEnds = exonEnds;
    }

    public static Optional<RefFlatRecord> fromRow(final String line) {
        final String[] fields = line.split(COLUMN_DELIMITER);
        final Optional<RefFlatRecord> record;
        if (fields.length == 11) {
            final String geneName = fields[0];
            final String name = fields[1];
            final String chromosomeName = fields[2];
            final Strand strand = Strand.toStrand(fields[3]);
            final int txStart = Integer.parseInt(fields[4]);
            final int txEnd = Integer.parseInt(fields[5]);
            final int cdsStart = Integer.parseInt(fields[6]);
            final int cdsEnd = Integer.parseInt(fields[7]);
            final int exonCount = Integer.parseInt(fields[8]);
            final int[] exonStarts = Arrays.stream(fields[9].split(COORDINATE_DELIMITER))
                    .mapToInt(Integer::parseInt)
                    .toArray();
            final int[] exonEnds = Arrays.stream(fields[10].split(COORDINATE_DELIMITER))
                    .mapToInt(Integer::parseInt)
                    .toArray();
            record = Optional.of(
                    new RefFlatRecord(geneName, name, chromosomeName, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds)
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
                geneName(), name(), chromosomeName(), strand,
                Integer.toString(transcriptStart()), Integer.toString(transcriptEnd()),
                Integer.toString(codingStart()), Integer.toString(codingEnd()),
                Integer.toString(exonCount()),
                String.join(
                        COORDINATE_DELIMITER,
                        exonStarts().mapToObj(Integer::toString).collect(Collectors.toList())
                ),
                String.join(
                        COORDINATE_DELIMITER,
                        exonEnds().mapToObj(Integer::toString).collect(Collectors.toList())
                ));
    }

    public String name() {
        return name;
    }

    public String geneName() {
        return geneName;
    }

    public String chromosomeName() {
        return chromosomeName;
    }

    public int exonCount() {
        return exonCount;
    }

    public IntStream exonStarts() {
        return Arrays.stream(exonStarts);
    }

    public IntStream exonEnds() {
        return Arrays.stream(exonEnds);
    }

    public int transcriptStart() {
        return transcriptStart;
    }

    public int transcriptEnd() {
        return transcriptEnd;
    }

    public int codingStart() {
        return codingStart;
    }

    public int codingEnd() {
        return codingEnd;
    }

    public Strand strand() {
        return strand;
    }
}
