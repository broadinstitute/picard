package picard.annotation;

import htsjdk.tribble.annotation.Strand;

import java.util.Arrays;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Class represents one GenePred record
 */
public class GenePredRecord {
    private static final String COLUMN_DELIMITER = "\t";
    private static final String COORDINATE_DELIMITER = ",";

    /**
     * Name of gene
     */
    private final String geneName;
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

    public GenePredRecord(String geneName, String chromosomeName, Strand strand, int transcriptStart, int transcriptEnd, int codingStart, int codingEnd, int exonCount, int[] exonStarts, int[] exonEnds) {
        this.geneName = geneName;
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

    public static Optional<GenePredRecord> fromRow(final String line) {
        final String[] fields = line.split(COLUMN_DELIMITER);
        final String name = fields[0];
        final String chromosomeName = fields[1];
        final Strand strand = Strand.toStrand(fields[2]);
        final int txStart = Integer.parseInt(fields[3]);
        final int txEnd = Integer.parseInt(fields[4]);
        final int cdsStart = Integer.parseInt(fields[5]);
        final int cdsEnd = Integer.parseInt(fields[6]);
        final int exonCount = Integer.parseInt(fields[7]);
        final int[] exonStarts = Arrays.stream(fields[8].split(COORDINATE_DELIMITER))
                .mapToInt(Integer::parseInt)
                .toArray();
        final int[] exonEnds = Arrays.stream(fields[9].split(COORDINATE_DELIMITER))
                .mapToInt(Integer::parseInt)
                .toArray();
        return Optional.of(
                new GenePredRecord(name, chromosomeName, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds)
        );
    }

    public String name() {
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

    @Override
    public String toString() {
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
                name(), chromosomeName(), strand,
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
                )
        );
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GenePredRecord that = (GenePredRecord) o;

        if (transcriptStart != that.transcriptStart) return false;
        if (transcriptEnd != that.transcriptEnd) return false;
        if (codingStart != that.codingStart) return false;
        if (codingEnd != that.codingEnd) return false;
        if (exonCount != that.exonCount) return false;
        if (geneName != null ? !geneName.equals(that.geneName) : that.geneName != null) return false;
        if (chromosomeName != null ? !chromosomeName.equals(that.chromosomeName) : that.chromosomeName != null)
            return false;
        if (strand != that.strand) return false;
        if (!Arrays.equals(exonStarts, that.exonStarts)) return false;
        return Arrays.equals(exonEnds, that.exonEnds);
    }

    @Override
    public int hashCode() {
        int result = geneName != null ? geneName.hashCode() : 0;
        result = 31 * result + (chromosomeName != null ? chromosomeName.hashCode() : 0);
        result = 31 * result + (strand != null ? strand.hashCode() : 0);
        result = 31 * result + transcriptStart;
        result = 31 * result + transcriptEnd;
        result = 31 * result + codingStart;
        result = 31 * result + codingEnd;
        result = 31 * result + exonCount;
        result = 31 * result + Arrays.hashCode(exonStarts);
        result = 31 * result + Arrays.hashCode(exonEnds);
        return result;
    }
}
