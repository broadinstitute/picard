package picard.annotation;

import htsjdk.tribble.annotation.Strand;

import java.util.*;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Class for convertation of RefFlat and GTF formats ecords to GenePred format records.
 */
public class Converter {
    public static final Stream<GenePredRecord> refFlatToGenePred(final Stream<RefFlatRecord> records) {
        return records.map(
                record -> new GenePredRecord(
                        record.name(),
                        record.chromosomeName(),
                        record.strand(),
                        record.transcriptStart(),
                        record.transcriptEnd(),
                        record.codingStart(),
                        record.codingEnd(),
                        record.exonCount(),
                        record.exonStarts().toArray(),
                        record.exonEnds().toArray()
                )
        );
    }

    public static Stream<GenePredRecord> gtfToGenePred(final Stream<GtfRecord> records) {
        final Spliterator<GtfRecord> gtfSpliterator = records.spliterator();
        final Iterator<GtfRecord> gtfIterator = Spliterators.iterator(gtfSpliterator);
        final int characteristics = gtfSpliterator.characteristics() & ~(Spliterator.SIZED | Spliterator.SORTED | Spliterator.SUBSIZED);
        final Iterator<GenePredRecord> iterator = new Iterator<GenePredRecord>() {
            private static final int UNDEFINED = -1;
            private GtfRecord firstInTranscript = gtfIterator.next();
            private boolean transcriptProcessed = false;

            @Override
            public boolean hasNext() {
                return gtfIterator.hasNext();
            }

            private GtfRecord take() {
                final GtfRecord record;
                if (transcriptProcessed) {
                    record = gtfIterator.next();
                } else {
                    record = firstInTranscript;
                    transcriptProcessed = true;
                }
                return record;
            }

            @Override
            public GenePredRecord next() {
                final String transcriptId = firstInTranscript.attributes().getOrDefault("transcript_id", "undefined");
                final String chromosomeName = firstInTranscript.sequenceName();
                final Strand strand = firstInTranscript.strand();
                final List<Integer> exonStarts = new ArrayList<>();
                final List<Integer> exonEnds = new ArrayList<>();
                int codingStart = UNDEFINED;
                int codingEnd = UNDEFINED;

                while (gtfIterator.hasNext()) {
                    final GtfRecord record = take();
                    if (!transcriptId.equals(record.attributes().getOrDefault("transcript_id", "undefined"))) {
                        firstInTranscript = record;
                        transcriptProcessed = false;
                        break;
                    }

                    final Feature feature = record.feature();
                    if (feature.equals(Feature.EXON)) {
                        exonStarts.add(record.start() - 1);
                        exonEnds.add(record.end());
                    }
                    if (strand.equals(Strand.POSITIVE)) {
                        if (feature.equals(Feature.START_CODON) && codingStart == UNDEFINED) {
                            codingStart = record.start();
                        }
                        if (feature.equals(Feature.STOP_CODON)) {
                            codingEnd = record.end();
                        }
                        if (feature.equals(Feature.CDS) && codingStart == UNDEFINED) {
                            codingStart = record.start();
                        }

                        if (feature.equals(Feature.CDS) && codingEnd == UNDEFINED) {
                            codingEnd = record.end();
                        }

                    } else {
                        if (feature.equals(Feature.STOP_CODON) && codingStart == UNDEFINED) {
                            codingStart = record.start();
                        }
                        if (feature.equals(Feature.START_CODON)) {
                            codingEnd = record.end();
                        }
                        if (feature.equals(Feature.CDS) && codingStart == UNDEFINED) {
                            codingStart = record.end();
                        }
                        if (feature.equals(Feature.CDS) && codingEnd == UNDEFINED) {
                            codingEnd = record.start();
                        }
                    }

                }

                codingStart = codingStart == UNDEFINED ? exonEnds.get(exonEnds.size() - 1) : codingStart - 1;
                codingEnd = codingEnd == UNDEFINED ? exonEnds.get(exonEnds.size() - 1) : codingEnd;

                final int transcriptStart = exonStarts.get(0);
                final int transcriptEnd = exonEnds.get(exonEnds.size()-1);
                final int exonCount = exonStarts.size();

                return new GenePredRecord(
                        transcriptId,
                        chromosomeName,
                        strand,
                        transcriptStart,
                        transcriptEnd,
                        codingStart,
                        codingEnd,
                        exonCount,
                        exonStarts.stream().mapToInt(Integer::intValue).toArray(),
                        exonEnds.stream().mapToInt(Integer::intValue).toArray()
                );
            }
        };

        return StreamSupport.stream(Spliterators.spliterator(iterator, -1, characteristics), false);
    }
}
