/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.annotation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.annotation.Strand;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Load gene annotations into an GeneOverlapDetectorLoader of Gene objects.
 */
public class GeneOverlapDetectorLoader {
    private static final Log log = Log.getInstance(GeneOverlapDetectorLoader.class);

    private final SAMSequenceDictionary dictionary;

    private GeneOverlapDetectorLoader(final SAMSequenceDictionary dictionary) {
        this.dictionary = dictionary;
    }

    /**
     * Calls GeneOverlapDetectorLoader.load()
     * @param records given records
     * @param dictionary given sequenceDictionary
     * @return OverlapDetector of Gene objects
     */
    public static OverlapDetector<Gene> load(final Stream<GenePredRecord> records, final SAMSequenceDictionary dictionary) {
        return new GeneOverlapDetectorLoader(dictionary).load(records);
    }

    private OverlapDetector<Gene> load(final Stream<GenePredRecord> records) {
        final OverlapDetector<Gene> overlapDetector = new OverlapDetector<>(0, 0);

        final Map<String, List<GenePredRecord>> geneGroupedRecords;
        geneGroupedRecords = records.filter(this::isSequenceRecognized)
                                    .collect(Collectors.groupingBy(GenePredRecord::name));
        int longestInterval = 0;
        int numIntervalsOver1MB = 0;
        for (final Map.Entry<String, List<GenePredRecord>> geneRows : geneGroupedRecords.entrySet()) {
            try {
                final String geneName = geneRows.getKey();
                final List<GenePredRecord> rows = geneRows.getValue();
                final Gene gene = makeGene(geneName, rows);
                overlapDetector.addLhs(gene, gene);
                if (gene.length() > longestInterval) longestInterval = gene.length();
                if (gene.length() > 1000000) numIntervalsOver1MB++;
            } catch (final AnnotationException exception) {
                log.debug(String.format("%s -- skipping", exception.getMessage()));
            }
        }
        log.debug(String.format("Longest gene: %d; number of genes > 1MB: %d", longestInterval, numIntervalsOver1MB));

        return overlapDetector;
    }

    private Gene makeGene(final String geneName, final List<GenePredRecord> records) {
        final GenePredRecord firstRecord = records.get(0);
        final boolean negative = firstRecord.strand().equals(Strand.NEGATIVE);
        final String chromosome = firstRecord.chromosomeName();

        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        for (final GenePredRecord record : records) {
            start = Math.min(start, record.transcriptStart());
            end = Math.max(start, record.transcriptEnd());
        }

        final Gene gene = new Gene(chromosome, start + 1, end, negative, geneName);

        for (final GenePredRecord record : records) {
            final Gene.Transcript transcript = gene.addTranscript(
                    record.name(),
                    record.transcriptStart() + 1,
                    record.transcriptEnd(),
                    record.codingStart() + 1,
                    record.codingEnd(),
                    record.exonCount()
            );

            final int[] starts = record.exonStarts().toArray();
            final int[] ends = record.exonEnds().toArray();
            for (int i = 0; i < record.exonCount(); i++) {
                transcript.addExon(starts[i] + 1, ends[i]);
            }
        }

        return gene;
    }

    private boolean isSequenceRecognized(final GenePredRecord record) {
        final String sequenceName = record.chromosomeName();
        final boolean recognized = Objects.nonNull(dictionary().getSequence(sequenceName));
        if (!recognized) log.debug(String.format("Skipping due to unrecognized sequence %s", sequenceName));
        return recognized;
    }

    private SAMSequenceDictionary dictionary() {
        return dictionary;
    }
}