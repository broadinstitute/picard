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
package net.sf.picard.annotation;

import net.sf.picard.util.Log;
import net.sf.picard.util.OverlapDetector;
import net.sf.picard.util.TabbedTextFileWithHeaderParser;
import net.sf.samtools.SAMSequenceDictionary;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Loads gene annotations from a refFlat file into an OverlapDetector<Gene>.  Discards annotations that are not
 * internally consistent, e.g. transcripts on different chromosomes or different strands.
 */
public class RefFlatReader {
    private static final Log LOG = Log.getInstance(RefFlatReader.class);
    // These are in the order that columns appear in refFlat format.
    public enum RefFlatColumns{GENE_NAME, TRANSCRIPT_NAME, CHROMOSOME, STRAND, TX_START, TX_END, CDS_START, CDS_END,
        EXON_COUNT, EXON_STARTS, EXON_ENDS}

    private static final String[] RefFlatColumnLabels = new String[RefFlatColumns.values().length];

    static {
        for (int i = 0; i < RefFlatColumnLabels.length; ++i) {
            RefFlatColumnLabels[i] = RefFlatColumns.values()[i].name();
        }
    }

    private final File refFlatFile;
    private final SAMSequenceDictionary sequenceDictionary;

    RefFlatReader(final File refFlatFile, final SAMSequenceDictionary sequenceDictionary) {
        this.refFlatFile = refFlatFile;
        this.sequenceDictionary = sequenceDictionary;
    }

    static OverlapDetector<Gene> load(final File refFlatFile, final SAMSequenceDictionary sequenceDictionary) {
        return new RefFlatReader(refFlatFile, sequenceDictionary).load();
    }

    OverlapDetector<Gene> load() {
        final OverlapDetector<Gene> overlapDetector = new OverlapDetector<Gene>(0, 0);

        final int expectedColumns = RefFlatColumns.values().length;
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(refFlatFile, RefFlatColumnLabels);
        final Map<String, List<TabbedTextFileWithHeaderParser.Row>> refFlatLinesByGene =
                new HashMap<String, List<TabbedTextFileWithHeaderParser.Row>>();

        for (final TabbedTextFileWithHeaderParser.Row row : parser) {
            final int lineNumber = parser.getCurrentLineNumber() - 1; // getCurrentLineNumber returns the number of the next line
            if (row.getFields().length != expectedColumns) {
                throw new AnnotationException("Wrong number of fields in refFlat file " + refFlatFile + " at line " +
                        lineNumber);
            }
            final String geneName = row.getField(RefFlatColumns.GENE_NAME.name());
            final String transcriptName = row.getField(RefFlatColumns.TRANSCRIPT_NAME.name());
            final String transcriptDescription = geneName + ":" + transcriptName;
            final String chromosome = row.getField(RefFlatColumns.CHROMOSOME.name());
            if (!isSequenceRecognized(chromosome)) {
                LOG.debug("Skipping " + transcriptDescription + " due to unrecognized sequence " + chromosome);
            } else {
                List<TabbedTextFileWithHeaderParser.Row> transcriptLines = refFlatLinesByGene.get(geneName);
                if (transcriptLines == null) {
                    transcriptLines = new ArrayList<TabbedTextFileWithHeaderParser.Row>();
                    refFlatLinesByGene.put(geneName, transcriptLines);
                }
                transcriptLines.add(row);
            }
        }

        int longestInterval = 0;
        int numIntervalsOver1MB = 0;

        for (final List<TabbedTextFileWithHeaderParser.Row> transcriptLines : refFlatLinesByGene.values()) {
            try {
                final Gene gene = makeGeneFromRefFlatLines(transcriptLines);
                overlapDetector.addLhs(gene, gene);
                if (gene.length() > longestInterval) longestInterval = gene.length();
                if (gene.length() > 1000000) ++numIntervalsOver1MB;
            } catch (AnnotationException e) {
                LOG.debug(e.getMessage() + " -- skipping");
            }
        }
        LOG.debug("Longest gene: " + longestInterval + "; number of genes > 1MB: " + numIntervalsOver1MB);
        return overlapDetector;
    }

    private boolean isSequenceRecognized(final String sequence) {
        return (sequenceDictionary.getSequence(sequence) != null);
    }

        private Gene makeGeneFromRefFlatLines(final List<TabbedTextFileWithHeaderParser.Row> transcriptLines) {
            final String geneName = transcriptLines.get(0).getField(RefFlatColumns.GENE_NAME.name());
            final String strandStr = transcriptLines.get(0).getField(RefFlatColumns.STRAND.name());
            final boolean negative = strandStr.equals("-");
            final String chromosome = transcriptLines.get(0).getField(RefFlatColumns.CHROMOSOME.name());
            final List<Gene.Transcript> transcripts = new ArrayList<Gene.Transcript>(transcriptLines.size());
            for (final TabbedTextFileWithHeaderParser.Row row: transcriptLines) {
                if (!strandStr.equals(row.getField(RefFlatColumns.STRAND.name()))) {
                    throw new AnnotationException("Strand disagreement in refFlat file for gene " + geneName);
                }
                if (!chromosome.equals(row.getField(RefFlatColumns.CHROMOSOME.name()))) {
                    throw new AnnotationException("Chromosome disagreement(" + chromosome + " != " + row.getField(RefFlatColumns.CHROMOSOME.name()) +
                            ") in refFlat file for gene " + geneName);
                }
                transcripts.add(makeTranscriptFromRefFlatLine(row));
            }
            int start = Integer.MAX_VALUE;
            int end = Integer.MIN_VALUE;
            for (final Gene.Transcript transcript : transcripts) {
                start = Math.min(start, transcript.start());
                end = Math.max(end, transcript.end());
            }
            return new Gene(chromosome, start, end, negative, geneName, transcripts);
        }

    /**
     * Conversion from 0-based half-open to 1-based inclusive intervals is done here.
     */
    private Gene.Transcript makeTranscriptFromRefFlatLine(final TabbedTextFileWithHeaderParser.Row row) {
            final String geneName = row.getField(RefFlatColumns.GENE_NAME.name());
            final String transcriptName = row.getField(RefFlatColumns.TRANSCRIPT_NAME.name());
            final String transcriptDescription = geneName + ":" + transcriptName;
            final int exonCount = Integer.parseInt(row.getField(RefFlatColumns.EXON_COUNT.name()));
            final String[] exonStarts = row.getField(RefFlatColumns.EXON_STARTS.name()).split(",");
            final String[] exonEnds = row.getField(RefFlatColumns.EXON_ENDS.name()).split(",");
            if (exonCount != exonStarts.length) {
                throw new AnnotationException("Number of exon starts does not agree with number of exons for " +
                        transcriptDescription);
            }
            if (exonCount != exonEnds.length) {
                throw new AnnotationException("Number of exon ends does not agree with number of exons for " +
                        transcriptDescription);
            }
            final Gene.Transcript.Exon[] exons = new Gene.Transcript.Exon[exonCount];
            for (int i = 0; i < exonCount; ++i) {
                exons[i] = new Gene.Transcript.Exon(Integer.parseInt(exonStarts[i]) + 1, Integer.parseInt(exonEnds[i]));
                if (exons[i].start >= exons[i].end) {
                    throw new AnnotationException("Exon has 0 or negative extent for " + transcriptDescription);
                }
                if (i > 0 && exons[i-1].end >= exons[i].start) {
                    throw new AnnotationException("Exons overlap for " + transcriptDescription);
                }
            }
            final int transcriptionStart = row.getIntegerField(RefFlatColumns.TX_START.name()) + 1;
            final int transcriptionEnd = row.getIntegerField(RefFlatColumns.TX_END.name());
            final int codingStart = row.getIntegerField(RefFlatColumns.CDS_START.name()) + 1;
            final int codingEnd = row.getIntegerField(RefFlatColumns.CDS_END.name());
            return new Gene.Transcript(transcriptName, transcriptionStart, transcriptionEnd, codingStart, codingEnd, exons);
        }

}
