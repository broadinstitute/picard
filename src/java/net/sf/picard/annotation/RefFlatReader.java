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
import net.sf.picard.util.TabbedInputParser;
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

    private final File refFlatFile;
    private final SAMSequenceDictionary sequenceDictionary;

    RefFlatReader(File refFlatFile, SAMSequenceDictionary sequenceDictionary) {
        this.refFlatFile = refFlatFile;
        this.sequenceDictionary = sequenceDictionary;
    }

    OverlapDetector<Gene> load() {
            final OverlapDetector<Gene> overlapDetector = new OverlapDetector<Gene>(0, 0);

            final int expectedColumns = RefFlatColumns.values().length;
            TabbedInputParser parser = new TabbedInputParser(false, refFlatFile);
            Map<String, List<String[]>> refFlatLinesByGene = new HashMap<String, List<String[]>>();
            for (final String[] fields : parser) {
                final int lineNumber = parser.getCurrentLineNumber() - 1; // getCurrentLineNumber returns the number of the next line
                if (fields.length != expectedColumns) {
                    throw new AnnotationException("Wrong number of fields in refFlat file " + refFlatFile + " at line " +
                            lineNumber);
                }
                final String geneName = fields[RefFlatColumns.GENE_NAME.ordinal()];
                final String transcriptName = fields[RefFlatColumns.TRANSCRIPT_NAME.ordinal()];
                final String transcriptDescription = geneName + ":" + transcriptName;
                String chromosome = fields[RefFlatColumns.CHROMOSOME.ordinal()];
                if (!isSequenceRecognized(chromosome)) {
                    LOG.info("Skipping " + transcriptDescription + " due to unrecognized sequence " + chromosome);
                } else {
                    List<String[]> transcriptLines = refFlatLinesByGene.get(geneName);
                    if (transcriptLines == null) {
                        transcriptLines = new ArrayList<String[]>();
                        refFlatLinesByGene.put(geneName, transcriptLines);
                    }
                    transcriptLines.add(fields);
                }
            }
            for (final List<String[]> transcriptLines : refFlatLinesByGene.values()) {
                try {
                    final Gene gene = makeGeneFromRefFlatLines(transcriptLines);
                    overlapDetector.addLhs(gene, gene);
                } catch (AnnotationException e) {
                    LOG.info(e.getMessage());
                }
            }
            return overlapDetector;
        }

    private boolean isSequenceRecognized(String sequence) {
        return (sequenceDictionary.getSequence(sequence) != null);
    }

        private Gene makeGeneFromRefFlatLines(List<String[]> transcriptLines) {
            final String geneName = transcriptLines.get(0)[RefFlatColumns.GENE_NAME.ordinal()];
            final String strandStr = transcriptLines.get(0)[RefFlatColumns.STRAND.ordinal()];
            final boolean negative = strandStr.equals("-");
            final String chromosome = transcriptLines.get(0)[RefFlatColumns.CHROMOSOME.ordinal()];
            final List<Transcript> transcripts = new ArrayList<Transcript>(transcriptLines.size());
            for (final String[] fields: transcriptLines) {
                if (!strandStr.equals(fields[RefFlatColumns.STRAND.ordinal()])) {
                    throw new AnnotationException("Strand disagreement in refFlat file for gene " + geneName + " -- ignoring gene");
                }
                if (!chromosome.equals(fields[RefFlatColumns.CHROMOSOME.ordinal()])) {
                    throw new AnnotationException("Chromosome disagreement(" + chromosome + " != " + fields[RefFlatColumns.CHROMOSOME.ordinal()] +
                            ") in refFlat file for gene " + geneName + " -- ignoring gene");
                }
                transcripts.add(makeTranscriptFromRefFlatLine(fields));
            }
            int start = Integer.MAX_VALUE;
            int end = Integer.MIN_VALUE;
            for (final Transcript transcript : transcripts) {
                start = Math.min(start, transcript.start());
                end = Math.max(end, transcript.end());
            }
            return new Gene(chromosome, start, end, negative, geneName, transcripts);
        }

        private Transcript makeTranscriptFromRefFlatLine(final String[] fields) {
            final String geneName = fields[RefFlatColumns.GENE_NAME.ordinal()];
            final String transcriptName = fields[RefFlatColumns.TRANSCRIPT_NAME.ordinal()];
            final String transcriptDescription = geneName + ":" + transcriptName;
            final int exonCount = Integer.parseInt(fields[RefFlatColumns.EXON_COUNT.ordinal()]);
            final String[] exonStarts = fields[RefFlatColumns.EXON_STARTS.ordinal()].split(",");
            final String[] exonEnds = fields[RefFlatColumns.EXON_ENDS.ordinal()].split(",");
            if (exonCount != exonStarts.length) {
                throw new AnnotationException("Number of exon starts does not agree with number of exons for " +
                        transcriptDescription + " -- ignoring gene");
            }
            if (exonCount != exonEnds.length) {
                throw new AnnotationException("Number of exon ends does not agree with number of exons for " +
                        transcriptDescription + " -- ignoring gene");
            }
            Exon[] exons = new Exon[exonCount];
            for (int i = 0; i < exonCount; ++i) {
                exons[i] = new Exon(Integer.parseInt(exonStarts[i]), Integer.parseInt(exonEnds[i]));
                if (exons[i].start >= exons[i].end) {
                    throw new AnnotationException("Exon has 0 or negative extent for " + transcriptDescription +
                            " -- ignoring gene");
                }
                if (i > 0 && exons[i-1].end >= exons[i].start) {
                    throw new AnnotationException("Exons overlap for " + transcriptDescription +
                            " -- ignoring gene");
                }
            }
            final int transcriptionStart = Integer.parseInt(fields[RefFlatColumns.TX_START.ordinal()]);
            final int transcriptionEnd = Integer.parseInt(fields[RefFlatColumns.TX_END.ordinal()]);
            final int codingStart = Integer.parseInt(fields[RefFlatColumns.CDS_START.ordinal()]);
            final int codingEnd = Integer.parseInt(fields[RefFlatColumns.CDS_END.ordinal()]);
            return new Transcript(transcriptName, transcriptionStart, transcriptionEnd, codingStart, codingEnd, exons);
        }

}
