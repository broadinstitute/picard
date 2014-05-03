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

import net.sf.picard.util.Interval;
import net.sf.samtools.util.CoordMath;
import sun.tools.tree.LengthExpression;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Holds annotation of a gene for storage in an OverlapDetector.  May hold multiple transcripts for the same gene.
 * The transcripts must all be relative to the same strand.
 */
public class Gene extends Interval implements Iterable<Gene.Transcript>  {
    private final Map<String, Transcript> transcripts = new HashMap<String, Transcript>();


    public Gene(final String sequence, final int start, final int end, final boolean negative, final String name) {
        super(sequence, start, end, negative, name);
    }

    public Transcript addTranscript(final String name, final int transcriptionStart, final int transcriptionEnd, final int codingStart, final int codingEnd, final int numExons) {
        if (transcripts.containsKey(name)) {
            throw new AnnotationException("Transcript " + name + " for gene " + this.getName() + " appears more than once");
        }
        else {
            final Transcript tx = new Transcript(name, transcriptionStart, transcriptionEnd, codingStart, codingEnd, numExons);
            transcripts.put(name, tx);
            return tx;
        }
    }

    public Iterator<Transcript> iterator() {
        return transcripts.values().iterator();
    }

    /**
     * A single transcript of a gene.  Sequence name is stored in the enclosing object (class Gene).
     */
    public class Transcript {
        public final String name;
        public final int transcriptionStart;
        public final int transcriptionEnd;
        public final int codingStart;
        public final int codingEnd;
        public final Exon[] exons;
        private int length; // the number of bases in the transcript

        /**
         * 1-based, inclusive representation of an exon.  The sequence name is stored in an enclosing object (class Gene).
         */
        public class Exon {
            public final int start;
            public final int end;

            public Exon(final int start, final int end) {
                this.start = start;
                this.end = end;
            }
        }

        public Transcript(final String name, final int transcriptionStart, final int transcriptionEnd, final int codingStart, final int codingEnd, final int numExons) {
            this.name = name;
            this.transcriptionStart = transcriptionStart;
            this.transcriptionEnd = transcriptionEnd;
            this.codingStart = codingStart;
            this.codingEnd = codingEnd;
            this.exons = new Exon[numExons];
        }

        public Exon addExon(final int start, final int end) {
            for (int i=0; i<this.exons.length; ++i) {
                if (exons[i] == null) {
                    exons[i] = new Exon(start, end);
                    this.length += CoordMath.getLength(start, end);
                    return exons[i];
                }
            }

            throw new IllegalStateException("Attempting to add more exons that exist for transcript.");
        }

        public int start() {
            return exons[0].start;
        }

        public int end() {
            return exons[exons.length -1].end;
        }

        public int length() {
            return this.length;
        }

        public boolean isSoloTranscript() {
            return Gene.this.transcripts.size() == 1;
        }

        public Gene getGene() {
            return Gene.this;
        }

        /**
         * Write into locusFunctions the function of each position from start to start + locusFunctions.length
         * relative to this transcript.  Does not overwrite an existing value in locusFunctions that is stronger
         * than the function for that locus in this transcript.
         * @param start 1-based genomic coordinate of the first position in locusFunctions.
         * @param locusFunctions
         */
        public void assignLocusFunctionForRange(final int start, final LocusFunction[] locusFunctions) {
            for (int i = Math.max(start, transcriptionStart);
                    i <= Math.min(transcriptionEnd, CoordMath.getEnd(start, locusFunctions.length)); ++i) {

                if (locusFunctions[i - start].ordinal() > LocusFunction.CODING.ordinal()) continue;

                final LocusFunction locusFunction;
                if (inExon(i)) {
                    if (utr(i)) locusFunction = LocusFunction.UTR;
                    else locusFunction = LocusFunction.CODING;
                } else locusFunction = LocusFunction.INTRONIC;
                if (locusFunction.ordinal() > locusFunctions[i - start].ordinal()) {
                    locusFunctions[i - start] = locusFunction;
                }
            }
        }

        /**
         *
         * @param genomeStart
         * @param genomeEnd
         * @param coverage
         */
        public void addCoverageCounts(final int genomeStart, final int genomeEnd, final int[] coverage) {
            for (int i=genomeStart; i<genomeEnd; ++i) {
                final int txBase = getTranscriptCoordinate(i);
                if (txBase > 0) coverage[txBase-1]++;
            }
        }

        /** Given a coordinate on the genome (same chromosome) give the corresponding coordinate in the transcript. */
        public int getTranscriptCoordinate(final int genomeCoordinate) {
            int exonOffset = 0;
            for (final Exon e : exons) {
                if (genomeCoordinate >= e.start && genomeCoordinate <=e.end) {
                    return (genomeCoordinate - e.start + 1) + exonOffset;
                }
                else {
                    exonOffset += CoordMath.getLength(e.start, e.end);
                }
            }

            return -1;
        }

        private boolean utr(final int locus) {
            return locus < codingStart || locus > codingEnd;
        }

        private boolean inExon(final int locus) {
            for (int i = 0; i < exons.length; ++i) {
                final Exon exon = exons[i];
                if (exon.start > locus) return false;
                if (inRange(exon.start, exon.end, locus)) return true;
            }
            return false;
        }

        private boolean inRange(final int start, final int end, final int locus) {
            return (locus >= start && locus <= end);
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final Transcript that = (Transcript) o;

            if (codingEnd != that.codingEnd) return false;
            if (codingStart != that.codingStart) return false;
            if (transcriptionEnd != that.transcriptionEnd) return false;
            if (transcriptionStart != that.transcriptionStart) return false;
            if (!name.equals(that.name)) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = name.hashCode();
            result = 31 * result + transcriptionStart;
            result = 31 * result + transcriptionEnd;
            result = 31 * result + codingStart;
            result = 31 * result + codingEnd;
            return result;
        }
    }
}
