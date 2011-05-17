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
package net.sf.picard.analysis;

import net.sf.picard.annotation.Gene;
import net.sf.picard.annotation.GeneAnnotationReader;
import net.sf.picard.annotation.LocusFunction;
import net.sf.picard.annotation.Transcript;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.Log;
import net.sf.picard.util.OverlapDetector;
import net.sf.samtools.*;
import net.sf.samtools.util.CoordMath;

import java.io.File;
import java.util.Collection;
import java.util.List;

public class CollectRnaSeqMetrics extends SinglePassSamProgram {
    private static final Log LOG = Log.getInstance(CollectRnaSeqMetrics.class);

    @Usage
    public final String USAGE = getStandardUsagePreamble() +
            "Program to collect metrics about the alignment of RNA to the various functional classes of loci in the genome:" +
            " coding, intronic, UTR, intragenic. ribosomal.\n" +
            "Also determines strand-specificity for strand-specific libraries.";

    public enum StrandSpecificity {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}


    @Option(doc="Gene annotations in refFlat form.")
    public File REF_FLAT;

    @Option(doc="Location of rRNA sequences in genome, in interval_list format.  " +
            "If not specified no bases will be identified as being ribosomal.", optional = true)
    public File RIBOSOMAL_INTERVALS;

    @Option(shortName = "STRAND", doc="For strand-specific library prep. " +
            "For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.")
    public StrandSpecificity STRAND_SPECIFICITY;

    private OverlapDetector<Gene> geneOverlapDetector;
    private final OverlapDetector<Interval> ribosomalSequenceOverlapDetector = new OverlapDetector<Interval>(0, 0);

    private RnaSeqMetrics metrics;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectRnaSeqMetrics().instanceMainWithExit(argv);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        geneOverlapDetector = GeneAnnotationReader.loadRefFlat(REF_FLAT, header.getSequenceDictionary());
        LOG.info("Loaded " + geneOverlapDetector.getAll().size() + " genes.");
        if (RIBOSOMAL_INTERVALS != null) {
            final IntervalList ribosomalIntervals = IntervalList.fromFile(RIBOSOMAL_INTERVALS);
            ribosomalIntervals.unique();
            final List<Interval> intervals = ribosomalIntervals.getIntervals();
            ribosomalSequenceOverlapDetector.addAll(intervals, intervals);
        }
        metrics = new RnaSeqMetrics();
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        if (rec.getReadUnmappedFlag() || rec.getReadFailsVendorQualityCheckFlag() || rec.getNotPrimaryAlignmentFlag()) return;
        final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
        final Collection<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(readInterval);
        final Collection<Interval> overlappingRibosomalIntervals = ribosomalSequenceOverlapDetector.getOverlaps(readInterval);
        final List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
        boolean overlapsExon = false;
        for (final AlignmentBlock alignmentBlock : alignmentBlocks) {
            // Get functional class for each position in the alignment block.
            final LocusFunction[] locusFunctions = new LocusFunction[alignmentBlock.getLength()];

            // By default, if base does not overlap with rRNA or gene, it is intergenic.
            for (int i = 0; i < locusFunctions.length; ++i) locusFunctions[i] = LocusFunction.INTERGENIC;

            for (final Gene gene : overlappingGenes) {
                for (final Transcript transcript : gene) {
                    transcript.getLocusFunctionForRange(alignmentBlock.getReferenceStart(), locusFunctions);
                }
            }

            // Do this last because ribosomal takes precedence over the other functions.
            for (final Interval ribosomalInterval : overlappingRibosomalIntervals) {
                assignValueForOverlappingRange(ribosomalInterval, alignmentBlock.getReferenceStart(),
                        locusFunctions, LocusFunction.RIBOSOMAL);
            }

            // Tally the function of each base in the alignment block.
            for (final LocusFunction locusFunction : locusFunctions) {
                ++metrics.ALIGNED_PF_BASES;
                switch (locusFunction) {
                    case INTERGENIC:
                        ++metrics.INTERGENIC_BASES;
                        break;
                    case INTRONIC:
                        ++metrics.INTRONIC_BASES;
                        break;
                    case UTR:
                        ++metrics.UTR_BASES;
                        overlapsExon = true;
                        break;
                    case CODING:
                        ++metrics.CODING_BASES;
                        overlapsExon = true;
                        break;
                    case RIBOSOMAL:
                        ++metrics.RIBOSOMAL_BASES;
                        break;
                }
            }
        }

        // Strand-specificity is tallied on read basis rather than base at a time.  A read that aligns to more than one
        // gene is not counted.
        if (overlapsExon && STRAND_SPECIFICITY != StrandSpecificity.NONE && overlappingGenes.size() == 1) {
            final boolean negativeTranscriptionStrand = overlappingGenes.iterator().next().isNegativeStrand();
            final boolean negativeReadStrand = rec.getReadNegativeStrandFlag();
            final boolean readOneOrUnpaired = !rec.getReadPairedFlag() || rec.getFirstOfPairFlag();
            // If the read strand is the same as the strand of the transcript, and the end is the one that is supposed to agree,
            // then the strand specificity for this read is correct.
            // -- OR --
            // If the read strand is not the same as the strand of the transcript, and the end is not the one that is supposed
            // to agree, then the strand specificity for this read is correct.
            if ((negativeReadStrand == negativeTranscriptionStrand) ==
                (readOneOrUnpaired == (STRAND_SPECIFICITY == StrandSpecificity.FIRST_READ_TRANSCRIPTION_STRAND))) {
                ++metrics.CORRECT_STRAND_READS;
            } else {
                ++metrics.INCORRECT_STRAND_READS;
            }
        }
    }

    @Override
    protected void finish() {
        if (metrics.ALIGNED_PF_BASES > 0) {
            metrics.PCT_RIBOSOMAL_BASES =  metrics.RIBOSOMAL_BASES  / (double) metrics.ALIGNED_PF_BASES;
            metrics.PCT_CODING_BASES =     metrics.CODING_BASES     / (double) metrics.ALIGNED_PF_BASES;
            metrics.PCT_UTR_BASES =        metrics.UTR_BASES        / (double) metrics.ALIGNED_PF_BASES;
            metrics.PCT_INTRONIC_BASES =   metrics.INTRONIC_BASES   / (double) metrics.ALIGNED_PF_BASES;
            metrics.PCT_INTERGENIC_BASES = metrics.INTERGENIC_BASES / (double) metrics.ALIGNED_PF_BASES;
        }
        if (metrics.CORRECT_STRAND_READS > 0 || metrics.INCORRECT_STRAND_READS > 0) {
            metrics.PCT_CORRECT_STRAND_READS = metrics.CORRECT_STRAND_READS/(double)(metrics.CORRECT_STRAND_READS + metrics.INCORRECT_STRAND_READS);
        }
        final MetricsFile<RnaSeqMetrics, Integer> file = getMetricsFile();
        file.addMetric(metrics);
        file.write(OUTPUT);
    }

    private void assignValueForOverlappingRange(final Interval interval, final int start,
                                                final LocusFunction[] locusFunctions, final LocusFunction valueToAssign) {
        for (int i = Math.max(start, interval.getStart());
                i <= Math.min(interval.getEnd(), CoordMath.getEnd(start, locusFunctions.length)); ++i) {
            locusFunctions[i - start] = valueToAssign;
        }
    }
}
