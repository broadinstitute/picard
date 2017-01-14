/*
 * The MIT License
 *
 * Copyright (c) 2017 Nils Homer
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

package picard.analysis;


import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import picard.annotation.Gene;
import picard.annotation.LocusFunction;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Optional;

/** A utility class for inferring strand specificity for RNA-Seq experiments. */
public class StrandSpecificityUtil {

    public enum StrandSpecificity {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}

    /**
     * This method can be used to infer the strand specificity of an RNA-Seq experiment.
     *
     * The strand specificity refers to how the sequencing was performed relative to the strand of transcription.  The
     * strand specificity of sequencing is determined from mapped reads while the strand of transcription is determined from
     * a set of genes whose transcription strand and regions are known.
     *
     * There are three possible outcomes:
     * 1. No strand specificity - occurs when there is no correlation between the strands of the reads and the strand of transcription.
     * 2. First read is on the transcription strand.
     * 3. First read is on the opposite strand of transcription.
     * For paired end reads, we assume the second read is on the opposite strand as that of the first (forward/reverse orientation).
     *
     * Vendor QC fail, secondary, supplementary, unmapped, low mapping quality, and duplicate reads will be ignored.  Additionally,
     * only F/R reads and those that have bases that map to UTR or coding sequence will be examined.
     *
     * Reads that overlap multiple genes, where those genes disagree on their transcription strand, cannot be classified, but
     * will count towards the total number of reads examined.
     *
     * @param reader the reader from which records are reads.
     * @param geneOverlapDetector the overlap detector of genes.
     * @param sampleSize the maximum number of reads to sample.
     * @param minimumMappingQuality the minimum mapping quality required to examine a read.
     * @param minFraction the minimum fraction of reads that must agree with a strand model.
     * @param log optionally write the results of the strand specificity to this log.
     * @return the inferred strand specificity.  If it cannot be determined, an Exception is thrown.
     */
    public static StrandSpecificity inferStrandSpecificity(final SamReader reader,
                                                           final OverlapDetector<Gene> geneOverlapDetector,
                                                           final int sampleSize,
                                                           final int minimumMappingQuality,
                                                           final double minFraction,
                                                           final Optional<Log> log
                                                      ) {
        int readOneAgree     = 0;
        int readOneDisagree  = 0;
        int readTwoAgree     = 0;
        int readTwoDisagree  = 0;
        int numReadsExamined = 0;

        for (final SAMRecord rec : reader) {
            if (numReadsExamined >= sampleSize) break;

            // Ignore some reads
            if (rec.getReadFailsVendorQualityCheckFlag() || rec.isSecondaryOrSupplementary()) continue;
            else if (rec.getReadUnmappedFlag() || (rec.getReadPairedFlag() && rec.getMateUnmappedFlag())) continue;
            else if (rec.getMappingQuality() < minimumMappingQuality) continue;
            else if (rec.getDuplicateReadFlag()) continue;
            else if (rec.getReadPairedFlag() && SamPairUtil.getPairOrientation(rec) != SamPairUtil.PairOrientation.FR) continue;

            // Grab information about the alignment and overlapping genes etc.
            final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
            final Collection<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(readInterval);

            // Ignore reads that overlap multiple genes
            if (overlappingGenes.size() == 0) continue;

            // Find if any alignment block overlaps a UTR or coding region in any gene
            final List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
            boolean overlapsExon = alignmentBlocks.stream().anyMatch(alignmentBlock -> {
                // Get functional class for each position in the alignment block.
                final LocusFunction[] locusFunctions = new LocusFunction[alignmentBlock.getLength()];

                // By default, if base does not overlap with rRNA or gene, it is intergenic.
                Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);

                // Assign the function for each transcripts
                for (final Gene gene : overlappingGenes) {
                    for (final Gene.Transcript transcript : gene) {
                        transcript.assignLocusFunctionForRange(alignmentBlock.getReferenceStart(), locusFunctions);
                    }
                }

                // Tally the function of each base in the alignment block.
                for (final LocusFunction locusFunction : locusFunctions) {
                    switch (locusFunction) {
                        case UTR:
                        case CODING:
                            return true;
                        case INTERGENIC:
                        case INTRONIC:
                        case RIBOSOMAL:
                            break;
                    }
                }

                return false;
            });

            // Only examine reads that overlap an exon
            if (overlapsExon) {
                // We would need both read pairs to check multiple genes, and therefore we ignore such cases and call
                // them unclassified later.
                if (overlappingGenes.stream().map(Gene::isPositiveStrand).distinct().count() == 1) {
                    final boolean negativeReadStrand            = rec.getReadNegativeStrandFlag();
                    final boolean negativeTranscriptionStrand   = overlappingGenes.iterator().next().isNegativeStrand();

                    if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                        if (negativeReadStrand == negativeTranscriptionStrand) readOneAgree    += 1;
                        else                                                   readOneDisagree += 1;
                    } else {
                        if (negativeReadStrand == negativeTranscriptionStrand) readTwoAgree    += 1;
                        else                                                   readTwoDisagree += 1;
                    }
                }
                numReadsExamined += 1;
            }
        }

        // NB: these two values may not add up to one if some reads cannot be classified due to overlapping multiple
        // genes that have different transcription strands
        double fractionFirstReadTranscriptionStrand  = (readOneAgree + readTwoDisagree) / (double) numReadsExamined;
        double fractionSecondReadTranscriptionStrand = (readOneDisagree + readTwoAgree) / (double) numReadsExamined;
        double fractionUnclassified = (numReadsExamined - readOneAgree - readOneDisagree - readTwoAgree - readTwoDisagree) / (double) numReadsExamined;

        final int numReads = numReadsExamined; // OK IntelliJ
        log.ifPresent(l -> {
            l.info("Examined " + numReads + " reads to infer strand specificity.");
            l.info("Fraction of reads consistent with read one on the transcription strand:     " + fractionFirstReadTranscriptionStrand);
            l.info("Fraction of reads consistent with read one not on the transcription strand: " + fractionSecondReadTranscriptionStrand);
            l.info("Unclassified reads: " + fractionUnclassified);
        });

        if (fractionFirstReadTranscriptionStrand >= minFraction && fractionSecondReadTranscriptionStrand >= minFraction) {
            throw new IllegalStateException("Could not infer strand specificity as both + ("
                    + fractionFirstReadTranscriptionStrand
                    + ") and - ("
                    + fractionSecondReadTranscriptionStrand
                    + ") are >= the threshold '" + minFraction + "'");
        }
        else if (fractionUnclassified >= minFraction) {
            throw new IllegalStateException("Could not infer strand specificity as the number of unclassified reads ("
                    + fractionUnclassified
                    + ") is >= the threshold '" + minFraction + "'");
        }
        else if (fractionFirstReadTranscriptionStrand >= minFraction) {
            log.ifPresent(l -> l.info("Inferred the first read is on the transcription strand."));
            return StrandSpecificity.FIRST_READ_TRANSCRIPTION_STRAND;
        }
        else if (fractionSecondReadTranscriptionStrand >= minFraction) {
            log.ifPresent(l -> l.info("Inferred the first read is on opposite strand as that of the transcription strand."));
            return StrandSpecificity.SECOND_READ_TRANSCRIPTION_STRAND;
        }
        else {
            log.ifPresent(l -> l.info("Inferred no strand specificity."));
            return StrandSpecificity.NONE;
        }
    }
}
