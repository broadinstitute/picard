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

import net.sf.picard.PicardException;
import net.sf.picard.annotation.Gene;
import net.sf.picard.annotation.Gene.Transcript;
import net.sf.picard.annotation.GeneAnnotationReader;
import net.sf.picard.annotation.LocusFunction;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.*;
import net.sf.samtools.*;
import net.sf.samtools.util.CoordMath;
import net.sf.samtools.util.SequenceUtil;

import java.io.File;
import java.util.*;

public class CollectRnaSeqMetrics extends SinglePassSamProgram {
    private static final Log LOG = Log.getInstance(CollectRnaSeqMetrics.class);

    @Usage
    public final String USAGE = getStandardUsagePreamble() +
            "Program to collect metrics about the alignment of RNA to various functional classes of loci in the genome:" +
            " coding, intronic, UTR, intergenic, ribosomal.\n" +
            "Also determines strand-specificity for strand-specific libraries.";

    public enum StrandSpecificity {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}


    @Option(doc="Gene annotations in refFlat form.  Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat")
    public File REF_FLAT;

    @Option(doc="Location of rRNA sequences in genome, in interval_list format.  " +
            "If not specified no bases will be identified as being ribosomal.  " +
            "Format described here: http://picard.sourceforge.net/javadoc/net/sf/picard/util/IntervalList.html", optional = true)
    public File RIBOSOMAL_INTERVALS;

    @Option(shortName = "STRAND", doc="For strand-specific library prep. " +
            "For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.")
    public StrandSpecificity STRAND_SPECIFICITY;

    @Option(doc="When calculating coverage based values (e.g. CV of coverage) only use transcripts of this length or greater.")
    public int MINIMUM_LENGTH = 500;

    @Option(doc="The PDF file to write out a plot of normalized position vs. coverage.", shortName="CHART", optional = true)
    public File CHART_OUTPUT;

    @Option(doc="If a read maps to a sequence specified with this option, all the bases in the read are counted as ignored bases.  " +
    "These reads are not counted as ")
    public Set<String> IGNORE_SEQUENCE = new HashSet<String>();

    @Option(doc="This percentage of the length of a fragment must overlap one of the ribosomal intervals for a read or read pair by this must in order to be considered rRNA.")
    public double RRNA_FRAGMENT_PERCENTAGE = 0.8;

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    private boolean calculateAll = false;
    private boolean calculateSample = false;
    private boolean calculateLibrary = false;
    private boolean calculateReadGroup = false;

    final RnaSeqMetricsCollector allReadsCollector  = new RnaSeqMetricsCollector(null, null, null, null);
    final Map<String,RnaSeqMetricsCollector> sampleCollectors    = new HashMap<String,RnaSeqMetricsCollector>();
    final Map<String,RnaSeqMetricsCollector> libraryCollectors   = new HashMap<String,RnaSeqMetricsCollector>();
    final Map<String,RnaSeqMetricsCollector> readGroupCollectors = new HashMap<String,RnaSeqMetricsCollector>();

    final private Set<Integer> ignoredSequenceIndices = new HashSet<Integer>();

    private OverlapDetector<Gene> geneOverlapDetector;
    private final OverlapDetector<Interval> ribosomalSequenceOverlapDetector = new OverlapDetector<Interval>(0, 0);

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectRnaSeqMetrics().instanceMainWithExit(argv);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {

        if (CHART_OUTPUT != null) IoUtil.assertFileIsWritable(CHART_OUTPUT);

        calculateAll = METRIC_ACCUMULATION_LEVEL.contains(MetricAccumulationLevel.ALL_READS);
        calculateSample = METRIC_ACCUMULATION_LEVEL.contains(MetricAccumulationLevel.SAMPLE);
        calculateLibrary = METRIC_ACCUMULATION_LEVEL.contains(MetricAccumulationLevel.LIBRARY);
        calculateReadGroup = METRIC_ACCUMULATION_LEVEL.contains(MetricAccumulationLevel.READ_GROUP);

        final Long ribosomalBasesInitialValue = RIBOSOMAL_INTERVALS != null ? 0L : null;
        allReadsCollector.metrics.RIBOSOMAL_BASES = ribosomalBasesInitialValue;
        for (SAMReadGroupRecord rg : header.getReadGroups()) {
            if (calculateSample) {
                if (!sampleCollectors.containsKey(rg.getSample())) {
                    sampleCollectors.put(rg.getSample(),
                        new RnaSeqMetricsCollector(rg.getSample(), null, null, ribosomalBasesInitialValue));
                }
            }
            if (calculateLibrary) {
                if (!libraryCollectors.containsKey(rg.getLibrary())) {
                    libraryCollectors.put(rg.getLibrary(),
                        new RnaSeqMetricsCollector(rg.getSample(), rg.getLibrary(), null, ribosomalBasesInitialValue));
                }
            }
            if (calculateReadGroup) {
                if (!readGroupCollectors.containsKey(rg.getPlatformUnit())) {
                    readGroupCollectors.put(rg.getPlatformUnit(),
                        new RnaSeqMetricsCollector(rg.getSample(), rg.getLibrary(), rg.getPlatformUnit(), ribosomalBasesInitialValue));
                }
            }
        }

        geneOverlapDetector = GeneAnnotationReader.loadRefFlat(REF_FLAT, header.getSequenceDictionary());
        LOG.info("Loaded " + geneOverlapDetector.getAll().size() + " genes.");

        if (RIBOSOMAL_INTERVALS != null) {

            final IntervalList ribosomalIntervals = IntervalList.fromFile(RIBOSOMAL_INTERVALS);
            try {
                SequenceUtil.assertSequenceDictionariesEqual(header.getSequenceDictionary(), ribosomalIntervals.getHeader().getSequenceDictionary());
            } catch (SequenceUtil.SequenceListsDifferException e) {
                throw new PicardException("Sequence dictionaries differ in " + INPUT.getAbsolutePath() + " and " + RIBOSOMAL_INTERVALS.getAbsolutePath(),
                        e);
            }
            ribosomalIntervals.unique();
            final List<Interval> intervals = ribosomalIntervals.getIntervals();
            ribosomalSequenceOverlapDetector.addAll(intervals, intervals);
        }

        for (final String sequenceName: IGNORE_SEQUENCE) {
            final SAMSequenceRecord sequenceRecord = header.getSequence(sequenceName);
            if (sequenceRecord == null) {
                throw new PicardException("Unrecognized sequence " + sequenceName + " passed as argument to IGNORE_SEQUENCE");
            }
            ignoredSequenceIndices.add(sequenceRecord.getSequenceIndex());
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        SAMReadGroupRecord rg = rec.getReadGroup();

        if (calculateAll) {
            allReadsCollector.acceptRead(rec);
        }

        if (calculateSample) {
            sampleCollectors.get(rg.getSample()).acceptRead(rec);
        }

        if (calculateLibrary) {
            libraryCollectors.get(rg.getLibrary()).acceptRead(rec);
        }

        if (calculateReadGroup) {
            readGroupCollectors.get(rg.getPlatformUnit()).acceptRead(rec);
        }
    }

    @Override
    protected void finish() {
        final MetricsFile<RnaSeqMetrics, Integer> file = getMetricsFile();

        if (calculateAll) {
            allReadsCollector.addMetricsToFile(file);
        }

        if (calculateSample) {
            for (RnaSeqMetricsCollector coll : sampleCollectors.values()) {
                coll.addMetricsToFile(file);
            }
        }

        if (calculateLibrary) {
            for (RnaSeqMetricsCollector coll : libraryCollectors.values()) {
                coll.addMetricsToFile(file);
            }
        }

        if (calculateReadGroup) {
            for (RnaSeqMetricsCollector coll : readGroupCollectors.values()) {
                coll.addMetricsToFile(file);
            }
        }
        file.write(OUTPUT);

        boolean atLeastOneHistogram = false;
        for (Histogram<Integer> histo : file.getAllHistograms()) {
            atLeastOneHistogram = atLeastOneHistogram || !histo.isEmpty();
        }
        // Generate the coverage by position plot
        if (CHART_OUTPUT != null && atLeastOneHistogram) {
            final int rResult = RExecutor.executeFromClasspath("net/sf/picard/analysis/rnaSeqCoverage.R",
                                                               OUTPUT.getAbsolutePath(),
                                                               CHART_OUTPUT.getAbsolutePath(),
                                                               INPUT.getName());

            if (rResult != 0) {
                throw new PicardException("Problem invoking R to generate plot.");
            }
        }
    }

    private class RnaSeqMetricsCollector {

        final RnaSeqMetrics metrics = new RnaSeqMetrics();
        private final Map<Transcript, int[]> coverageByTranscript = new HashMap<Transcript, int[]>();

        public RnaSeqMetricsCollector(final String sample,
                                      final String library,
                                      final String readGroup,
                                      final Long ribosomalBasesInitialValue) {
            this.metrics.SAMPLE = sample;
            this.metrics.LIBRARY = library;
            this.metrics.READ_GROUP = readGroup;
            this.metrics.RIBOSOMAL_BASES = ribosomalBasesInitialValue;
        }

        public void acceptRead(SAMRecord rec) {
            // Filter out some reads, and collect the total number of PF bases
            if (rec.getReadFailsVendorQualityCheckFlag() || rec.getNotPrimaryAlignmentFlag()) return;

            this.metrics.PF_BASES += rec.getReadLength();
            if (rec.getReadUnmappedFlag()) return;

            if (ignoredSequenceIndices.contains(rec.getReferenceIndex())) {
                ++this.metrics.IGNORED_READS;
                return;
            }

            // Grab information about the alignment and overlapping genes etc.
            final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());

            // Attempt to get an interval for the entire fragment (if paired read) else just use the read itself.
            // If paired read is chimeric or has one end unmapped, don't create an interval.
            final Interval fragmentInterval;
            if (!rec.getReadPairedFlag()) {
                fragmentInterval = readInterval;
            } else if (rec.getMateUnmappedFlag() || rec.getReferenceIndex() != rec.getMateReferenceIndex()) {
                fragmentInterval = null;
            } else {
                final int fragmentStart = Math.min(rec.getAlignmentStart(), rec.getMateAlignmentStart());
                final int fragmentEnd = CoordMath.getEnd(fragmentStart, Math.abs(rec.getInferredInsertSize()));
                fragmentInterval = new Interval(rec.getReferenceName(), fragmentStart, fragmentEnd);
            }
            if (fragmentInterval != null) {
                final Collection<Interval> overlappingRibosomalIntervals = ribosomalSequenceOverlapDetector.getOverlaps(fragmentInterval);
                int intersectionLength = 0;
                for (final Interval overlappingInterval : overlappingRibosomalIntervals) {
                    final int thisIntersectionLength = overlappingInterval.getIntersectionLength(fragmentInterval);
                    intersectionLength = Math.max(intersectionLength, thisIntersectionLength);
                }
                if (intersectionLength/(double)fragmentInterval.length() >= RRNA_FRAGMENT_PERCENTAGE) {
                    // Assume entire read is ribosomal.
                    // TODO: Should count reads, not bases?
                    metrics.RIBOSOMAL_BASES += rec.getReadLength();
                    int numAlignedBases = 0;
                    for (final AlignmentBlock alignmentBlock : rec.getAlignmentBlocks()) {
                        numAlignedBases += alignmentBlock.getLength();
                    }
                    metrics.PF_ALIGNED_BASES += numAlignedBases;
                    return;
                }
            }

            final Collection<Gene> overlappingGenes                  = geneOverlapDetector.getOverlaps(readInterval);
            final List<AlignmentBlock> alignmentBlocks               = rec.getAlignmentBlocks();
            boolean overlapsExon = false;

            for (final AlignmentBlock alignmentBlock : alignmentBlocks) {
                // Get functional class for each position in the alignment block.
                final LocusFunction[] locusFunctions = new LocusFunction[alignmentBlock.getLength()];

                // By default, if base does not overlap with rRNA or gene, it is intergenic.
                Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);

                for (final Gene gene : overlappingGenes) {
                    for (final Gene.Transcript transcript : gene) {
                        transcript.assignLocusFunctionForRange(alignmentBlock.getReferenceStart(), locusFunctions);

                        // Add coverage to our coverage counter for this transcript
                        int[] coverage = this.coverageByTranscript.get(transcript);
                        if (coverage == null) {
                            coverage = new int[transcript.length()];
                            this.coverageByTranscript.put(transcript, coverage);
                        }
                        transcript.addCoverageCounts(alignmentBlock.getReferenceStart(),
                                                     CoordMath.getEnd(alignmentBlock.getReferenceStart(), alignmentBlock.getLength()),
                                                     coverage);
                    }
                }

                // Tally the function of each base in the alignment block.
                for (final LocusFunction locusFunction : locusFunctions) {
                    ++metrics.PF_ALIGNED_BASES;
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
                final boolean readAndTranscriptStrandsAgree = negativeReadStrand == negativeTranscriptionStrand;
                final boolean readOneOrUnpaired = !rec.getReadPairedFlag() || rec.getFirstOfPairFlag();
                final boolean firstReadExpectedToAgree = STRAND_SPECIFICITY == StrandSpecificity.FIRST_READ_TRANSCRIPTION_STRAND;
                final boolean thisReadExpectedToAgree = readOneOrUnpaired == firstReadExpectedToAgree;
                // If the read strand is the same as the strand of the transcript, and the end is the one that is supposed to agree,
                // then the strand specificity for this read is correct.
                // -- OR --
                // If the read strand is not the same as the strand of the transcript, and the end is not the one that is supposed
                // to agree, then the strand specificity for this read is correct.
                if (readAndTranscriptStrandsAgree == thisReadExpectedToAgree) {
                    ++metrics.CORRECT_STRAND_READS;
                } else {
                    ++metrics.INCORRECT_STRAND_READS;
                }
            }

        }

        public void addMetricsToFile(final MetricsFile<RnaSeqMetrics,Integer> file) {

            if (metrics.PF_ALIGNED_BASES > 0) {
                if (metrics.RIBOSOMAL_BASES != null) {
                    metrics.PCT_RIBOSOMAL_BASES =  metrics.RIBOSOMAL_BASES  / (double) metrics.PF_ALIGNED_BASES;
                }
                metrics.PCT_CODING_BASES =     metrics.CODING_BASES     / (double) metrics.PF_ALIGNED_BASES;
                metrics.PCT_UTR_BASES =        metrics.UTR_BASES        / (double) metrics.PF_ALIGNED_BASES;
                metrics.PCT_INTRONIC_BASES =   metrics.INTRONIC_BASES   / (double) metrics.PF_ALIGNED_BASES;
                metrics.PCT_INTERGENIC_BASES = metrics.INTERGENIC_BASES / (double) metrics.PF_ALIGNED_BASES;
                metrics.PCT_MRNA_BASES =       metrics.PCT_CODING_BASES + metrics.PCT_UTR_BASES;
                metrics.PCT_USABLE_BASES =     (metrics.CODING_BASES + metrics.UTR_BASES) / (double) metrics.PF_BASES;
            }

            if (metrics.CORRECT_STRAND_READS > 0 || metrics.INCORRECT_STRAND_READS > 0) {
                metrics.PCT_CORRECT_STRAND_READS = metrics.CORRECT_STRAND_READS/(double)(metrics.CORRECT_STRAND_READS + metrics.INCORRECT_STRAND_READS);
            }

            // Compute metrics based on coverage of top 1000 genes
            final Histogram<Integer> normalizedCovByPos = computeCoverageMetrics();
            file.addMetric(metrics);
            file.addHistogram(normalizedCovByPos);
        }

        /**
         * Computes a set of coverage based metrics on the mostly highly expressed genes' most highly
         * expressed transcripts.
         */
        private Histogram<Integer> computeCoverageMetrics() {
            final Histogram<Double> cvs = new Histogram<Double>();
            final Histogram<Double> fivePrimeSkews = new Histogram<Double>();
            final Histogram<Double> threePrimeSkews = new Histogram<Double>();
            final Histogram<Double> gapBasesPerKb = new Histogram<Double>();
            final Histogram<Double> fiveToThreeSkews = new Histogram<Double>();
            String prefix = null;
            if (this.metrics.READ_GROUP != null) {
                prefix = this.metrics.READ_GROUP + ".";
            }
            else if (this.metrics.LIBRARY != null) {
                prefix = this.metrics.LIBRARY + ".";
            }
            else if (this.metrics.SAMPLE != null) {
                prefix = this.metrics.SAMPLE + ".";
            }
            else {
                prefix = "All_Reads.";
            }

            final Histogram<Integer> normalizedCoverageByNormalizedPosition = new Histogram<Integer>("normalized_position", prefix + "normalized_coverage");

            final Map<Transcript,int[]> transcripts = pickTranscripts(coverageByTranscript);
            final double transcriptCount = transcripts.size();

            for (final Map.Entry<Transcript,int[]> entry : transcripts.entrySet()) {
                final Transcript tx = entry.getKey();
                final double[] coverage;
                {
                    final double[] tmp = MathUtil.promote(entry.getValue());
                    if (tx.getGene().isPositiveStrand())  coverage = tmp;
                    else coverage = copyAndReverse(tmp);
                }
                final double mean = MathUtil.mean(coverage, 0, coverage.length);

                // Calculate the CV of coverage for this tx
                final double stdev = MathUtil.stddev(coverage, 0, coverage.length, mean);
                final double cv    = stdev / mean;
                cvs.increment(cv);

                // Calculate the 5' and 3' biases
                {
                    final int PRIME_BASES = 100;
                    final double fivePrimeCoverage = MathUtil.mean(coverage, 0, PRIME_BASES);
                    final double threePrimeCoverage = MathUtil.mean(coverage, coverage.length - PRIME_BASES, coverage.length);

                    fivePrimeSkews.increment(fivePrimeCoverage / mean);
                    threePrimeSkews.increment(threePrimeCoverage / mean);
                    fiveToThreeSkews.increment(fivePrimeCoverage / threePrimeCoverage);
                }

                // Calculate normalized coverage vs. normalized position
                {
                    final int lastIndex = coverage.length - 1;

                    for (int percent=0; percent<=100; ++percent) {
                        final double p = percent / 100d;
                        final int start  = (int) Math.max(0,         lastIndex * (p-0.005));
                        final int end    = (int) Math.min(lastIndex, lastIndex * (p+0.005));
                        final int length = end - start + 1;

                        double sum = 0;
                        for (int i=start; i<=end; ++i) sum += coverage[i];
                        final double normalized = (sum / length) / mean;
                        normalizedCoverageByNormalizedPosition.increment(percent, normalized / transcriptCount);
                    }
                }

                // Calculate gap bases per kilobase
//            {
//                int gapBases = 0;
//                final double minCoverage = mean * 0.1;
//                for (int i=0; i<coverage.length; ++i) {
//                    if (coverage[i] < minCoverage) ++gapBases;
//                }
//                gapBasesPerKb.increment(gapBases / (coverage.length / 1000d));
//            }
            }

            this.metrics.MEDIAN_CV_COVERAGE = cvs.getMedian();
            this.metrics.MEDIAN_5PRIME_BIAS = fivePrimeSkews.getMedian();
            this.metrics.MEDIAN_3PRIME_BIAS = threePrimeSkews.getMedian();
            this.metrics.MEDIAN_5PRIME_TO_3PRIME_BIAS = fiveToThreeSkews.getMedian();

            return normalizedCoverageByNormalizedPosition;
        }

        /** Little method to copy an array and reverse it at the same time. */
        private double[] copyAndReverse(final double[] in) {
            final double[] out = new double[in.length];
            for (int i=0, j=in.length-1; i<in.length; ++i, --j) out[j] = in[i];
            return out;
        }

        /** Picks the set of transcripts on which the coverage metrics are to be calculated. */
        public Map<Transcript, int[]> pickTranscripts(final Map<Transcript, int[]> transcriptCoverage) {
            final Map<Transcript, Double> bestPerGene = new HashMap<Transcript, Double>();

            // Make a map of the best transcript per gene to it's mean coverage
            for (final Gene gene : geneOverlapDetector.getAll()) {
                Transcript best = null;
                double bestMean = 0;

                for (final Transcript tx : gene) {
                    final int[] cov = transcriptCoverage.get(tx);

                    if (tx.length() < Math.max(MINIMUM_LENGTH, 100)) continue;
                    if (cov == null) continue;

                    final double mean = MathUtil.mean(MathUtil.promote(cov), 0, cov.length);
                    if (mean < 1d) continue;
                    if (best == null || mean > bestMean) {
                        best = tx;
                        bestMean = mean;
                    }
                }

                if (best != null) bestPerGene.put(best, bestMean);
            }

            // Find the 1000th best coverage value
            final double[] coverages = new double[bestPerGene.size()];
            int i=0;
            for (final double d : bestPerGene.values()) coverages[i++] = d;
            Arrays.sort(coverages);
            final double min = coverages.length == 0 ? 0 : coverages[Math.max(0, coverages.length - 1001)];

            // And finally build the output map
            final Map<Transcript, int[]> retval = new HashMap<Transcript, int[]>();
            for (final Map.Entry<Transcript,Double> entry : bestPerGene.entrySet()) {
                final Transcript tx = entry.getKey();
                final double coverage = entry.getValue();

                if (coverage >= min) {
                    retval.put(tx, transcriptCoverage.get(tx));
                }
            }

            return retval;
        }

    }

}
