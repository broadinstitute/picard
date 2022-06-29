/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
import picard.metrics.PerUnitMetricCollector;
import picard.metrics.SAMRecordAndReference;
import picard.metrics.SAMRecordAndReferenceMultiLevelCollector;
import picard.util.MathUtil;

import java.util.List;
import java.util.Set;

public class AlignmentSummaryMetricsCollector extends SAMRecordAndReferenceMultiLevelCollector<AlignmentSummaryMetrics, Integer> {
    // If we have a reference sequence, collect metrics on how well we aligned to it
    private final boolean doRefMetrics;

    //Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.
    private final int maxInsertSize;

    //Paired-end reads that do not have this expected orientation will be considered chimeric.
    private final Set<PairOrientation> expectedOrientations;

    //Whether the SAM or BAM file consists of bisulfite sequenced reads.
    private final boolean isBisulfiteSequenced;

    //The minimum mapping quality a base has to meet in order to be considered high quality
    private final int MAPPING_QUALITY_THRESHOLD = 20;

    //The minimum quality a base has to meet in order to be consider hq_20
    private static final int BASE_QUALITY_THRESHOLD = 20;

    //the adapter utility class
    private final AdapterUtility adapterUtility;

    public AlignmentSummaryMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords,
                                            final boolean doRefMetrics, final List<String> adapterSequence, final int maxInsertSize,
                                            final Set<PairOrientation> expectedOrientations, final boolean isBisulfiteSequenced) {
        this.doRefMetrics = doRefMetrics;
        this.adapterUtility = new AdapterUtility(adapterSequence);
        this.maxInsertSize = maxInsertSize;
        this.expectedOrientations = expectedOrientations;
        this.isBisulfiteSequenced = isBisulfiteSequenced;
        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected PerUnitMetricCollector<AlignmentSummaryMetrics, Integer, SAMRecordAndReference> makeChildCollector(String sample, String library, String readGroup) {
        return new GroupAlignmentSummaryMetricsPerUnitMetricCollector(sample, library, readGroup);
    }

    @Override
    public void acceptRecord(final SAMRecord rec, final ReferenceSequence ref) {
        if (!rec.isSecondaryOrSupplementary()) {
            super.acceptRecord(rec, ref);
        }
    }

    /**
     * returns The sum of lengths of a particular cigar operator in the provided cigar
     *
     * @param cigar The input Cigar of the read
     * @param op    The operator that is being looked for
     * @return Sum of lengths of the Cigar elements in cigar that are of the operator op
     */
    static private int getTotalCigarOperatorCount(final Cigar cigar, final CigarOperator op) {
        return cigar.getCigarElements().stream()
                .filter(e -> e.getOperator().equals(op))
                .mapToInt(CigarElement::getLength)
                .reduce(Integer::sum).orElse(0);
    }


    /**
     * returns the length of the soft clip on the 3' end
     *
     * If there are no-non-clipping operators, method will return 0 as it is unclear which clips should be considered on the
     * "3'" end.
     *
     * @param cigar          The input Cigar of the read
     * @param negativeStrand the negativeStrandFlag of the read
     * @return the amount of soft-clipping that the read has on its 3' end (the later read cycles)
     */
    @VisibleForTesting
    static protected int get3PrimeSoftClippedBases(final Cigar cigar, final boolean negativeStrand) {

        final List<CigarElement> cigarElements;
        if (!negativeStrand) {
            cigarElements = cigar.getCigarElements();
        } else {
            // flip the order of the operators, so that we can always just look for the softclip at the end
            cigarElements = Lists.reverse(cigar.getCigarElements());
        }

        boolean foundNonSoftClipOperator = false;
        int softclipsFound = 0;
        // sum up the last softclips as long they are after non-clipping ops
        for (CigarElement cigarElement : cigarElements) {
            if (!cigarElement.getOperator().isClipping()) {
                foundNonSoftClipOperator = true;
                continue;
            }
            if (foundNonSoftClipOperator && cigarElement.getOperator().equals(CigarOperator.SOFT_CLIP)) {
                softclipsFound += cigarElement.getLength();
            }
        }

        return softclipsFound;
    }

    public class GroupAlignmentSummaryMetricsPerUnitMetricCollector implements PerUnitMetricCollector<AlignmentSummaryMetrics, Integer, SAMRecordAndReference> {
        final IndividualAlignmentSummaryMetricsCollector unpairedCollector;
        final IndividualAlignmentSummaryMetricsCollector firstOfPairCollector;
        final IndividualAlignmentSummaryMetricsCollector secondOfPairCollector;
        final IndividualAlignmentSummaryMetricsCollector pairCollector;
        final String sample;
        final String library;
        final String readGroup;

        public GroupAlignmentSummaryMetricsPerUnitMetricCollector(final String sample, final String library, final String readGroup) {
            this.sample = sample;
            this.library = library;
            this.readGroup = readGroup;
            unpairedCollector = new IndividualAlignmentSummaryMetricsCollector(AlignmentSummaryMetrics.Category.UNPAIRED, sample, library, readGroup);
            firstOfPairCollector = new IndividualAlignmentSummaryMetricsCollector(AlignmentSummaryMetrics.Category.FIRST_OF_PAIR, sample, library, readGroup);
            secondOfPairCollector = new IndividualAlignmentSummaryMetricsCollector(AlignmentSummaryMetrics.Category.SECOND_OF_PAIR, sample, library, readGroup);
            pairCollector = new IndividualAlignmentSummaryMetricsCollector(AlignmentSummaryMetrics.Category.PAIR, sample, library, readGroup);
        }

        public void acceptRecord(final SAMRecordAndReference args) {

            if (args.getSamRecord().getReadPairedFlag()) {
                if (args.getSamRecord().getFirstOfPairFlag()) {
                    firstOfPairCollector.acceptRecord(args);
                } else {
                    secondOfPairCollector.acceptRecord(args);
                }

                pairCollector.acceptRecord(args);
            } else {
                unpairedCollector.acceptRecord(args);
            }
        }

        @Override
        public void finish() {
            // Let the collectors do any summary computations etc.
            unpairedCollector.finish();
            firstOfPairCollector.finish();
            secondOfPairCollector.finish();
            pairCollector.finish();
        }

        @Override
        public void addMetricsToFile(final MetricsFile<AlignmentSummaryMetrics, Integer> file) {
            if (firstOfPairCollector.getMetrics().TOTAL_READS > 0) {
                // override how bad cycle is determined for paired reads, it should be
                // the sum of first and second reads
                pairCollector.getMetrics().BAD_CYCLES = firstOfPairCollector.getMetrics().BAD_CYCLES +
                        secondOfPairCollector.getMetrics().BAD_CYCLES;

                firstOfPairCollector.addMetricsToFile(file);
                secondOfPairCollector.addMetricsToFile(file);
                pairCollector.addMetricsToFile(file);
            }

            // if there are no reads in any category then we will returned an unpaired alignment summary metric with all zero values
            if (unpairedCollector.getMetrics().TOTAL_READS > 0 || firstOfPairCollector.getMetrics().TOTAL_READS == 0) {
                unpairedCollector.addMetricsToFile(file);
            }
        }
    }

    /**
     * Class that counts reads that match various conditions
     */
    public class IndividualAlignmentSummaryMetricsCollector implements PerUnitMetricCollector<AlignmentSummaryMetrics, Integer, SAMRecordAndReference> {
        private long numPositiveStrand;
        private final Histogram<Integer> readLengthHistogram = new Histogram<>("count", "readLength");
        private final Histogram<Integer> alignedReadLengthHistogram = new Histogram<>("count", "alignedReadLength");

        private final AlignmentSummaryMetrics metrics;
        private long chimeras;
        private long chimerasDenominator;
        private long adapterReads;
        private long indels;

        private long numSoftClipped;
        private long num3PrimeSoftClippedBases;
        private long numReadsWith3PrimeSoftClips;

        private long numHardClipped;

        private long nonBisulfiteAlignedBases;
        private long hqNonBisulfiteAlignedBases;
        private final Histogram<Long> mismatchHistogram = new Histogram<>();
        private final Histogram<Long> hqMismatchHistogram = new Histogram<>();
        private final Histogram<Integer> badCycleHistogram = new Histogram<>();

        public IndividualAlignmentSummaryMetricsCollector(final AlignmentSummaryMetrics.Category pairingCategory,
                                                          final String sample,
                                                          final String library,
                                                          final String readGroup) {
            metrics = new AlignmentSummaryMetrics();
            metrics.CATEGORY = pairingCategory;
            metrics.SAMPLE = sample;
            metrics.LIBRARY = library;
            metrics.READ_GROUP = readGroup;
        }

        public void acceptRecord(final SAMRecordAndReference samRecordAndReference) {
            final SAMRecord record = samRecordAndReference.getSamRecord();
            final ReferenceSequence ref = samRecordAndReference.getReferenceSequence();

            if (record.isSecondaryAlignment()) {
                // only want 1 count per read so skip non-primary alignments
                return;
            }

            collectReadData(record);
            collectQualityData(record, ref);
        }

        @Override
        public void finish() {
            //summarize read data
            if (metrics.TOTAL_READS > 0) {
                metrics.PCT_PF_READS = (double) metrics.PF_READS / (double) metrics.TOTAL_READS;
                metrics.PCT_ADAPTER = adapterReads / (double) metrics.PF_READS;
                metrics.MEAN_READ_LENGTH = readLengthHistogram.getMean();
                metrics.SD_READ_LENGTH = readLengthHistogram.getStandardDeviation();
                metrics.MEDIAN_READ_LENGTH = readLengthHistogram.getMedian();
                metrics.MAD_READ_LENGTH = readLengthHistogram.getMedianAbsoluteDeviation();
                metrics.MIN_READ_LENGTH = readLengthHistogram.getMin();
                metrics.MAX_READ_LENGTH = readLengthHistogram.getMax();

                //Calculate BAD_CYCLES
                metrics.BAD_CYCLES = 0;
                for (final Histogram.Bin<Integer> cycleBin : badCycleHistogram.values()) {
                    final double badCyclePercentage = cycleBin.getValue() / metrics.TOTAL_READS;
                    if (badCyclePercentage >= 0.8) {
                        metrics.BAD_CYCLES++;
                    }
                }

                if (doRefMetrics) {
                    final double totalBases = readLengthHistogram.getSum();
                    metrics.PCT_PF_READS_ALIGNED = MathUtil.divide(metrics.PF_READS_ALIGNED, (double) metrics.PF_READS);
                    metrics.PCT_READS_ALIGNED_IN_PAIRS = MathUtil.divide(metrics.READS_ALIGNED_IN_PAIRS, (double) metrics.PF_READS_ALIGNED);
                    metrics.PCT_PF_READS_IMPROPER_PAIRS = MathUtil.divide(metrics.PF_READS_IMPROPER_PAIRS, (double) metrics.PF_READS_ALIGNED);
                    metrics.MEAN_ALIGNED_READ_LENGTH = alignedReadLengthHistogram.getMean();
                    metrics.STRAND_BALANCE = MathUtil.divide(numPositiveStrand, (double) metrics.PF_READS_ALIGNED);
                    metrics.PCT_CHIMERAS = MathUtil.divide(chimeras, (double) chimerasDenominator);
                    metrics.PF_INDEL_RATE = MathUtil.divide(indels, (double) metrics.PF_ALIGNED_BASES);
                    metrics.PF_MISMATCH_RATE = MathUtil.divide(mismatchHistogram.getSum(), (double) nonBisulfiteAlignedBases);
                    metrics.PF_HQ_ERROR_RATE = MathUtil.divide(hqMismatchHistogram.getSum(), (double) hqNonBisulfiteAlignedBases);

                    metrics.PCT_HARDCLIP = MathUtil.divide(numHardClipped, totalBases);
                    metrics.PCT_SOFTCLIP = MathUtil.divide(numSoftClipped, totalBases);
                    metrics.AVG_POS_3PRIME_SOFTCLIP_LENGTH = MathUtil.divide(num3PrimeSoftClippedBases, (double) numReadsWith3PrimeSoftClips);

                    metrics.PF_HQ_MEDIAN_MISMATCHES = hqMismatchHistogram.getMedian();
                }
            }
        }

        @Override
        public void addMetricsToFile(final MetricsFile<AlignmentSummaryMetrics, Integer> file) {
            file.addMetric(metrics);
        }

        /**
         * returns The number of read bases that are not clipped, from the cigar
         *
         * @param cigar The input Cigar of the read
         * @return Number of read bases that are not clipped
         */
        private int getUnclippedBaseCount(final Cigar cigar) {
            return cigar.getCigarElements().stream()
                    .filter(e -> e.getOperator().consumesReadBases())
                    .filter(e -> !e.getOperator().isClipping())
                    .mapToInt(CigarElement::getLength)
                    .reduce(Integer::sum).orElse(0);
        }


        private void collectReadData(final SAMRecord record) {
            // NB: for read count metrics, do not include supplementary records, but for base count metrics, do include supplementary records.
            if (record.getSupplementaryAlignmentFlag()) {
                return;
            }

            metrics.TOTAL_READS++;

            if (!record.getReadFailsVendorQualityCheckFlag()) {
                metrics.PF_READS++;
                if (isNoiseRead(record)) {
                    metrics.PF_NOISE_READS++;
                }

                readLengthHistogram.increment(record.getReadBases().length);
                alignedReadLengthHistogram.increment(getUnclippedBaseCount(record.getCigar()));

                // See if the read is an adapter sequence
                if (adapterUtility.isAdapter(record)) {
                    adapterReads++;
                }
                // count clipped bases
                numHardClipped += getTotalCigarOperatorCount(record.getCigar(), CigarOperator.HARD_CLIP);

                if (!record.getReadUnmappedFlag()) {
                    numSoftClipped += getTotalCigarOperatorCount(record.getCigar(), CigarOperator.SOFT_CLIP);

                    final int threePrimeSoftClippedBases = get3PrimeSoftClippedBases(record.getCigar(), record.getReadNegativeStrandFlag());
                    if (threePrimeSoftClippedBases > 0) {
                        num3PrimeSoftClippedBases += threePrimeSoftClippedBases;
                        numReadsWith3PrimeSoftClips++;
                    }
                    if (doRefMetrics) {

                        metrics.PF_READS_ALIGNED++;
                        if (record.getReadPairedFlag() && !record.getProperPairFlag()) {
                            metrics.PF_READS_IMPROPER_PAIRS++;
                        }
                        if (!record.getReadNegativeStrandFlag()) {
                            numPositiveStrand++;
                        }
                        if (record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {
                            metrics.READS_ALIGNED_IN_PAIRS++;

                            // Check that both ends have mapq > minimum
                            final Integer mateMq = record.getIntegerAttribute(SAMTag.MQ.toString());
                            if (mateMq == null || mateMq >= MAPPING_QUALITY_THRESHOLD && record.getMappingQuality() >= MAPPING_QUALITY_THRESHOLD) {
                                ++chimerasDenominator;

                                // With both reads mapped we can see if this pair is chimeric
                                if (ChimeraUtil.isChimeric(record, maxInsertSize, expectedOrientations)) {
                                    ++chimeras;
                                }
                            }
                        } else { // fragment reads or read pairs with one end that maps
                            // Consider chimeras that occur *within* the read using the SA tag
                            if (record.getMappingQuality() >= MAPPING_QUALITY_THRESHOLD) {
                                ++chimerasDenominator;
                                if (record.getAttribute(SAMTag.SA.toString()) != null) {
                                    ++chimeras;
                                }
                            }
                        }
                    }
                }
            }
        }

        private void collectQualityData(final SAMRecord record, final ReferenceSequence reference) {
            // NB: for read count metrics, do not include supplementary records, but for base count metrics, do include supplementary records.

            // If the read isn't an aligned PF read then look at the read for no-calls
            if (record.getReadUnmappedFlag() || record.getReadFailsVendorQualityCheckFlag() || !doRefMetrics) {
                final byte[] readBases = record.getReadBases();
                for (int i = 0; i < readBases.length; i++) {
                    if (SequenceUtil.isNoCall(readBases[i])) {
                        badCycleHistogram.increment(CoordMath.getCycle(record.getReadNegativeStrandFlag(), readBases.length, i));
                    }
                }
            } else if (!record.getReadFailsVendorQualityCheckFlag()) {
                final boolean highQualityMapping = isHighQualityMapping(record);
                if (highQualityMapping && !record.getSupplementaryAlignmentFlag()) {
                    metrics.PF_HQ_ALIGNED_READS++;
                }

                final byte[] readBases = record.getReadBases();
                final byte[] refBases = reference == null ? null : reference.getBases();
                final int refLength = reference == null ? Integer.MAX_VALUE : refBases.length;
                final byte[] qualities = record.getBaseQualities();
                long mismatchCount = 0;
                long hqMismatchCount = 0;

                for (final AlignmentBlock alignmentBlock : record.getAlignmentBlocks()) {
                    final int readIndex = alignmentBlock.getReadStart() - 1;
                    final int refIndex = alignmentBlock.getReferenceStart() - 1;
                    final int length = alignmentBlock.getLength();

                    for (int i = 0; i < length && refIndex + i < refLength; ++i) {
                        final int readBaseIndex = readIndex + i;
                        boolean mismatch = refBases != null && !SequenceUtil.basesEqual(readBases[readBaseIndex], refBases[refIndex + i]);

                        final boolean bisulfiteMatch = refBases != null && isBisulfiteSequenced && SequenceUtil.bisulfiteBasesEqual(record.getReadNegativeStrandFlag(), readBases[readBaseIndex], refBases[readBaseIndex]);

                        final boolean bisulfiteBase = mismatch && bisulfiteMatch;
                        mismatch = mismatch && !bisulfiteMatch;

                        if (mismatch) {
                            mismatchCount++;
                        }

                        metrics.PF_ALIGNED_BASES++;
                        if (!bisulfiteBase) {
                            nonBisulfiteAlignedBases++;
                        }

                        if (highQualityMapping) {
                            metrics.PF_HQ_ALIGNED_BASES++;
                            if (!bisulfiteBase) {
                                hqNonBisulfiteAlignedBases++;
                            }
                            if (qualities[readBaseIndex] >= BASE_QUALITY_THRESHOLD) {
                                metrics.PF_HQ_ALIGNED_Q20_BASES++;
                            }
                            if (mismatch) {
                                hqMismatchCount++;
                            }
                        }

                        if (mismatch || SequenceUtil.isNoCall(readBases[readBaseIndex])) {
                            badCycleHistogram.increment(CoordMath.getCycle(record.getReadNegativeStrandFlag(), readBases.length, i));
                        }
                    }
                }

                mismatchHistogram.increment(mismatchCount);
                hqMismatchHistogram.increment(hqMismatchCount);


                // Add any insertions and/or deletions to the global count
                for (final CigarElement elem : record.getCigar().getCigarElements()) {
                    final CigarOperator op = elem.getOperator();
                    if (op == CigarOperator.INSERTION || op == CigarOperator.DELETION) {
                        ++indels;
                    }
                }
            }
        }


        private boolean isNoiseRead(final SAMRecord record) {
            final Object noiseAttribute = record.getAttribute(ReservedTagConstants.XN);
            return (noiseAttribute != null && noiseAttribute.equals(1));
        }

        private boolean isHighQualityMapping(final SAMRecord record) {
            return !record.getReadFailsVendorQualityCheckFlag() &&
                    record.getMappingQuality() >= MAPPING_QUALITY_THRESHOLD;
        }

        public AlignmentSummaryMetrics getMetrics() {
            return metrics;
        }

        public Histogram<Integer> getReadHistogram() {
            return readLengthHistogram;
        }

        public Histogram<Integer> getAlignedReadHistogram() {
            return alignedReadLengthHistogram;
        }
    }
}
