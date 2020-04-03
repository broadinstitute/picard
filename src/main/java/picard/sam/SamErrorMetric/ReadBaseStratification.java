/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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
package picard.sam.SamErrorMetric;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.SamLocusAndReferenceIterator.SAMLocusAndReference;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.sam.util.Pair;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.PhysicalLocationInt;

import java.util.Collection;
import java.util.LinkedList;
import java.util.concurrent.ExecutionException;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Stream;

/**
 * Classes, methods, and enums that deal with the stratification of read bases and reference information.
 */
public class ReadBaseStratification {

    private static final Log log = Log.getInstance(CollectSamErrorMetrics.class);

    // This variable has to be set using setLongHomopolymer _before_ the first call to binnedHomopolymerStratifier.get()
    // if you need different binned Homopolymers with different LONG_HOMOPOLYMER values, you'll need to make your own
    // instances.
    private static int LONG_HOMOPOLYMER = 6;
    private static int GC_CACHE_SIZE = 1000;


    /* ***** SETTERS ********** */

    public static void setGcCacheSize(int gcCacheSize) {
        GC_CACHE_SIZE = gcCacheSize;
    }

    /**
     * defaults to 6
     **/
    public static void setLongHomopolymer(int longHomopolymer) {
        LONG_HOMOPOLYMER = longHomopolymer;
    }

    /* ******* general-use classes, for defining and creating new stratifiers ***********/

    /**
     * The main interface for a stratifier.
     * The Stratifier needs to be able to ake a {@link htsjdk.samtools.util.SamLocusIterator.RecordAndOffset RecordAndOffset}
     * and return a value of type T. It also needs to have a suffix for automatically generating metrics file suffixes.
     */
    public interface RecordAndOffsetStratifier<T extends Comparable<T>> {
        // The method that stratifies a base in a read (provided by RecordAndOffset) into a type T.
        T stratify(final RecordAndOffset recordAndOffset, final SAMLocusAndReference locusInfo);

        // The string suffix that will be used to generate the extension of the metric file.
        String getSuffix();
    }

    /**
     * A simpler stratifier for cases when only the record suffices
     */
    abstract static class RecordStratifier<T extends Comparable<T>> implements RecordAndOffsetStratifier<T> {
        @Override
        public T stratify(RecordAndOffset recordAndOffset, SAMLocusAndReference locusInfo) {
            return stratify(recordAndOffset.getRecord());
        }

        abstract T stratify(final SAMRecord sam);
    }

    /**
     * A factory for generating stateless stratifier instances given a static function and a string
     *
     * @param staticStratify the (static) function (from RecordAndOffset and SAMLocusAndReference) to stratify with
     * @param suffix         the suffix to use
     * @param <T>            the type into which the stratification happens
     * @return an instance of a stratification class that will stratify accordingly
     */
    private static <T extends Comparable<T>> RecordAndOffsetStratifier<T> wrapStaticFunction(BiFunction<RecordAndOffset, SAMLocusAndReference, T> staticStratify, String suffix) {
        return new RecordAndOffsetStratifier<T>() {
            @Override
            public T stratify(RecordAndOffset recordAndOffset, SAMLocusAndReference locusInfo) {
                return staticStratify.apply(recordAndOffset, locusInfo);
            }

            @Override
            public String getSuffix() {
                return suffix;
            }
        };
    }

    /**
     * A factory for generating stateless stratifier instances given a static function and a string
     *
     * @param staticStratify the (static) function (from RecordAndOffset) to stratify with
     * @param suffix         the suffix to use
     * @param <T>            the type into which the stratification happens
     * @return an instance of a stratification class that will stratify accordingly
     */
    private static <T extends Comparable<T>> RecordAndOffsetStratifier<T> wrapStaticFunction(Function<RecordAndOffset, T> staticStratify, String suffix) {
        return wrapStaticFunction((rao, ignored) -> staticStratify.apply(rao), suffix);
    }

    /**
     * A simpler factory for generating stateless stratifier instances given a static function and a string
     *
     * @param staticStratify the (static) function (from SAMRecord alone) to stratify with
     * @param suffix         the suffix to us
     * @param <T>            the type into which the stratification happens
     * @return an instance of a stratification class that will stratify accordingly
     */
    private static <T extends Comparable<T>> RecordStratifier<T> wrapStaticReadFunction(Function<SAMRecord, T> staticStratify, String suffix) {
        return new RecordStratifier<T>() {
            @Override
            public T stratify(final SAMRecord sam) {
                return staticStratify.apply(sam);
            }

            @Override
            public String getSuffix() {
                return suffix;
            }
        };
    }

    /**
     * A PairStratifier is a stratifier that uses two other stratifiers to inform the stratification.
     * For a given input, the result is the {@link Pair} of outputs that the two stratifiers return.
     * The suffix for this stratifier is generated from the suffixes of the two provided stratifiers.
     */
    public static class PairStratifier<T extends Comparable<T>, R extends Comparable<R>> implements RecordAndOffsetStratifier<Pair<T, R>> {
        public PairStratifier(final RecordAndOffsetStratifier<T> a, final RecordAndOffsetStratifier<R> b) {
            this.a = a;
            this.b = b;
        }

        final RecordAndOffsetStratifier<T> a;
        final RecordAndOffsetStratifier<R> b;

        @Override
        public Pair<T, R> stratify(final RecordAndOffset recordAndOffset, final SAMLocusAndReference locusInfo) {

            final T a = this.a.stratify(recordAndOffset, locusInfo);
            final R b = this.b.stratify(recordAndOffset, locusInfo);
            if (a == null || b == null) {
                return null;
            }

            return new Pair<>(a, b);
        }

        @Override
        public String getSuffix() {
            return a.getSuffix() + "_and_" + b.getSuffix();
        }
    }

    /**
     * A CollectionStratifier is a stratifier that uses a collection of stratifiers to inform the stratification.
     * For a given input, the result is the "Tuple" (actually repeated pairs) of outputs that the stratifiers return.
     * The suffix for this stratifier is generated from the suffixes of the provided stratifiers.
     */
    public static class CollectionStratifier implements RecordAndOffsetStratifier {
        public CollectionStratifier(final Collection<RecordAndOffsetStratifier<?>> stratifiers) {

            if (stratifiers.isEmpty()) {
                throw new IllegalArgumentException("Must construct with a non-empty collection of stratifiers.");
            }

            final LinkedList<RecordAndOffsetStratifier<?>> linkedListStratifiers = new LinkedList<>(stratifiers);

            while (linkedListStratifiers.size() > 1) {
                ReadBaseStratification.RecordAndOffsetStratifier<?> first = linkedListStratifiers.remove(0);
                ReadBaseStratification.RecordAndOffsetStratifier<?> second = linkedListStratifiers.remove(0);
                ReadBaseStratification.PairStratifier<? extends Comparable, ? extends Comparable> newPair = new ReadBaseStratification.PairStratifier<>(first, second);
                linkedListStratifiers.add(0, newPair);
            }

            stratifier = linkedListStratifiers.remove();
        }

        final RecordAndOffsetStratifier stratifier;

        @Override
        public Comparable stratify(final RecordAndOffset recordAndOffset, final SAMLocusAndReference locusInfo) {

            return stratifier.stratify(recordAndOffset, locusInfo);
        }

        @Override
        public String getSuffix() {
            return stratifier.getSuffix();
        }
    }

    /**
     * A factory for generating "pair" stratifier instances from two stratifiers and a string
     *
     * @param <T>             the type of the left Stratifier
     * @param <S>             the type of the right Stratifier
     * @param leftStratifier  a {@link RecordAndOffsetStratifier} to use
     * @param rightStratifier a {@link RecordAndOffsetStratifier} to use
     * @param suffix          the suffix to use for the new stratifier
     * @return an instance of {@link  PairStratifier} that will stratify according to both <code>leftStratifier</code>
     * and <code>rightStratifier</code>
     */
    public static <T extends Comparable<T>, S extends Comparable<S>> PairStratifier<T, S> PairStratifierFactory(
            final RecordAndOffsetStratifier<T> leftStratifier,
            final RecordAndOffsetStratifier<S> rightStratifier,
            final String suffix) {
        return new PairStratifier<T, S>(leftStratifier, rightStratifier) {
            @Override
            public String getSuffix() {
                return suffix;
            }
        };
    }


    /* ****************  Actual Stratification Classes *****************/

    /**
     * A stratifier that uses GC (of the read) to stratify.
     * Since the reads are expected to be seen over and over again, the GC is cached
     * <p>
     * Stratification happens into (integer) percents.
     */
    public static class GCContentStratifier extends RecordStratifier<Double> {

        // a cache to keep the GC of each read, since we will be visiting each read multiple times
        final Cache<SAMRecord, Double> gcCache = CacheBuilder.newBuilder().maximumSize(GC_CACHE_SIZE).build();

        @Override
        public Double stratify(final SAMRecord sam) {

            try {
                return gcCache.get(sam, () -> (double) Math.round(100 * SequenceUtil.calculateGc(sam.getReadBases())) / 100D);
            } catch (final ExecutionException e) {
                e.printStackTrace();
                throw new RuntimeException(e);
            }
        }

        @Override
        public String getSuffix() {
            return "gc";
        }
    }

    /**
     * Stratify bases according to the type of Homopolymer that they belong to (repeating element, final reference base and
     * whether the length is "long" or not). Read direction and only the preceding bases are taken into account.
     */

    public static class LongShortHomopolymerStratifier implements RecordAndOffsetStratifier<LongShortHomopolymer> {
        final int longHomopolymer;

        @Override
        public LongShortHomopolymer stratify(final RecordAndOffset recordAndOffset,
                                             final SAMLocusAndReference locusInfo) {

            final Integer hpLength = homoPolymerLengthStratifier.stratify(recordAndOffset, locusInfo);
            if (hpLength == null) {
                return null;
            }

            return hpLength < longHomopolymer ? LongShortHomopolymer.SHORT_HOMOPOLYMER : LongShortHomopolymer.LONG_HOMOPOLYMER;

        }

        LongShortHomopolymerStratifier(final int longHomopolymer) {
            this.longHomopolymer = longHomopolymer;
        }

        @Override
        public String getSuffix() {
            return "long_short_homopolymer";
        }
    }

    /**
     * Stratify by tags used during duplex and single index consensus calling.
     */
    public static class ConsensusStratifier extends RecordStratifier<Consensus> {
        // FIRST_STRAND_TAG refers to the number of copies of the first UMI, SECOND_STRAND_TAG is the number of the second
        // UMI.  BOTH_STRANDS_TAG is the sum of FIRST_STRAND_TAG and SECOND_STRAND_TAG in duplex consensus reads, and is the actual number
        // of copies of a read in single index consensus calling.

        /**
         * Tag that stores the number of duplicates of the first strand used in consensus calling.
         */
        static final String FIRST_STRAND_TAG = "aD";
        /**
         * Tag that stores the number of duplicates of the second strand used in consensus calling.
         */
        static final String SECOND_STRAND_TAG = "bD";
        /**
         * Tag that stores the total number of duplicates used in consensus calling.
         */
        static final String BOTH_STRANDS_TAG = "cD";

        @Override
        public Consensus stratify(final SAMRecord sam) {
            final int copiesOfFirstStrand, copiesOfSecondStrand, copiesOfBothStrands;

            // FIRST_STRAND_TAG and SECOND_STRAND_TAG should always be paired.  It makes no sense to have one and not the other.
            if (sam.hasAttribute(FIRST_STRAND_TAG) && sam.hasAttribute(SECOND_STRAND_TAG)) {
                copiesOfFirstStrand = sam.getIntegerAttribute(FIRST_STRAND_TAG);
                copiesOfSecondStrand = sam.getIntegerAttribute(SECOND_STRAND_TAG);
            } else {
                copiesOfFirstStrand = 0;
                copiesOfSecondStrand = 0;
            }

            // It's possible to have BOTH_STRANDS_TAG, but not FIRST_STRAND_TAG or SECOND_STRAND_TAG so we check here independently.
            if (sam.hasAttribute(BOTH_STRANDS_TAG)) {
                copiesOfBothStrands = sam.getIntegerAttribute(BOTH_STRANDS_TAG);
            } else {
                copiesOfBothStrands = 0;
            }

            if (copiesOfBothStrands == 1) {
                // Only one copy of one strand observed
                return Consensus.SIMPLEX_SINGLETON;
            } else if (copiesOfSecondStrand == 0 && copiesOfBothStrands > 1) {
                // Multiple copies of one strand observed, with no observations of the second strand
                return Consensus.SIMPLEX_CONSENSUS;
            } else if (copiesOfFirstStrand > 0 && copiesOfSecondStrand > 0 && (copiesOfFirstStrand == 1 || copiesOfSecondStrand == 1)) {
                // Both strands were observed at least once, but one strand was observed only once
                return Consensus.DUPLEX_SINGLETON;
            } else if (copiesOfFirstStrand > 1 && copiesOfSecondStrand > 1) {
                // Both strands were observed more than once
                return Consensus.DUPLEX_CONSENSUS;
            }

            // Either the read doesn't have consensus tags, or the tags don't make sense.
            return Consensus.UNKNOWN;
        }

        @Override
        public String getSuffix() {
            return "consensus";
        }
    }

    /**
     * Stratify by the number of Ns found in the read.
     * This is particularly useful for data that has been consensus-called (a process that can add 'N' bases when there is no consensus)
     */
    public static class NsInReadStratifier extends RecordStratifier<Integer> {
        private static String numberOfNsTag = "numberOfNs";

        @Override
        public Integer stratify(final SAMRecord sam) {
            int numberOfNsInRead = 0;

            // Check to see if we've already seen this record before.  If not, count the number
            // of Ns in the read, and store it in a transient attribute.  If we have seen it
            // before, simply retrieve the value and avoid doing the computation again.
            if (sam.getTransientAttribute(numberOfNsTag) != null) {
                numberOfNsInRead = (Integer) sam.getTransientAttribute(numberOfNsTag);
            } else {
                final byte[] bases = sam.getReadBases();
                for (final byte base : bases) {
                    if (SequenceUtil.isNoCall(base)) {
                        numberOfNsInRead++;
                    }
                }
                sam.setTransientAttribute(numberOfNsTag, numberOfNsInRead);
            }
            return numberOfNsInRead;
        }

        @Override
        public String getSuffix() {
            return "ns_in_read";
        }
    }

    /**
     * Stratifies into quintiles of read cycle.
     */
    public static class BinnedReadCycleStratifier implements RecordAndOffsetStratifier<CycleBin> {
        @Override
        public CycleBin stratify(final RecordAndOffset recordAndOffset,
                                 final SAMLocusAndReference locusInfo) {

            final int readCycle = stratifyCycle(recordAndOffset);
            final double relativePosition = (double) readCycle / recordAndOffset.getRecord().getReadLength();
            return CycleBin.valueOf(relativePosition);
        }

        @Override
        public String getSuffix() {
            return "binned_cycle";
        }
    }

    /**
     * Stratifies according to the overall mismatches (from {@link SAMTag#NM}) that the read has against the reference, NOT
     * including the current base.
     */
    public static class MismatchesInReadStratifier implements RecordAndOffsetStratifier<Integer> {
        @Override
        public Integer stratify(final RecordAndOffset recordAndOffset,
                                final SAMLocusAndReference locusInfo) {
            Integer numberMismatches = recordAndOffset.getRecord().getIntegerAttribute(SAMTag.NM.name());

            // Record may not contain an NM tag in which case we cannot stratify over it.
            if (numberMismatches == null) {
                return null;
            }

            // We are only interested in the number of errors on the read
            // not including the current base.
            if (recordAndOffset.getReadBase() != locusInfo.getReferenceBase()) {
                return numberMismatches - 1;
            }
            return numberMismatches;
        }

        @Override
        public String getSuffix() {
            return "mismatches_in_read";
        }
    }

    /**
     * Stratifies base into their read's tile which is parsed from the read-name.
     */
    public static class FlowCellTileStratifier extends RecordStratifier<Integer> {
        private static OpticalDuplicateFinder opticalDuplicateFinder = new OpticalDuplicateFinder();

        @Override
        public Integer stratify(final SAMRecord sam) {
            try {
                final PhysicalLocation location = new PhysicalLocationInt();
                opticalDuplicateFinder.addLocationInformation(sam.getReadName(), location);
                return (int) location.getTile();
            } catch (final IllegalArgumentException ignored) {
                return null;
            }
        }

        @Override
        public String getSuffix() {
            return "tile";
        }
    }

    /**
     * Stratifies according to the number of matching cigar operators (from CIGAR string) that the read has.
     */
    public static class CigarOperatorsInReadStratifier extends RecordStratifier<Integer> {

        private CigarOperator operator;

        public CigarOperatorsInReadStratifier(final CigarOperator op) {
            operator = op;
        }

        @Override
        public Integer stratify(final SAMRecord samRecord) {
            return stratifyCigarOperatorsInRead(samRecord, operator);
        }

        @Override
        public String getSuffix() {
            return "cigar_elements_" + operator.name() + "_in_read";
        }
    }

    /**
     * Stratifies according to the number of indel bases (from CIGAR string) that the read has.
     */
    public static class IndelsInReadStratifier extends RecordStratifier<Integer> {

        /**
         * Returns the number of bases associated with I and D CIGAR elements.
         * @param samRecord The read to investigate
         * @return The number of bases associated with I and D CIGAR elements, or null if the evaluation of either
         *         operation caused an error
         */
        @Override
        public Integer stratify(final SAMRecord samRecord) {
            final Integer insertedBasesInRead = stratifyCigarOperatorsInRead(samRecord, CigarOperator.I);
            final Integer deletedBasesInRead = stratifyCigarOperatorsInRead(samRecord, CigarOperator.D);
            if (insertedBasesInRead == null || deletedBasesInRead == null) {
                return null;
            }
            return insertedBasesInRead + deletedBasesInRead;
        }

        @Override
        public String getSuffix() {
            return "indels_in_read";
        }
    }

    /**
     * Stratifies according to the length of an insertion or deletion.
     */
    public static class IndelLengthStratifier implements RecordAndOffsetStratifier {

        @Override
        public Integer stratify(final RecordAndOffset recordAndOffset, final SAMLocusAndReference locusInfo) {
            return stratifyIndelLength(recordAndOffset, locusInfo);
        }

        @Override
        public String getSuffix() {
            return "indel_length";
        }
    }


    /* ************ Instances of stratifiers so that they do not need to be allocated every time ****************/

    /**
     * Stratifies bases into the current (uppercase) base as it was read from the sequencer (i.e. complemented if needed)
     */
    public static final RecordAndOffsetStratifier<Character> currentReadBaseStratifier = wrapStaticFunction(rao -> stratifyReadBase(rao, 0), "read_base");

    /**
     * Stratifies bases into the previous (uppercase) base as it was read from the sequencer (i.e. complemented if needed)
     */
    public static final RecordAndOffsetStratifier<Character> previousReadBaseStratifier = wrapStaticFunction(rao -> stratifyReadBase(rao, -1), "prev_base");
    /**
     * Stratifies bases into the following (uppercase) base as it was read from the sequencer (i.e. complemented if needed)
     */
    public static final RecordAndOffsetStratifier<Character> nextReadBaseStratifier = wrapStaticFunction(rao -> stratifyReadBase(rao, 1), "next_base");

    /**
     * Stratifies a base onto the reference base that it covers, possibly reverse complemented if the read
     * has been reversed by the aligner.
     */
    public static final RecordAndOffsetStratifier<Character> referenceBaseStratifier =
            wrapStaticFunction(ReadBaseStratification::stratifyReferenceBase, "ref_base");

    /**
     * Stratifies a base onto the reference base that it covers and the following base, possibly reverse complemented if the read
     * has been reversed by the aligner.
     */
    public static final PairStratifier<Character, Character> postDiNucleotideStratifier =
            PairStratifierFactory(referenceBaseStratifier, nextReadBaseStratifier, "post_dinuc");

    /**
     * Stratifies a base onto the reference base that it covers and the preceding base, possibly reverse complemented if the read
     * has been reversed by the aligner.
     */
    public static final PairStratifier<Character, Character> preDiNucleotideStratifier =
            PairStratifierFactory(previousReadBaseStratifier, referenceBaseStratifier, "pre_dinuc");

    /**
     * Stratifies a base onto the length of the homopolymer preceding it (in the read direction).
     * Read direction and only the preceding bases are taken into account (i.e. ignoring reference and current base).
     */
    public static final RecordAndOffsetStratifier<Integer> homoPolymerLengthStratifier =
            wrapStaticFunction(ReadBaseStratification::stratifyHomopolymerLength, "homopolymer_length");

    /**
     * Stratifies a base onto the make up (repeating base and following reference base) and length of the homopolymer from whence it came.
     * Read direction and only the preceding bases are taken into account.
     */
    public static final PairStratifier<Integer, Pair<Character, Character>> homopolymerStratifier =
            PairStratifierFactory(homoPolymerLengthStratifier, preDiNucleotideStratifier, "homopolymer_and_following_ref_base");

    /**
     * Stratifies a base onto the make up (repeating base and following reference base) and length-scale of the homopolymer from whence it came.
     * Read direction and only the preceding bases are taken into account.
     * <p>
     * This is done with a Lazy wrapper since it requires access to {@link #LONG_HOMOPOLYMER} which can be set with {@link #setLongHomopolymer(int)}.
     */
    public static final Lazy<PairStratifier<LongShortHomopolymer, Pair<Character, Character>>> binnedHomopolymerStratifier = new Lazy<>(() -> PairStratifierFactory(new LongShortHomopolymerStratifier(LONG_HOMOPOLYMER), preDiNucleotideStratifier,
            "binned_length_homopolymer_and_following_ref_base"));

    /**
     * Stratifies into the one base context of the base, which includes one read-base on each side and the reference base from where the current base is
     * the bases will be reverse-complemented so that the bases are in the original order they were read from the sequencer
     */
    public static final RecordAndOffsetStratifier<String> oneBasePaddedContextStratifier =
            wrapStaticFunction((rao, lar) -> stratifySurroundingContext(rao, lar, 1, 1), "one_base_padded_context");
    /**
     * Stratifies into the two base context of the base, which includes two read-bases on each side and the reference base from where the current base is
     * the bases will be reverse-complemented so that the bases are in the original order they were read from the sequencer
     */
    public static final RecordAndOffsetStratifier<String> twoBasePaddedContextStratifier =
            wrapStaticFunction((rao, lar) -> stratifySurroundingContext(rao, lar, 2, 2), "two_base_padded_context");

    /**
     * A constant stratifier which places all the reads into a single stratum.
     */
    public static final RecordStratifier<String> nonStratifier = wrapStaticReadFunction(sam -> "all", "all");

    /**
     * A stratifier that uses GC (of the read) to stratify.
     * Since the reads are expected to be seen over and over again, the GC is cached
     * <p>
     * Stratification happens into (integer) percents.
     */
    public static final GCContentStratifier gcContentStratifier = new GCContentStratifier();

    /**
     * Stratifies base into their read's tile which is parsed from the read-name.
     */
    public static final FlowCellTileStratifier flowCellTileStratifier = new FlowCellTileStratifier();

    /**
     * Stratifies to the readgroup id of the read.
     */
    public static final RecordStratifier<String> readgroupStratifier = wrapStaticReadFunction(ReadBaseStratification::stratifyReadGroup, "read_group");

    /**
     * Stratifies bases into their read's Ordinality (i.e. First or Second)
     */
    public static final RecordAndOffsetStratifier<ReadOrdinality> readOrdinalityStratifier = wrapStaticReadFunction(ReadOrdinality::of, "read_ordinality");

    /**
     * Stratifies bases into their read's Proper-pairedness
     */
    public static final RecordAndOffsetStratifier<ProperPaired> readPairednessStratifier = wrapStaticReadFunction(ProperPaired::of, "pair_proper");

    /**
     * Stratifies bases into their read's Direction (i.e. forward or reverse)
     */
    public static final RecordAndOffsetStratifier<ReadDirection> readDirectionStratifier = wrapStaticReadFunction(ReadDirection::of, "read_direction");

    /**
     * Stratifies bases into their read-pair's Orientation (i.e. F1R2, F2R1, F1F2 or R1R2)
     */
    public static final RecordAndOffsetStratifier<PairOrientation> readOrientationStratifier = wrapStaticReadFunction(PairOrientation::of, "pair_orientation");

    /**
     * Stratifies into quintiles of read cycle.
     */
    public static final BinnedReadCycleStratifier binnedReadCycleStratifier = new BinnedReadCycleStratifier();

    /**
     * Get the one-based cycle number of the base, taking the direction of the read into account
     */
    public static final RecordAndOffsetStratifier<Integer> baseCycleStratifier = wrapStaticFunction(ReadBaseStratification::stratifyCycle, "cycle");

    /**
     * Stratifies into the read-pairs estimated insert-length, as long as it isn't larger than 10x the length of the read
     */
    public static final RecordAndOffsetStratifier<Integer> insertLengthStratifier = wrapStaticReadFunction(ReadBaseStratification::stratifyInsertLength, "insert_length");

    /**
     * Stratifies into the number of soft-clipped bases that the read has in its alignment, or {@value NOT_ALIGNED_ERROR} if not aligned.
     */
    public static final RecordAndOffsetStratifier<Integer> softClipsLengthStratifier = wrapStaticReadFunction(ReadBaseStratification::stratifySoftClippedBases, "softclipped_bases");

    /**
     * Stratifies into the base-quality of the base under consideration
     */
    public static final RecordAndOffsetStratifier<Byte> baseQualityStratifier = wrapStaticFunction(ReadBaseStratification::stratifyBaseQuality, "base_quality");

    /**
     * Stratifies into the mapping-quality of the read under consideration
     */
    public static final RecordAndOffsetStratifier<Integer> mappingQualityStratifier = wrapStaticReadFunction(ReadBaseStratification::stratifyMappingQuality, "mapping_quality");

    /**
     * Stratifies according to the overall mismatches (from NM) that the read has against the reference, NOT
     * including the current base.
     */
    public static final MismatchesInReadStratifier mismatchesInReadStratifier = new MismatchesInReadStratifier();

    /**
     * Stratify by tags used during duplex and single index consensus calling.
     */
    public static final ConsensusStratifier consensusStratifier = new ConsensusStratifier();

    /**
     * Stratify by the number of Ns found in the read.
     * This is particularly useful for data that has been consensus-called (a process that can add 'N' bases when there is no consensus)
     */
    public static final NsInReadStratifier nsInReadStratifier = new NsInReadStratifier();

    /**
     * Stratify by Insertions in the read cigars.
     */
    public static final CigarOperatorsInReadStratifier insertionsInReadStratifier = new CigarOperatorsInReadStratifier(CigarOperator.I);

    /**
     * Stratify by Deletions in the read cigars.
     */
    public static final CigarOperatorsInReadStratifier deletionsInReadStratifier = new CigarOperatorsInReadStratifier(CigarOperator.D);

    /**
     * Stratify by Indels in the read cigars.
     */
    public static final IndelsInReadStratifier indelsInReadStratifier = new IndelsInReadStratifier();

    /**
     * Stratifies into the number of bases in an insertion
     */
    public static final IndelLengthStratifier indelLengthStratifier = new IndelLengthStratifier();

    /* *************** enums **************/

    /**
     * An Enum that is used to generate stratifiers from strings
     * <p>
     * To use this given a String 'str':
     * <p>
     * Stratifier.valueOf(str).makeStratifier()
     * <p>
     * This is used in {@link CollectSamErrorMetrics} to convert an input argument to a fully functional {@link BaseErrorAggregation} object.
     */
    enum Stratifier implements CommandLineParser.ClpEnum {
        ALL(() -> nonStratifier, "Puts all bases in the same stratum."),
        GC_CONTENT(() -> gcContentStratifier, "The GC-content of the read."),
        READ_ORDINALITY(() -> readOrdinalityStratifier, "The read ordinality (i.e. first or second)."),
        READ_BASE(() -> currentReadBaseStratifier, "the base in the original reading direction."),
        READ_DIRECTION(() -> readDirectionStratifier, "The alignment direction of the read (encoded as + or -)."),
        PAIR_ORIENTATION(() -> readOrientationStratifier, "The read-pair's orientation (encoded as '[FR]1[FR]2')."),
        PAIR_PROPERNESS(() -> readPairednessStratifier, "The properness of the read-pair's alignment. Looks for indications of chimerism."),
        REFERENCE_BASE(() -> referenceBaseStratifier, "The reference base in the read's direction."),
        PRE_DINUC(() -> preDiNucleotideStratifier, "The read base at the previous cycle, and the current reference base."),
        POST_DINUC(() -> postDiNucleotideStratifier, "The read base at the subsequent cycle, and the current reference base."),
        HOMOPOLYMER_LENGTH(() -> homoPolymerLengthStratifier, "The length of homopolymer the base is part of (only accounts for bases that were read prior to the current base)."),
        HOMOPOLYMER(() -> homopolymerStratifier, "The length of homopolymer, the base that the homopolymer is comprised of, and the reference base."),
        //using a lazy initializer to enable the value of LONG_HOMOPOLYMER to be used;
        BINNED_HOMOPOLYMER(binnedHomopolymerStratifier::get, "The scale of homopolymer (long or short), the base that the homopolymer is comprised of, and the reference base."),
        FLOWCELL_TILE(() -> flowCellTileStratifier, "The flowcell and tile where the base was read (taken from the read name)."),
        READ_GROUP(() -> readgroupStratifier, "The read-group id of the read."),
        CYCLE(() -> baseCycleStratifier, "The machine cycle during which the base was read."),
        BINNED_CYCLE(() -> binnedReadCycleStratifier, "The binned machine cycle. Similar to CYCLE, but binned into 5 evenly spaced ranges across the size of the read.  This stratifier may produce confusing results when used on datasets with variable sized reads."),
        SOFT_CLIPS(() -> softClipsLengthStratifier, "The number of softclipped bases the read has."),
        INSERT_LENGTH(() -> insertLengthStratifier, "The insert-size they came from (taken from the TLEN field.)"),
        BASE_QUALITY(() -> baseQualityStratifier, "The base quality."),
        MAPPING_QUALITY(() -> mappingQualityStratifier, "The read's mapping quality."),
        MISMATCHES_IN_READ(() -> mismatchesInReadStratifier, "The number of bases in the read that mismatch the reference, excluding the current base.  This stratifier requires the NM tag."),
        ONE_BASE_PADDED_CONTEXT(() -> oneBasePaddedContextStratifier, "The current reference base and a one base padded region from the read resulting in a 3-base context."),
        TWO_BASE_PADDED_CONTEXT(() -> twoBasePaddedContextStratifier, "The current reference base and a two base padded region from the read resulting in a 5-base context."),
        CONSENSUS(() -> consensusStratifier, "Whether or not duplicate reads were used to form a consensus read.  This stratifier makes use of the aD, bD, and cD tags for duplex consensus reads.  If the reads are single index consensus, only the cD tags are used."),
        NS_IN_READ(() -> nsInReadStratifier, "The number of Ns in the read."),
        INSERTIONS_IN_READ(() -> insertionsInReadStratifier, "The number of Insertions in the read cigar."),
        DELETIONS_IN_READ(() -> deletionsInReadStratifier, "The number of Deletions in the read cigar."),
        INDELS_IN_READ(() -> indelsInReadStratifier, "The number of INDELs in the read cigar."),
        INDEL_LENGTH(() -> indelLengthStratifier, "The number of bases in an indel");

        private final String docString;

        @Override
        public String getHelpDoc() {
            return docString + " Suffix is '" + stratifier.get().getSuffix() + "'.";
        }

        private final Supplier<RecordAndOffsetStratifier<?>> stratifier;

        Stratifier(final Supplier<RecordAndOffsetStratifier<?>> stratifier, final String docString) {
            this.stratifier = stratifier;
            this.docString = docString;
        }

        public RecordAndOffsetStratifier<?> makeStratifier() {
            return stratifier.get();
        }
    }

    /**
     * An enum to hold the ordinality of a read
     */
    public enum ReadOrdinality {
        FIRST,
        SECOND;

        public static ReadOrdinality of(final SAMRecord sam) {
            if (!sam.getReadPairedFlag()) {
                return null;
            }
            return sam.getFirstOfPairFlag() ? ReadOrdinality.FIRST : ReadOrdinality.SECOND;
        }
    }

    /**
     * An enum to hold information about the "properness" of a read pair
     */
    public enum ProperPaired {
        // Read has no supplementary alignments, is aligned to the same reference as its mate,
        // and has the "properly aligned" flag 0x2 set.
        PROPER,
        // Read has no supplementary alignments, is aligned to the same reference as its mate,
        // and has the "properly aligned" flag 0x2 UN-set.
        IMPROPER,
        // Read is softclipped, is aligned to the same reference as its mate, and
        // has at least one supplementary alignment
        CHIMERIC,
        // Read is aligned to a different reference than its mate
        DISCORDANT,
        // Read or Mate are not aligned, or read has no mate and cannot be declared to be CHIMERIC.
        UNKNOWN;

        public static ProperPaired of(final SAMRecord sam) {

            if (sam.getReadPairedFlag() &&
                    !sam.getMateUnmappedFlag() && !sam.getReadUnmappedFlag() &&
                    !sam.getMateReferenceIndex().equals(sam.getReferenceIndex())) {
                return DISCORDANT;
            }

            if (!sam.getReadUnmappedFlag() && sam.getCigar().isClipped() && sam.hasAttribute(SAMTag.SA.toString())) {
                return CHIMERIC;
            }

            if (sam.getReadUnmappedFlag() || sam.getReadPairedFlag() && sam.getMateUnmappedFlag()) {
                return UNKNOWN;
            }

            if (!sam.getProperPairFlag()) {
                return IMPROPER;
            }

            return PROPER;
        }
    }

    /**
     * An enum designed to hold a binned version of any probability-like number (between 0 and 1)
     * in quintiles
     */
    public enum CycleBin {
        // CycleBins are groupings of inputs by their binned quintile (20%)
        // Each bin represents 20% of the distribution.
        QUINTILE_1(0, 0.2),
        QUINTILE_2(0.2, 0.4),
        QUINTILE_3(0.4, 0.6),
        QUINTILE_4(0.6, 0.8),
        QUINTILE_5(0.8, 1.0);

        final double lower;
        final double upper;

        CycleBin(double lower, double upper) {
            this.lower = lower;
            this.upper = upper;
        }

        static CycleBin valueOf(double value) {
            return Stream.of(CycleBin.values())
                    .filter(e -> value >= e.lower && value <= e.upper)
                    .findFirst()
                    .orElseThrow(() -> new IllegalArgumentException(String.format("Value for CycleBin must be between 0 and 1 (inclusive), found: %g", value)));
        }
    }

    /**
     * An enum for holding the direction for a read (positive strand or negative strand
     */
    public enum ReadDirection {
        POSITIVE("+"),
        NEGATIVE("-");

        private final String outputString;

        ReadDirection(final String output) {
            outputString = output;
        }

        @Override
        public String toString() {
            return outputString;
        }

        public static ReadDirection of(final SAMRecord sam) {
            return sam.getReadNegativeStrandFlag() ? ReadDirection.NEGATIVE : ReadDirection.POSITIVE;
        }
    }

    /**
     * An enum for holding a reads read-pair's Orientation (i.e. F1R2, F2R1, or TANDEM) indicating
     * the strand (positive or negative) that each of the two mated reads are aligned to. In connection with
     * READ_BASE and similar stratifiers this can be used to observe oxoG-type sequencing artifacts.
     * <p>
     * Note that this is not related to FR/RF/TANDEM classification as per {@link htsjdk.samtools.SamPairUtil.PairOrientation PairOrientation}.
     */
    public enum PairOrientation {
        F1R2,
        F2R1,
        F1F2,
        R1R2;

        private static PairOrientation ofExplicit(final boolean firstPositive, final boolean secondPositive) {
            if (firstPositive == secondPositive) {
                return firstPositive ? PairOrientation.R1R2 : PairOrientation.F1F2;
            } else {
                return firstPositive ? PairOrientation.F1R2 : PairOrientation.F2R1;
            }
        }

        public static PairOrientation of(final SAMRecord sam) {
            final ReadOrdinality ordinality = ReadOrdinality.of(sam);
            final ReadDirection direction = ReadDirection.of(sam);

            // Make sure that the read is paired
            if (!sam.getReadPairedFlag()) {
                return null;
            }

            // if read is unmapped, read isn't mapped, or mate is unmapped, return null
            if (direction == null ||
                    sam.getReadUnmappedFlag() ||
                    sam.getMateUnmappedFlag()) {
                return null;
            }

            final boolean matePositiveStrand = !sam.getMateNegativeStrandFlag();

            if (ordinality == ReadOrdinality.FIRST) {
                return PairOrientation.ofExplicit(direction == ReadDirection.POSITIVE, matePositiveStrand);
            } else {
                return PairOrientation.ofExplicit(matePositiveStrand, direction == ReadDirection.POSITIVE);
            }
        }
    }

    public enum LongShortHomopolymer {
        SHORT_HOMOPOLYMER,
        LONG_HOMOPOLYMER
    }

    /**
     * Types of consensus reads as determined by the number of duplicates used from
     * first and second strands.
     * <li>{@link #SIMPLEX_SINGLETON}</li>
     * <li>{@link #SIMPLEX_CONSENSUS}</li>
     * <li>{@link #DUPLEX_SINGLETON}</li>
     * <li>{@link #DUPLEX_CONSENSUS}</li>
     * <li>{@link #UNKNOWN}</li>
     */
    public enum Consensus {
        /**
         * Read which only one observation is made.  In a sense, these are not consensus reads.
         */
        SIMPLEX_SINGLETON,
        /**
         * Read which had multiple observations, but all from the same strand (or the strand could not be determined).
         */
        SIMPLEX_CONSENSUS,
        /**
         * Read which has two observations, one from each direction.
         */
        DUPLEX_SINGLETON,
        /**
         * Consensus read which has two or more observations in both directions.
         */
        DUPLEX_CONSENSUS,
        /**
         * Read whose consensus status cannot be determined.
         */
        UNKNOWN
    }

    /**
     * Returns the number of bases associated with a specific cigar operator within a read
     * @param samRecord The read to investigate
     * @param operator The operator that the counted bases should be associated with
     */
    private static Integer stratifyCigarOperatorsInRead(final SAMRecord samRecord, final CigarOperator operator) {
        try {
            return samRecord.getCigar().getCigarElements().stream()
                    .filter(ce -> ce.getOperator().equals(operator))
                    .mapToInt(CigarElement::getLength)
                    .sum();
        } catch (final Exception ex) {
            return null;
        }
    }

    /**
     * A function for stratifying a recordAndOffset and a SAMLocus into a context of bases surrounding the base in question
     * Read direction is taken into consideration (all bases are reverse-complemented if needed so that the result is in the
     * same order as was originally read) and the base under consideration is replaced with the reference base at that location
     * All other bases are taken from the read, not the reference.
     */
    private static String stratifySurroundingContext(final RecordAndOffset recordAndOffset, final SAMLocusAndReference locusInfo, final int basesBefore, final int basesAfter) {
        final StringBuilder stringBuilder = new StringBuilder(basesAfter + basesBefore + 1);

        for (int offset = -basesBefore; offset <= basesAfter; offset++) {
            if (offset == 0) {
                // Use the reference base as the center of the context
                stringBuilder.append(stratifyReferenceBase(recordAndOffset, locusInfo));
            } else {
                // Use the surrounding read bases as the surrounding context
                final Character surroundingBase = stratifyReadBase(recordAndOffset, offset);
                if (surroundingBase == null) {
                    return null;
                } else {
                    stringBuilder.append(surroundingBase);
                }
            }
        }

        return stringBuilder.toString();
    }

    /**
     * Get the one-based cycle number of the base, taking the direction of the read into account
     */
    private static int stratifyCycle(final RecordAndOffset recordAndOffset) {
        final SAMRecord rec = recordAndOffset.getRecord();
        final int offset = recordAndOffset.getOffset();
        // Get either the offset into the array or the distance from the end depending on whether the read is
        // on the negative strand.
        int retval = rec.getReadNegativeStrandFlag() ? (rec.getReadLength() - offset - 1) : offset;
        // add 1 to move to a one-based system
        retval += 1;

        return retval;
    }

    private static Integer stratifyHomopolymerLength(final RecordAndOffset recordAndOffset, final SAMLocusAndReference locusInfo) {

        final ReadDirection direction = ReadDirection.of(recordAndOffset.getRecord());
        final byte readBases[] = recordAndOffset.getRecord().getReadBases();
        if (SequenceUtil.isNoCall(locusInfo.getReferenceBase())) {
            return null;
        }
        int runLengthOffset = recordAndOffset.getOffset();

        if (runLengthOffset < 0 || runLengthOffset >= recordAndOffset.getRecord().getReadLength()) {
            return null;
        }

        if (direction == ReadDirection.POSITIVE) {
            while (--runLengthOffset >= 0) {
                if (readBases[runLengthOffset] != readBases[recordAndOffset.getOffset() - 1]) {
                    break;
                }
            }
            return recordAndOffset.getOffset() - runLengthOffset - 1;
        } else {
            while (++runLengthOffset < recordAndOffset.getRecord().getReadLength()) {
                if (readBases[runLengthOffset] != readBases[recordAndOffset.getOffset() + 1]) {
                    break;
                }
            }
            return runLengthOffset - recordAndOffset.getOffset() - 1;
        }
    }

    private static Character stratifyReadBase(final RecordAndOffset recordAndOffset, int offset) {
        final ReadDirection direction = ReadDirection.of(recordAndOffset.getRecord());

        final int requestedOffset = recordAndOffset.getOffset() + offset * (direction == ReadDirection.POSITIVE ? 1 : -1);

        if (requestedOffset < 0 || requestedOffset >= recordAndOffset.getRecord().getReadLength()) {
            return null;
        } else {
            return stratifySequenceBase(recordAndOffset.getRecord().getReadBases()[requestedOffset], direction == ReadDirection.NEGATIVE);
        }
    }

    private static Character stratifyReferenceBase(final RecordAndOffset recordAndOffset,
                                                   final SAMLocusAndReference locusInfo) {
        final ReadDirection direction = ReadDirection.of(recordAndOffset.getRecord());

        if (SequenceUtil.isNoCall(locusInfo.getReferenceBase())) {
            return null;
        }

        return stratifySequenceBase(locusInfo.getReferenceBase(), direction == ReadDirection.NEGATIVE);
    }

    private static char stratifySequenceBase(final byte input, final Boolean getComplement) {
        return (char) SequenceUtil.upperCase(getComplement ? SequenceUtil.complement(input) : input);
    }

    private static Integer stratifyInsertLength(final SAMRecord sam) {
        return Math.min(
                sam.getReadLength() * 10,
                Math.abs(sam.getInferredInsertSize()));
    }

    public static final int NOT_ALIGNED_ERROR = -1;

    private static Integer stratifySoftClippedBases(final SAMRecord sam) {
        final Cigar cigar = sam.getCigar();
        if (cigar == null) {
            return NOT_ALIGNED_ERROR;
        }
        return cigar.getCigarElements().stream()
                .filter(e -> e.getOperator() == CigarOperator.S)
                .mapToInt(CigarElement::getLength).sum();
    }

    private static Byte stratifyBaseQuality(final RecordAndOffset recordAndOffset) {
        return recordAndOffset.getBaseQuality();
    }

    private static int stratifyMappingQuality(final SAMRecord sam) {
        return sam.getMappingQuality();
    }

    private static String stratifyReadGroup(final SAMRecord sam) {
        return sam.getReadGroup().getReadGroupId();
    }

    private static Integer stratifyIndelLength(final RecordAndOffset recordAndOffset, final SAMLocusAndReference locusInfo) {
        // If the base is not an indel, stratify it as a length of 0
        if (recordAndOffset.getAlignmentType() != AbstractRecordAndOffset.AlignmentType.Insertion &&
                recordAndOffset.getAlignmentType() != AbstractRecordAndOffset.AlignmentType.Deletion) {
            return 0;
        }

        final CigarElement cigarElement = getIndelElement(recordAndOffset);
        if (cigarElement == null) {
            // No CIGAR operation for given position.
            return null;
        }
        if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Insertion &&
                cigarElement.getOperator() != CigarOperator.I) {
            throw new IllegalStateException("Wrong CIGAR operator for the given position.");
        }
        if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Deletion &&
                cigarElement.getOperator() != CigarOperator.D) {
            throw new IllegalStateException("Wrong CIGAR operator for the given position.");
        }
        return cigarElement.getLength();
    }

    public static CigarElement getIndelElement(final RecordAndOffset recordAndOffset) {
        final SAMRecord record = recordAndOffset.getRecord();
        final int offset = recordAndOffset.getOffset();

        if (recordAndOffset.getAlignmentType() != AbstractRecordAndOffset.AlignmentType.Insertion &&
                recordAndOffset.getAlignmentType() != AbstractRecordAndOffset.AlignmentType.Deletion) {
            log.warn("This method is not supported for matching bases.");
            // But could be included by handling matches the same way as insertions?
            return null;
        }

        if (record == null) {
            throw new IllegalArgumentException("record must not be null.");
        }

        // -1 is still a valid input for a deletion (returns the first cigar element). Everything below that is an error.
        if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Insertion && offset < 0) {
            throw new IllegalArgumentException("offset must greater than zero for an insertion.");
        }

        if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Deletion && offset < -1) {
            throw new IllegalArgumentException("offset must greater than -1 for a deletion.");
        }

        final Cigar cigar = record.getCigar();
        if (cigar.isEmpty())
            return null;

        if (offset == -1) {
            return cigar.getCigarElement(0);
        }

        int readPosition = 0;
        for (CigarElement cigarElement : cigar) {
            if (readPosition > offset + 1) {
                // We somehow went past the desired location
                return null;
            }

            if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Insertion) {
                // If it doesn't consume bases, skip to the next
                if (cigarElement.getOperator().consumesReadBases() && readPosition == offset) {
                    return cigarElement;
                }
            }
            else if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Deletion) {
                if (readPosition == offset + 1) {
                    return cigarElement;
                }
            }

            readPosition += cigarElement.getOperator().consumesReadBases() ? cigarElement.getLength() : 0;
        }
        return null;
    }
}