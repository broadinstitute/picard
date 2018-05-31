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
package picard.sam.BamErrorMetric;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.PhysicalLocationInt;

import java.util.concurrent.ExecutionException;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Stream;

/**
 * Classes, methods, and enums that deal with the stratification of read bases and reference information.
 *
 * @author Yossi Farjoun
 */
public class ReadBaseStratification {
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
     *
     * @param <T>
     */
    public interface RecordAndOffsetStratifier<T extends Comparable<T>> {
        // The method that stratifies a base in a read (provided by RecordAndOffset) into a type T.
        T stratify(final SamLocusIterator.RecordAndOffset recordAndOffset, final SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo);

        // The string suffix that will be used to generate the extension of the metric file.
        String getSuffix();
    }

    /**
     * A simpler stratifier for cases when only the record suffices
     */
    abstract static class RecordStratifier<T extends Comparable<T>> implements RecordAndOffsetStratifier<T> {
        @Override
        public T stratify(SamLocusIterator.RecordAndOffset recordAndOffset, SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo) {
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
    private static <T extends Comparable<T>> RecordAndOffsetStratifier<T> wrapStaticFunction(BiFunction<SamLocusIterator.RecordAndOffset, SAMLocusAndReferenceIterator.SAMLocusAndReference, T> staticStratify, String suffix) {
        return new RecordAndOffsetStratifier<T>() {
            @Override
            public T stratify(SamLocusIterator.RecordAndOffset recordAndOffset, SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo) {
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
    private static <T extends Comparable<T>> RecordAndOffsetStratifier<T> wrapStaticFunction(Function<SamLocusIterator.RecordAndOffset, T> staticStratify, String suffix) {
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
     * The suffix for this stratifiers is automatically generated from the suffixes of the two provided
     * stratifiers.
     */
    public static class PairStratifier<T extends Comparable<T>, R extends Comparable<R>> implements RecordAndOffsetStratifier<Pair<T, R>> {
        public PairStratifier(final RecordAndOffsetStratifier<T> a, final RecordAndOffsetStratifier<R> b) {
            this.a = a;
            this.b = b;
        }

        final RecordAndOffsetStratifier<T> a;
        final RecordAndOffsetStratifier<R> b;

        @Override
        public Pair<T, R> stratify(final SamLocusIterator.RecordAndOffset recordAndOffset, final SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo) {

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
     * A factory for generating "pair" stratifier instances from two stratifiers and a string
     *
     * @param suffix the suffix to use
     * @param <T>    the type into which the stratification happens
     * @return an instance of a stratification class that will stratify accordingly
     */
    public static <T extends Comparable<T>, S extends Comparable<S>> PairStratifier<T, S> PairStratifierFactory(
            final RecordAndOffsetStratifier<T> stratifier1,
            final RecordAndOffsetStratifier<S> stratifier2,
            final String suffix) {
        return new PairStratifier<T, S>(stratifier1, stratifier2) {
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
        public LongShortHomopolymer stratify(final SamLocusIterator.RecordAndOffset recordAndOffset,
                                             final SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo) {

            final Integer hpLength = homoPolymerLengthStratifier.stratify(recordAndOffset, locusInfo);
            if (hpLength == null) return null;

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
        public CycleBin stratify(final SamLocusIterator.RecordAndOffset recordAndOffset,
                                 final SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo) {

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
     * Stratifies according to the overall mismatches (from NM) that the read has against the reference, NOT
     * including the current base.
     */
    public static class MismatchesInReadStratifier implements RecordAndOffsetStratifier<Integer> {
        @Override
        public Integer stratify(final SamLocusIterator.RecordAndOffset recordAndOffset,
                                final SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo) {
            Integer numberMismatches = recordAndOffset.getRecord().getIntegerAttribute(SAMTag.NM.name());

            // Record may not contain an NM tag in which case we cannot stratify over it.
            if (numberMismatches == null) {
                return null;
            }

            // We are only interested in the number of errors on the read
            // not including the current base.
            if (recordAndOffset.getReadBase() != locusInfo.referenceBase) {
                return numberMismatches - 1;
            } else {
                return numberMismatches;
            }
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
            final int tile;
            try {
                PhysicalLocation location = new PhysicalLocationInt();
                opticalDuplicateFinder.addLocationInformation(sam.getReadName(), location);
                tile = location.getTile();
                return tile;
            } catch (final IllegalArgumentException ignored) {
                return null;
            }
        }

        @Override
        public String getSuffix() {
            return "tile";
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
    public static final Lazy<PairStratifier<LongShortHomopolymer, Pair<Character, Character>>> binnedHomopolymerStratifier = new Lazy<>(() -> PairStratifierFactory( new LongShortHomopolymerStratifier(LONG_HOMOPOLYMER), preDiNucleotideStratifier,
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
     * Get the cycle number of the base, taking the direction of the read into account
     */
    public static final RecordAndOffsetStratifier<Integer> baseCycleStratifier = wrapStaticFunction(ReadBaseStratification::stratifyCycle, "cycle");

    /**
     * Stratifies into the read-pairs estimated insert-length, as long as it isn't larger than 10x the length of the read
     */
    public static final RecordAndOffsetStratifier<Integer> insertLengthStratifier = wrapStaticReadFunction(ReadBaseStratification::stratifyInsertLength, "insert_length");

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

    /* *************** enums **************/

    /**
     * An Enum that is used to generate stratifiers from strings
     * <p>
     * To use this given a String 'str':
     * <p>
     * Stratifiers.valueOf(str).makeStratifier()
     * <p>
     * This is used in {@link CollectBamErrorMetrics} to convert an input argument to a fully functional {@link BaseErrorAggregation} object.
     */
    enum Stratifiers implements CommandLineParser.ClpEnum {
        ALL(() -> nonStratifier, "Puts all bases in the same stratum."),
        GC_CONTENT(() -> gcContentStratifier, "Stratifies bases according to the gc content of their read."),
        READ_ORDINALITY(() -> readOrdinalityStratifier, "Stratifies bases according to their read ordinality (i.e. first or second)."),
        READ_BASE(() -> currentReadBaseStratifier, "Stratifies bases by the read in the original reading direction."),
        READ_DIRECTION(() -> readDirectionStratifier, "Stratifies bases to +/- based on the alignment direction of the read."),
        PAIR_ORIENTATION(() -> readOrientationStratifier, "Stratifies bases into F1R2/F2R1 according to their reads orientation and ordinality. Assumes reads are \"innies\"."),
        REFERENCE_BASE(() -> referenceBaseStratifier, "Stratifies bases by the read-directed reference base."),
        PRE_DINUC(() -> preDiNucleotideStratifier, "Stratifies bases by the read base at the previous cycle, and the current reference base."),
        POST_DINUC(() -> postDiNucleotideStratifier, "Stratifies bases by the read base at the previous cycle, and the current reference base."),
        HOMOPOLYMER_LENGTH(() -> homoPolymerLengthStratifier, "Stratifies bases based on the length of homopolymer they are part of (only accounts for bases that were read prior to the current base)."),
        HOMOPOLYMER(() -> homopolymerStratifier, "Stratifies bases based on the length of homopolymer, the base that the homopolymer is comprised of, and the reference base."),
        //using a lazy initializer to enable the value of LONG_HOMOPOLYMER to be used;
        BINNED_HOMOPOLYMER(binnedHomopolymerStratifier::get, "Stratifies bases based on the scale of homopolymer (long or short), the base that the homopolymer is comprised of, and the reference base."),
        FLOWCELL_TILE(() -> flowCellTileStratifier, "Stratifies according to the flowcell-tile where the base was read (taken from the read name)."),
        READ_GROUP(() -> readgroupStratifier, "Stratifies bases according to their read-group id."),
        CYCLE(() -> baseCycleStratifier, "Stratifies bases to the machine cycle during which they were read."),
        BINNED_CYCLE(() -> binnedReadCycleStratifier, "Stratifies bases according to the binned machine cycle in the read similar to CYCLE, but binned into 5 evenly spaced ranges across the size of the read.  This stratifier may produce confusing results when used on datasets with variable sized reads."),
        INSERT_LENGTH(() -> insertLengthStratifier, "Stratifies bases according to the insert-size they came from (taken from the TLEN field.)"),
        BASE_QUALITY(() -> baseQualityStratifier, "Stratifies bases according to their base quality."),
        MAPPING_QUALITY(() -> mappingQualityStratifier, "Stratifies bases according to their read's mapping quality."),
        MISMATCHES_IN_READ(() -> mismatchesInReadStratifier, "Stratifies bases according to the number of bases in the read that mismatch the reference excluding the current base.  This stratifier requires the NM tag."),
        ONE_BASE_PADDED_CONTEXT(() -> oneBasePaddedContextStratifier, "Stratifies bases according the current reference base and a one base padded region from the read resulting in a 3-base context."),
        TWO_BASE_PADDED_CONTEXT(() -> twoBasePaddedContextStratifier, "Stratifies bases according the current reference base and a two base padded region from the read resulting in a 5-base context."),
        CONSENSUS(() -> consensusStratifier, "Stratifies bases according to whether or not duplicate reads were used to form a consensus read.  This stratifier makes use of the aD, bD, and cD tags for duplex consensus reads.  If the reads are single index consensus, only the cD tags are used."),
        NS_IN_READ(() -> nsInReadStratifier, "Stratifies bases according to the number of Ns in the read.");

        private final String docString;

        @Override
        public String getHelpDoc() {
            return docString + " Suffix is " + stratifier.get().getSuffix();
        }

        private final Supplier<RecordAndOffsetStratifier<?>> stratifier;

        Stratifiers(final Supplier<RecordAndOffsetStratifier<?>> stratifier, final String docString) {
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
     * An enum designed to hold a binned version of
     * any probability-like number (between 0 and 1)
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
     * An enum for holding a reads read-pair's Orientation (i.e. F1R2, F2R1, F1F2 or R1R2)
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

            // if read is unmapped, read isn't mapped, or mate is unmapped, return null
            if (direction == null ||
                    sam.getReadUnmappedFlag() ||
                    sam.getMateUnmappedFlag()) return null;

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
     * A function for stratifying a recordAndOffset and a SAMLocus into a context of bases surrounding the base in question
     * Read direction is taken into consideration (all bases are reverse-complemented if needed so that the result is in the
     * same order as was originally read) and the base under consideration is replaced with the reference base at that location
     * All other bases are taken from the read, not the reference.
     */

    private static String stratifySurroundingContext(final SamLocusIterator.RecordAndOffset recordAndOffset, final SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo, final int basesBefore, final int basesAfter) {
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
     * Get the cycle number of the base, taking the direction of the read into account
     */
    private static int stratifyCycle(final SamLocusIterator.RecordAndOffset recordAndOffset) {
        final SAMRecord rec = recordAndOffset.getRecord();
        final int offset = recordAndOffset.getOffset();
        // get either the offset into the array or the distance from the end depending on whether the read is
        // on the negative strand.
        int retval = rec.getReadNegativeStrandFlag() ? (rec.getReadLength() - offset - 1) : offset;
        // add 1 to move to a 1-based system
        retval += 1;

        return retval;
    }

    private static Integer stratifyHomopolymerLength(final SamLocusIterator.RecordAndOffset recordAndOffset, final SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo) {

        final ReadDirection direction = ReadDirection.of(recordAndOffset.getRecord());
        final byte readBases[] = recordAndOffset.getRecord().getReadBases();
        if (SequenceUtil.isNoCall(locusInfo.referenceBase)) {
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

    private static Character stratifyReadBase(final SamLocusIterator.RecordAndOffset recordAndOffset, int offset) {
        final ReadDirection direction = ReadDirection.of(recordAndOffset.getRecord());

        final int requestedOffset = recordAndOffset.getOffset() + offset * (direction == ReadDirection.POSITIVE ? 1 : -1);

        if (requestedOffset < 0 || requestedOffset >= recordAndOffset.getRecord().getReadLength()) {
            return null;
        } else {
            return stratifySequenceBase(recordAndOffset.getRecord().getReadBases()[requestedOffset], direction == ReadDirection.NEGATIVE);
        }
    }

    private static Character stratifyReferenceBase(final SamLocusIterator.RecordAndOffset recordAndOffset,
                                                   final SAMLocusAndReferenceIterator.SAMLocusAndReference locusInfo) {
        final ReadDirection direction = ReadDirection.of(recordAndOffset.getRecord());

        if (SequenceUtil.isNoCall(locusInfo.referenceBase)) {
            return null;
        }

        return stratifySequenceBase(locusInfo.referenceBase, direction == ReadDirection.NEGATIVE);
    }

    private static char stratifySequenceBase(final byte input, final Boolean getComplement) {
        return (char) SequenceUtil.upperCase(getComplement ? SequenceUtil.complement(input) : input);
    }

    private static Integer stratifyInsertLength(final SAMRecord sam) {
        return Math.min(
                sam.getReadLength() * 10,
                Math.abs(sam.getInferredInsertSize()));
    }

    private static Byte stratifyBaseQuality(final SamLocusIterator.RecordAndOffset recordAndOffset) {
        return recordAndOffset.getBaseQuality();
    }

    private static int stratifyMappingQuality(final SAMRecord sam) {
        return sam.getMappingQuality();
    }

    private static String stratifyReadGroup(final SAMRecord sam) {
        return sam.getReadGroup().getReadGroupId();
    }
}
