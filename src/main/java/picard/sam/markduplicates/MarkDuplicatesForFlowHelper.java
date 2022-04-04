package picard.sam.markduplicates;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SortingLongCollection;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;
import picard.sam.markduplicates.util.RepresentativeReadIndexerCodec;
import picard.sam.util.RepresentativeReadIndexer;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * MarkDuplicates calculation helper class for flow based mode
 *
 * The class extends the behavior of MarkDuplicates which contains the complete
 * code for the non-flow based mode. When in flow mode, additional parameters
 * may control the establishment of read ends (start/end) for the purpose of
 * determining duplication status. Additionally, the logic used to gather reads into
 * (duplicate) buckets (chunks) is enhanced with an optional mechanism of read end
 * uncertainty threshold. When active, reads are considered to belong to the same chunk if
 * for each read in the chunk there exists at least one other read with the uncertainty
 * distance on the read end.
 */
public class MarkDuplicatesForFlowHelper implements MarkDuplicatesHelper {

    private final Log log = Log.getInstance(MarkDuplicatesForFlowHelper.class);

    private static final int END_INSIGNIFICANT_VALUE = 0;
    private static final String ATTR_DUPLICATE_SCORE = "ForFlowDuplicateScore";

    // constants for clippingTagContainsAny
    public static final String CLIPPING_TAG_NAME = "tm";
    public static final char[]  CLIPPING_TAG_CONTAINS_A = {'A'};
    public static final char[]  CLIPPING_TAG_CONTAINS_AQ = {'A', 'Q'};
    public static final char[]  CLIPPING_TAG_CONTAINS_QZ = {'Q', 'Z'};

    // instance of hosting MarkDuplicates
    private final MarkDuplicates  md;

    public MarkDuplicatesForFlowHelper(final MarkDuplicates md) {
        this.md = md;

        validateFlowParameteres();
    }

    private void validateFlowParameteres() {

        if ( md.flowBasedArguments.UNPAIRED_END_UNCERTAINTY != 0 && !md.flowBasedArguments.USE_END_IN_UNPAIRED_READS ) {
            throw new IllegalArgumentException("invalid parameter combination. UNPAIRED_END_UNCERTAINTY can not be specified when USE_END_IN_UNPAIRED_READS not specified");
        }
    }

    /**
     * This method is identical in function to generateDuplicateIndexes except that it accomodates for
     * the possible significance of the end side of the reads (w/ or w/o uncertainty). This is only
     * applicable for flow mode invocation.
     */
    public void generateDuplicateIndexes(final boolean useBarcodes, final boolean indexOpticalDuplicates) {
        final int entryOverhead;
        if (md.TAG_DUPLICATE_SET_MEMBERS) {
            // Memory requirements for RepresentativeReadIndexer:
            // three int entries + overhead: (3 * 4) + 4 = 16 bytes
            entryOverhead = 16;
        } else {
            entryOverhead = SortingLongCollection.SIZEOF;
        }
        // Keep this number from getting too large even if there is a huge heap.
        int maxInMemory = (int) Math.min((Runtime.getRuntime().maxMemory() * 0.25) / entryOverhead, (double) (Integer.MAX_VALUE - 5));
        // If we're also tracking optical duplicates, reduce maxInMemory, since we'll need two sorting collections
        if (indexOpticalDuplicates) {
            maxInMemory /= ((entryOverhead + SortingLongCollection.SIZEOF) / entryOverhead);
            md.opticalDuplicateIndexes = new SortingLongCollection(maxInMemory, md.TMP_DIR.toArray(new File[md.TMP_DIR.size()]));
        }
        log.info("Will retain up to " + maxInMemory + " duplicate indices before spilling to disk.");
        md.duplicateIndexes = new SortingLongCollection(maxInMemory, md.TMP_DIR.toArray(new File[md.TMP_DIR.size()]));
        if (md.TAG_DUPLICATE_SET_MEMBERS) {
            final RepresentativeReadIndexerCodec representativeIndexCodec = new RepresentativeReadIndexerCodec();
            md.representativeReadIndicesForDuplicates = SortingCollection.newInstance(RepresentativeReadIndexer.class,
                    representativeIndexCodec,
                    Comparator.comparing(read -> read.readIndexInFile),
                    maxInMemory,
                    md.TMP_DIR);
        }

        // this code does support pairs at this time
        if ( md.pairSort.iterator().hasNext() ) {
            throw new IllegalArgumentException("flow based code does not support paired reads");
        }
        md.pairSort.cleanup();
        md.pairSort = null;

        /**
         *  Now deal with the fragments
         *
         *  The end processing semantics depends on the following factors:
         *  1. Whether the end is marked as significant (as specified by USE_END_IN_UNPAIRED_READS)
         *  2. Whether end certainty is specified (by UNPAIRED_END_UNCERTAINTY)
         *
         *  - If ends are insignificant, they are ignored
         *  - If ends are significant and uncertainty is set to 0 - they must be equal for fragments to be considered same
         *  - Otherwise, fragments are accumulated (into the same bucket) as long as they are with the
         *  specified uncertainty from at least one existing fragment. Note that using this strategy the effective
         *  range of end locations associated with fragments in a bucket may grow, but only in 'uncertainty' steps.
         */
        log.info("Traversing fragment information and detecting duplicates.");
        ReadEndsForMarkDuplicates firstOfNextChunk = null;
        int nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
        int nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;
        final List<ReadEndsForMarkDuplicates> nextChunk = new ArrayList<>(200);
        boolean containsPairs = false;
        boolean containsFrags = false;

        for (final ReadEndsForMarkDuplicates next : md.fragSort) {
            if (firstOfNextChunk != null && areComparableForDuplicatesWithEndSignificance(firstOfNextChunk, next, useBarcodes,
                    nextChunkRead1Coordinate2Min, nextChunkRead1Coordinate2Max)) {
                nextChunk.add(next);
                containsPairs = containsPairs || next.isPaired();
                containsFrags = containsFrags || !next.isPaired();
                if ( next.read2Coordinate != END_INSIGNIFICANT_VALUE) {
                    nextChunkRead1Coordinate2Min = Math.min(nextChunkRead1Coordinate2Min, next.read2Coordinate);
                    nextChunkRead1Coordinate2Max = Math.max(nextChunkRead1Coordinate2Max, next.read2Coordinate);

                    if ( firstOfNextChunk.read2Coordinate == END_INSIGNIFICANT_VALUE)
                        firstOfNextChunk = next;
                }
            } else {
                if (nextChunk.size() > 1 && containsFrags) {
                    md.markDuplicateFragments(nextChunk, containsPairs);
                }
                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                if ( next.read2Coordinate != END_INSIGNIFICANT_VALUE)
                    nextChunkRead1Coordinate2Min = nextChunkRead1Coordinate2Max = next.read2Coordinate;
                else {
                    nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
                    nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;
                }
                containsPairs = next.isPaired();
                containsFrags = !next.isPaired();
            }
        }
        md.markDuplicateFragments(nextChunk, containsPairs);
        md.fragSort.cleanup();
        md.fragSort = null;

        log.info("Sorting list of duplicate records.");
        md.duplicateIndexes.doneAddingStartIteration();
        if (md.opticalDuplicateIndexes != null) {
            md.opticalDuplicateIndexes.doneAddingStartIteration();
        }
        if (md.TAG_DUPLICATE_SET_MEMBERS) {
            md.representativeReadIndicesForDuplicates.doneAdding();
        }
    }

    /**
     * Builds a read ends object that represents a single read - for flow based read
     */
    @Override
    public ReadEndsForMarkDuplicates buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec, final boolean useBarcodes) {
        final ReadEndsForMarkDuplicates ends = md.buildReadEnds(header, index, rec, useBarcodes);

        // this code only supported unpaired reads
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            throw new IllegalArgumentException("FLOW_MODE does not support paired reads. offending read: " + rec);
        }

        // adjust start/end coordinates
        ends.read1Coordinate = getReadEndCoordinate(rec, !rec.getReadNegativeStrandFlag(), true, md.flowBasedArguments);
        if (md.flowBasedArguments.USE_END_IN_UNPAIRED_READS) {
            ends.read2Coordinate = getReadEndCoordinate(rec, rec.getReadNegativeStrandFlag(), false, md.flowBasedArguments);
        }

        // adjust score
        if ( md.flowBasedArguments.FLOW_QUALITY_SUM_STRATEGY ) {
            ends.score = computeFlowDuplicateScore(rec, ends.read2Coordinate);
        }

        return ends;
    }

    /**
     * update score for pairedEnds
     */
    @Override
    public short getReadDuplicateScore(final SAMRecord rec, final ReadEndsForMarkDuplicates pairedEnds) {
        if (md.flowBasedArguments.FLOW_QUALITY_SUM_STRATEGY ) {
            return computeFlowDuplicateScore(rec, pairedEnds.read2Coordinate);
        } else {
            return md.getReadDuplicateScore(rec, pairedEnds);
        }
    }

    /**
     * This method is identical in function to areComparableForDuplicates except that it accomodates for
     * the possible significance of the end side of the reads (w/ or wo/ uncertainty). This is only
     * applicable for flow mode invocation.
     */
    private boolean areComparableForDuplicatesWithEndSignificance(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs, final boolean useBarcodes,
                                                                  final int lhsRead1Coordinate2Min, final int lhsRead1Coordinate2Max) {
        boolean areComparable = md.areComparableForDuplicates(lhs, rhs, false, useBarcodes);

        if (areComparable) {
            areComparable = (!endCoorSignificant(lhs.read2Coordinate, rhs.read2Coordinate) ||
                    endCoorInRangeWithUncertainty(lhsRead1Coordinate2Min, lhsRead1Coordinate2Max,
                            rhs.read2Coordinate, md.flowBasedArguments.UNPAIRED_END_UNCERTAINTY));
        }

        return areComparable;
    }

    private boolean endCoorSignificant(final int lhsCoor, final int rhsCoor) {
        return lhsCoor != END_INSIGNIFICANT_VALUE && rhsCoor != END_INSIGNIFICANT_VALUE;
    }

    private boolean endCoorInRangeWithUncertainty(final int lhsCoorMin, final int lhsCoorMax, final int rhsCoor, final int uncertainty) {
        return (rhsCoor >= (lhsCoorMin - uncertainty)) && (rhsCoor <= (lhsCoorMax + uncertainty));
    }

    /**
     * A quality summing scoring strategy used for flow based reads.
     *
     * The method walks on the bases of the read, in the synthesis direction. For each base, the effective
     * quality value is defined as the value on the first base on the hmer to which the base belongs to. The score
     * is defined to be the sum of all effective values above a given threshold.
     *
     * @param rec - SAMRecord to get a score for
     * @param threshold - threshold above which effective quality is included
     * @return - calculated score (see method description)
     */
    static protected int getFlowSumOfBaseQualities(final SAMRecord rec, final int threshold) {
        int score = 0;

        // access qualities and bases
        final byte[] quals = rec.getBaseQualities();
        final byte[]  bases = rec.getReadBases();

        // create iteration range and direction
        final int startingOffset = !rec.getReadNegativeStrandFlag() ? 0 : bases.length;
        final int endOffset = !rec.getReadNegativeStrandFlag() ? bases.length : 0;
        final int  iterIncr = !rec.getReadNegativeStrandFlag() ? 1 : -1;

        // loop on bases, extract qual related to homopolymer from start of homopolymer
        byte lastBase = 0;
        byte effectiveQual = 0;
        for ( int i = startingOffset ; i != endOffset ; i += iterIncr ) {
            final byte base = bases[i];
            if ( base != lastBase ) {
                effectiveQual = quals[i];
            }
            if ( effectiveQual >= threshold) {
                score += effectiveQual;
            }
            lastBase = base;
        }

        return score;
    }

    private short computeFlowDuplicateScore(final SAMRecord rec, final int end) {

        if ( end == END_INSIGNIFICANT_VALUE)
            return -1;

        Short storedScore = (Short)rec.getTransientAttribute(ATTR_DUPLICATE_SCORE);
        if ( storedScore == null ) {
            short score = 0;

            score += (short) Math.min(getFlowSumOfBaseQualities(rec, md.flowBasedArguments.FLOW_EFFECTIVE_QUALITY_THRESHOLD), Short.MAX_VALUE / 2);

            score += rec.getReadFailsVendorQualityCheckFlag() ? (short) (Short.MIN_VALUE / 2) : 0;
            storedScore = score;
            rec.setTransientAttribute(ATTR_DUPLICATE_SCORE, storedScore);
        }

        return storedScore;
    }

    @VisibleForTesting
    protected static int getReadEndCoordinate(final SAMRecord rec, final boolean startEnd, final boolean certain, final MarkDuplicatesForFlowArgumentCollection flowBasedArguments) {
        final FlowOrder flowOrder = new FlowOrder(rec);
        final int unclippedCoor = startEnd ? rec.getUnclippedStart() : rec.getUnclippedEnd();
        final int alignmentCoor = startEnd ? rec.getAlignmentStart() : rec.getAlignmentEnd();

        // this code requires a valid flow order
        if ( flowOrder.isValid() ) {

            // simple case
            if ( flowBasedArguments.USE_UNPAIRED_CLIPPED_END ) {
                return alignmentCoor;
            }

            // "skipping" case
            if (certain && flowBasedArguments.FLOW_SKIP_FIRST_N_FLOWS != 0) {
                final byte[] bases = rec.getReadBases();
                byte hmerBase = startEnd ? bases[0] : bases[bases.length - 1];
                int hmersLeft = flowBasedArguments.FLOW_SKIP_FIRST_N_FLOWS;      // number of hmer left to trim

                // advance flow order to base
                while (flowOrder.current() != hmerBase) {
                    flowOrder.advance();
                    hmersLeft--;
                }

                int hmerSize;
                for (hmerSize = 1; hmerSize < bases.length; hmerSize++) {
                    if ((startEnd ? bases[hmerSize] : bases[bases.length - 1 - hmerSize]) != hmerBase) {
                        if (--hmersLeft <= 0) {
                            break;
                        } else {
                            hmerBase = startEnd ? bases[hmerSize] : bases[bases.length - 1 - hmerSize];
                            flowOrder.advance();
                            while (flowOrder.current() != hmerBase) {
                                flowOrder.advance();
                                hmersLeft--;
                            }
                            if (hmersLeft <= 0) {
                                break;
                            }
                        }
                    }
                }
                final int coor = unclippedCoor + (startEnd ? hmerSize : -hmerSize);
                return flowBasedArguments.USE_UNPAIRED_CLIPPED_END
                        ? (startEnd ? Math.max(coor, alignmentCoor) : Math.min(coor, alignmentCoor))
                        : coor;
            }

            // "know end" case
            if (flowBasedArguments.FLOW_Q_IS_KNOWN_END ? isAdapterClipped(rec) : isAdapterClippedWithQ(rec)) {
                return unclippedCoor;
            }

            // "uncertain quality clipped" case
            if (!certain && isQualityClipped(rec)) {
                return END_INSIGNIFICANT_VALUE;
            }
        }

        // if here, return a default
        return unclippedCoor;
    }

    public static boolean isAdapterClipped(final SAMRecord rec) {
        return clippingTagContainsAny(rec, CLIPPING_TAG_CONTAINS_A);
    }

    public static boolean isAdapterClippedWithQ(final SAMRecord rec) {
        return clippingTagContainsAny(rec, CLIPPING_TAG_CONTAINS_AQ);
    }

    public static boolean isQualityClipped(final SAMRecord rec) {
        return clippingTagContainsAny(rec, CLIPPING_TAG_CONTAINS_QZ);
    }

    private static boolean clippingTagContainsAny(final SAMRecord rec, final char[] chars) {
        final String clippingTagValue = rec.getStringAttribute(CLIPPING_TAG_NAME);

        if ( clippingTagValue == null ) {
            return false;
        } else {
            for ( final char ch : chars ) {
                if ( clippingTagValue.indexOf(ch) >= 0 ) {
                    return true;
                }
            }
            return false;
        }
    }

    /**
     * private class used to represent use a SAMRecord's flow order, if such is present
     */
    static private class FlowOrder {

        byte[] flowOrder; // the flow order byte string
        int flowIndex = 0; // the current position on the flow order

        private FlowOrder(final SAMRecord rec) {

            // access flow order from record's read group
            if ( rec.getReadGroup() != null && rec.getReadGroup().getFlowOrder() != null ) {
                flowOrder = rec.getReadGroup().getFlowOrder().getBytes();
                return;
            }

            // fallback on finding a flow order elsewhere
            final SAMFileHeader header = rec.getHeader();
            for ( final SAMReadGroupRecord rg : header.getReadGroups() ) {
                if (rg.getFlowOrder() != null) {
                    flowOrder = rg.getFlowOrder().getBytes();
                    return;
                }
            }

            // otherwise, no flow order
            flowOrder = null;
        }

        private boolean isValid() {
            return flowOrder != null;
        }

        private void advance() {
            if (++flowIndex >= flowOrder.length) {
                flowIndex = 0;
            }
        }

        private byte current() {
            return flowOrder[flowIndex];
        }
    }
}
