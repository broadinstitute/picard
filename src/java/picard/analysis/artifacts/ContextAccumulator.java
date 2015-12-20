package picard.analysis.artifacts;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ListMap;
import htsjdk.samtools.util.SequenceUtil;
import picard.PicardException;
import picard.analysis.artifacts.SequencingArtifactMetrics.BaitBiasDetailMetrics;
import picard.analysis.artifacts.SequencingArtifactMetrics.DetailPair;
import picard.analysis.artifacts.SequencingArtifactMetrics.PreAdapterDetailMetrics;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Keeps track of the AlignmentAccumulators for each artifact / context of interest.
 */
class ContextAccumulator {

    // are the PE reads expected to face the same direction?
    private final boolean expectedTandemReads;

    // mapping from contexts to the accumulators
    private final Map<String, AlignmentAccumulator[]> artifactMap;

    public ContextAccumulator(final Set<String> contexts, final boolean expectedTandemReads) {
        this.expectedTandemReads = expectedTandemReads;
        this.artifactMap = new HashMap<>();
        for (final String context : contexts) {

            // sanity check that the context length is odd
            if ((context.length() & 1) == 0) throw new PicardException("Contexts cannot have an even number of bases: " + context);

            final AlignmentAccumulator[] accumulators = new AlignmentAccumulator[Transition.Base.values().length];
            for (int i = 0; i < Transition.Base.values().length; i++) {
                accumulators[i] = new AlignmentAccumulator();
            }
            this.artifactMap.put(context, accumulators);
        }
    }

    public void countRecord(final String refContext, final char calledBase, final SAMRecord rec) {
        artifactMap.get(refContext)[Transition.baseIndexMap[calledBase]].countRecord(rec);
    }

    /**
     * Core method to compute detailed (i.e. context-by-context) metrics from this accumulator.
     */
    public ListMap<Transition, DetailPair> calculateMetrics(final String sampleAlias, final String library) {
        final ListMap<Transition, DetailPair> detailMetricsMap = new ListMap<>();
        for (final String context : new TreeSet<>(artifactMap.keySet())) {

            // sanity check that the context length is odd
            if ((context.length() & 1) == 0) throw new PicardException("Contexts cannot have an even number of bases: " + context + ".  This should never happen here!");
            final char refBase = context.charAt(context.length() / 2);

            for (final Transition.Base altBase : Transition.Base.values()) {
                final Transition transition = Transition.transitionOf(refBase, (char)altBase.base);

                // each combination of artifact + context represents a single metric row
                final PreAdapterDetailMetrics preAdapterDetailMetrics = new PreAdapterDetailMetrics();
                final BaitBiasDetailMetrics baitBiasDetailMetrics = new BaitBiasDetailMetrics();

                // populate basic fields
                preAdapterDetailMetrics.SAMPLE_ALIAS = sampleAlias;
                preAdapterDetailMetrics.LIBRARY = library;
                preAdapterDetailMetrics.CONTEXT = context;
                preAdapterDetailMetrics.REF_BASE = transition.ref();
                preAdapterDetailMetrics.ALT_BASE = transition.call();

                baitBiasDetailMetrics.SAMPLE_ALIAS = sampleAlias;
                baitBiasDetailMetrics.LIBRARY = library;
                baitBiasDetailMetrics.CONTEXT = context;
                baitBiasDetailMetrics.REF_BASE = transition.ref();
                baitBiasDetailMetrics.ALT_BASE = transition.call();

                // retrieve all the necessary alignment counters.
                final AlignmentAccumulator[] accumulators = artifactMap.get(context);
                final AlignmentAccumulator[] reverseCompAccumulators = artifactMap.get(SequenceUtil.reverseComplement(context));

                final AlignmentAccumulator fwdRefAlignments = accumulators[Transition.baseIndexMap[transition.ref()]];
                final AlignmentAccumulator fwdAltAlignments = accumulators[Transition.baseIndexMap[transition.call()]];
                final AlignmentAccumulator revRefAlignments = reverseCompAccumulators[Transition.baseIndexMap[transition.complement().ref()]];
                final AlignmentAccumulator revAltAlignments = reverseCompAccumulators[Transition.baseIndexMap[transition.complement().call()]];

                // categorize observations of pre-adapter artifacts
                if (expectedTandemReads) {
                    // if both ends are sequenced on the same strand, then read1/read2 should exhibit the same bias
                    preAdapterDetailMetrics.PRO_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R2_POS + revRefAlignments.R1_NEG + revRefAlignments.R2_NEG;
                    preAdapterDetailMetrics.PRO_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R2_POS + revAltAlignments.R1_NEG + revAltAlignments.R2_NEG;
                    preAdapterDetailMetrics.CON_REF_BASES = fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_NEG + revRefAlignments.R1_POS + revRefAlignments.R2_POS;
                    preAdapterDetailMetrics.CON_ALT_BASES = fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_NEG + revAltAlignments.R1_POS + revAltAlignments.R2_POS;
                } else {
                    // if ends are sequenced on opposite strands, then read1/read2 should exhibit opposite biases
                    preAdapterDetailMetrics.PRO_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R2_NEG + revRefAlignments.R1_NEG + revRefAlignments.R2_POS;
                    preAdapterDetailMetrics.PRO_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R2_NEG + revAltAlignments.R1_NEG + revAltAlignments.R2_POS;
                    preAdapterDetailMetrics.CON_REF_BASES = fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + revRefAlignments.R1_POS + revRefAlignments.R2_NEG;
                    preAdapterDetailMetrics.CON_ALT_BASES = fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + revAltAlignments.R1_POS + revAltAlignments.R2_NEG;
                }

                // categorize observations of bait bias artifacts
                baitBiasDetailMetrics.FWD_CXT_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + fwdRefAlignments.R2_NEG;
                baitBiasDetailMetrics.FWD_CXT_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + fwdAltAlignments.R2_NEG;
                baitBiasDetailMetrics.REV_CXT_REF_BASES = revRefAlignments.R1_POS + revRefAlignments.R1_NEG + revRefAlignments.R2_POS + revRefAlignments.R2_NEG;
                baitBiasDetailMetrics.REV_CXT_ALT_BASES = revAltAlignments.R1_POS + revAltAlignments.R1_NEG + revAltAlignments.R2_POS + revAltAlignments.R2_NEG;

                // calculate error rates + Q-scores
                preAdapterDetailMetrics.calculateDerivedStatistics();
                baitBiasDetailMetrics.calculateDerivedStatistics();

                // add the finalized metrics to the map
                detailMetricsMap.add(transition, new DetailPair(preAdapterDetailMetrics, baitBiasDetailMetrics));
            }
        }
        return detailMetricsMap;
    }

    /**
     * Little class for breaking down alignments by read1/read2 and positive/negative strand.
     */
    private static class AlignmentAccumulator {
        private long R1_POS = 0;
        private long R1_NEG = 0;
        private long R2_POS = 0;
        private long R2_NEG = 0;

        private void countRecord(final SAMRecord rec) {
            final boolean isNegativeStrand = rec.getReadNegativeStrandFlag();
            final boolean isReadTwo = rec.getReadPairedFlag() && rec.getSecondOfPairFlag();
            if (isReadTwo) {
                if (isNegativeStrand) this.R2_NEG++;
                else this.R2_POS++;
            } else {
                if (isNegativeStrand) this.R1_NEG++;
                else this.R1_POS++;
            }
        }
    }
}
