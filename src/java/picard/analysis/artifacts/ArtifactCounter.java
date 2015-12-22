package picard.analysis.artifacts;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ListMap;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.analysis.artifacts.SequencingArtifactMetrics.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Keeps track of artifact counts, and extracts metrics once accumulation is finished.
 */
class ArtifactCounter {
    private final String sampleAlias;
    private final String library;

    private final Map<String, RefContext> contextMap = new HashMap<>();

    private final ContextAccumulator fullContextAccumulator;
    private final ContextAccumulator halfContextAccumulator;
    private final ContextAccumulator zeroContextAccumulator;

    private final List<PreAdapterSummaryMetrics> preAdapterSummaryMetricsList;
    private final List<PreAdapterDetailMetrics> preAdapterDetailMetricsList;
    private final List<BaitBiasSummaryMetrics> baitBiasSummaryMetricsList;
    private final List<BaitBiasDetailMetrics> baitBiasDetailMetricsList;

    private final Set<String> leadingContexts = new HashSet<>();
    private final Set<String> trailingContexts = new HashSet<>();

    // tuple to keep track of the different types of sub-contexts from a given reference context
    protected final class RefContext {
        final String ref, leading, trailing, zero;

        public RefContext(final String ref, final String leading, final String trailing, final String zero) {
            this.ref = ref;
            this.leading = leading;
            this.trailing = trailing;
            this.zero = zero;
        }
    }

    public ArtifactCounter(final String sampleAlias, final String library, final int contextSize, final boolean expectedTandemReads) {
        this.sampleAlias = sampleAlias;
        this.library = library;

        // define the contexts
        final HashSet<String> fullContexts = new HashSet<>();
        for (final byte[] kmer : SequenceUtil.generateAllKmers(2 * contextSize + 1)) {
            fullContexts.add(StringUtil.bytesToString(kmer));
        }

        final Set<String> zeroContexts = new HashSet<>();

        // the half contexts specify either leading or trailing bases. the zero context is just the center.
        // NB: we use N to represent a wildcard base, rather than an ambiguous base. It's assumed that all of the input
        // contexts are unambiguous, and that any actual N's in the data have been dealt with elsewhere.
        final String padding = StringUtil.repeatCharNTimes('N', contextSize);
        for (final String context : fullContexts) {
            final char centralBase = context.charAt(contextSize);
            final String leading = context.substring(0, contextSize) + centralBase + padding;
            final String trailing = padding + centralBase + context.substring(contextSize + 1, context.length());
            final String zero = padding + centralBase + padding;
            contextMap.put(context, new RefContext(context, leading, trailing, zero));

            leadingContexts.add(leading);
            trailingContexts.add(trailing);
            zeroContexts.add(zero);
        }

        final Set<String> halfContexts = new HashSet<>(leadingContexts);
        halfContexts.addAll(trailingContexts);

        this.fullContextAccumulator = new ContextAccumulator(fullContexts, expectedTandemReads);
        this.halfContextAccumulator = new ContextAccumulator(halfContexts, expectedTandemReads);
        this.zeroContextAccumulator = new ContextAccumulator(zeroContexts, expectedTandemReads);

        // these will get populated in the final step
        preAdapterSummaryMetricsList = new ArrayList<PreAdapterSummaryMetrics>();
        preAdapterDetailMetricsList = new ArrayList<PreAdapterDetailMetrics>();
        baitBiasSummaryMetricsList = new ArrayList<BaitBiasSummaryMetrics>();
        baitBiasDetailMetricsList = new ArrayList<BaitBiasDetailMetrics>();
    }

    /**
     * Add a record to all the accumulators.
     */
    public void countRecord(final String refContext, final char calledBase, final SAMRecord rec) {
        if (this.contextMap.containsKey(refContext)) {
            final RefContext contexts = contextMap.get(refContext);
            this.fullContextAccumulator.countRecord(contexts.ref, calledBase, rec);
            this.halfContextAccumulator.countRecord(contexts.leading, calledBase, rec);
            this.halfContextAccumulator.countRecord(contexts.trailing, calledBase, rec);
            this.zeroContextAccumulator.countRecord(contexts.zero, calledBase, rec);
        }
    }

    /**
     * Stop counting, tally things up, and extract metrics.
     */
    public void finish() {
        final ListMap<Transition, DetailPair> allDetailMetrics = getDetailMetrics();
        final Map<Transition, SummaryPair> allSummaryMetrics = getSummaryMetrics();

        for (final Transition transition : Transition.altValues()) {
            final SummaryPair summary = allSummaryMetrics.get(transition);
            final List<DetailPair> details = allDetailMetrics.get(transition);
            preAdapterSummaryMetricsList.add(summary.preAdapterMetrics);
            baitBiasSummaryMetricsList.add(summary.baitBiasMetrics);
            for (final DetailPair detail : details) {
                preAdapterDetailMetricsList.add(detail.preAdapterMetrics);
                baitBiasDetailMetricsList.add(detail.baitBiasMetrics);
            }
        }
    }

    public List<PreAdapterSummaryMetrics> getPreAdapterSummaryMetrics() { return preAdapterSummaryMetricsList; }
    public List<PreAdapterDetailMetrics> getPreAdapterDetailMetrics() { return preAdapterDetailMetricsList; }
    public List<BaitBiasSummaryMetrics> getBaitBiasSummaryMetrics() { return baitBiasSummaryMetricsList; }
    public List<BaitBiasDetailMetrics> getBaitBiasDetailMetrics() { return baitBiasDetailMetricsList; }

    /**
     * Core method to compute summary metrics. For each transition, we report:
     * 1. the total Q-score across all contexts
     * 2. the worst full context and its Q-score
     * 3. the worst leading context and its Q-score
     * 4. the worst trailing context and its Q-score
     *
     */
    private Map<Transition, SummaryPair> getSummaryMetrics() {
        final Map<Transition, SummaryPair> summaryMetricsMap = new HashMap<Transition, SummaryPair>();

        // extract the detail metrics from each accumulator
        final ListMap<Transition, DetailPair> fullMetrics = this.fullContextAccumulator.calculateMetrics(sampleAlias, library);
        final ListMap<Transition, DetailPair> halfMetrics = this.halfContextAccumulator.calculateMetrics(sampleAlias, library);
        final ListMap<Transition, DetailPair> zeroMetrics = this.zeroContextAccumulator.calculateMetrics(sampleAlias, library);

        // compute the summary metrics - one row for each transition
        for (final Transition transition : Transition.altValues()) {
            final List<DetailPair> fullMetricsForTransition = fullMetrics.get(transition);
            final List<DetailPair> zeroMetricsForTransition = zeroMetrics.get(transition);
            if (zeroMetricsForTransition.size() != 1) {
                throw new PicardException("Should have exactly one context-free metric pair for transition: " + transition);
            }

            // we want to report on leading / trailing contexts separately
            final List<DetailPair> leadingMetricsForTransition = new ArrayList<DetailPair>();
            final List<DetailPair> trailingMetricsForTransition = new ArrayList<DetailPair>();
            for (final DetailPair metrics : halfMetrics.get(transition)) {
                // first make sure they're the same context
                if (!metrics.preAdapterMetrics.CONTEXT.equals(metrics.baitBiasMetrics.CONTEXT)) {
                    throw new PicardException("Input detail metrics are not matched up properly - contexts differ.");
                }
                final boolean isLeading = leadingContexts.contains(metrics.preAdapterMetrics.CONTEXT);
                final boolean isTrailing = trailingContexts.contains(metrics.preAdapterMetrics.CONTEXT);
                // if the original contextSize is 0, there's no difference between leading and trailing, so add it to both
                if (isLeading) leadingMetricsForTransition.add(metrics);
                if (isTrailing) trailingMetricsForTransition.add(metrics);
            }

            // get the worst cases
            final DetailPair totalMetric = zeroMetricsForTransition.get(0);
            final DetailPair worstFullMetric = getWorstMetrics(fullMetricsForTransition);
            final DetailPair worstLeadingMetric = getWorstMetrics(leadingMetricsForTransition);
            final DetailPair worstTrailingMetric = getWorstMetrics(trailingMetricsForTransition);

            // construct the actual summary metrics - a combination of all the data we've just extracted
            final PreAdapterSummaryMetrics preAdapterSummaryMetrics = new PreAdapterSummaryMetrics();
            final BaitBiasSummaryMetrics baitBiasSummaryMetrics = new BaitBiasSummaryMetrics();

            preAdapterSummaryMetrics.SAMPLE_ALIAS = this.sampleAlias;
            preAdapterSummaryMetrics.LIBRARY = this.library;
            preAdapterSummaryMetrics.REF_BASE = transition.ref();
            preAdapterSummaryMetrics.ALT_BASE = transition.call();
            preAdapterSummaryMetrics.TOTAL_QSCORE = totalMetric.preAdapterMetrics.QSCORE;
            preAdapterSummaryMetrics.WORST_CXT = worstFullMetric.preAdapterMetrics.CONTEXT;
            preAdapterSummaryMetrics.WORST_CXT_QSCORE = worstFullMetric.preAdapterMetrics.QSCORE;
            preAdapterSummaryMetrics.WORST_PRE_CXT = worstLeadingMetric.preAdapterMetrics.CONTEXT;
            preAdapterSummaryMetrics.WORST_PRE_CXT_QSCORE = worstLeadingMetric.preAdapterMetrics.QSCORE;
            preAdapterSummaryMetrics.WORST_POST_CXT = worstTrailingMetric.preAdapterMetrics.CONTEXT;
            preAdapterSummaryMetrics.WORST_POST_CXT_QSCORE = worstTrailingMetric.preAdapterMetrics.QSCORE;
            preAdapterSummaryMetrics.inferArtifactName();

            baitBiasSummaryMetrics.SAMPLE_ALIAS = this.sampleAlias;
            baitBiasSummaryMetrics.LIBRARY = this.library;
            baitBiasSummaryMetrics.REF_BASE = transition.ref();
            baitBiasSummaryMetrics.ALT_BASE = transition.call();
            baitBiasSummaryMetrics.TOTAL_QSCORE = totalMetric.baitBiasMetrics.QSCORE;
            baitBiasSummaryMetrics.WORST_CXT = worstFullMetric.baitBiasMetrics.CONTEXT;
            baitBiasSummaryMetrics.WORST_CXT_QSCORE = worstFullMetric.baitBiasMetrics.QSCORE;
            baitBiasSummaryMetrics.WORST_PRE_CXT = worstLeadingMetric.baitBiasMetrics.CONTEXT;
            baitBiasSummaryMetrics.WORST_PRE_CXT_QSCORE = worstLeadingMetric.baitBiasMetrics.QSCORE;
            baitBiasSummaryMetrics.WORST_POST_CXT = worstTrailingMetric.baitBiasMetrics.CONTEXT;
            baitBiasSummaryMetrics.WORST_POST_CXT_QSCORE = worstTrailingMetric.baitBiasMetrics.QSCORE;
            baitBiasSummaryMetrics.inferArtifactName();

            // add the finalized metrics to the map
            summaryMetricsMap.put(transition, new SummaryPair(preAdapterSummaryMetrics, baitBiasSummaryMetrics));
        }
        return summaryMetricsMap;
    }

    private ListMap<Transition, DetailPair> getDetailMetrics() {
        return this.fullContextAccumulator.calculateMetrics(this.sampleAlias, this.library);
    }

    /**
     * Given a list of detail metrics, get the worst pre-adapter metrics, and independently from that get the worst bait bias metrics
     * (in terms of Q-score).
     */
    private DetailPair getWorstMetrics(final List<DetailPair> metrics) {
        PreAdapterDetailMetrics worstPreAdapterMetrics = null;
        BaitBiasDetailMetrics worstBaitBiasMetrics = null;
        for (final DetailPair m : metrics) {

            //The comparator first comparse by QSCORE and then uses other fields to guarrantee a deterministic order
            if (worstPreAdapterMetrics == null || m.preAdapterMetrics.compareTo(worstPreAdapterMetrics) < 0) worstPreAdapterMetrics = m.preAdapterMetrics;
            if (worstBaitBiasMetrics   == null || m.baitBiasMetrics.  compareTo(worstBaitBiasMetrics)   < 0) worstBaitBiasMetrics   = m.baitBiasMetrics;
        }
        return new DetailPair(worstPreAdapterMetrics, worstBaitBiasMetrics);
    }
}
