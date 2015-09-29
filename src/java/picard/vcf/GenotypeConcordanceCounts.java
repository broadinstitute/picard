package picard.vcf;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.Histogram;
import picard.PicardException;
import picard.vcf.GenotypeConcordanceStates.*;

/**
 * A class to store the counts for various truth and call state classifications relative to a reference.  With these counts and a provided
 * scheme, summary metrics can be returned.
 * @author nhomer
 */
public class GenotypeConcordanceCounts {

    /**
     * Pre-defined sets based on if the caller wishes to return the sensitivity given the common homozygous reference, heterozygous, and homozygous variant cases.
     * Var truth states can also be used to calculate specificity. Hom Ref and Missing can be added to var truth states because they are needed
     * for specificity and do not contain TP or FN, so they will not effect sensitivity.
     */
    static final TruthState[] HOM_REF_TRUTH_STATES = {TruthState.HOM_REF};
    static final TruthState[] HET_TRUTH_STATES = {TruthState.HET_REF_VAR1, TruthState.HET_VAR1_VAR2};
    static final TruthState[] HOM_VAR_TRUTH_STATES = {TruthState.HOM_VAR1};
    static final TruthState[] VAR_TRUTH_STATES = {TruthState.HET_REF_VAR1, TruthState.HET_VAR1_VAR2, TruthState.HOM_VAR1,
            TruthState.HOM_REF, TruthState.MISSING};

    /**
     * Pre-defined sets based on if the caller wishes to return the PPV given the common homozygous reference, heterozygous, and homozygous variant cases.
     */
    static final CallState[] HOM_REF_CALL_STATES = {CallState.HOM_REF};
    static final CallState[] HET_CALL_STATES = {CallState.HET_REF_VAR1, CallState.HET_REF_VAR2, CallState.HET_REF_VAR3,
            CallState.HET_VAR1_VAR2, CallState.HET_VAR1_VAR3, CallState.HET_VAR3_VAR4};
    static final CallState[] HOM_VAR_CALL_STATES = {CallState.HOM_VAR1, CallState.HOM_VAR2, CallState.HOM_VAR3};
    static final CallState[] VAR_CALL_STATES = {CallState.HET_REF_VAR1, CallState.HET_REF_VAR2, CallState.HET_REF_VAR3,
            CallState.HET_VAR1_VAR2, CallState.HET_VAR1_VAR3, CallState.HET_VAR3_VAR4,
            CallState.HOM_VAR1, CallState.HOM_VAR2, CallState.HOM_VAR3};

    /** The underlying counts table */
    private final Histogram<TruthAndCallStates> counter = new Histogram<TruthAndCallStates>();

    /**
     * Increments a count for the truth/call state tuple.
     * @param truthAndCallStates
     */
    public void increment(final TruthAndCallStates truthAndCallStates) {
        this.counter.increment(truthAndCallStates);
    }

    public void increment(final TruthAndCallStates truthAndCallStates, final double count){
        this.counter.increment(truthAndCallStates, count);
    }

    public double getCounterSize() {
        return this.counter.getCount();
    }

    /**
     * Validates that there are no counts for NA states in the underlying scheme
     */
    public void validateCountsAgainstScheme(final GenotypeConcordanceScheme scheme) {
        final Set<ContingencyState> naContingencyStates = getContingencyStateSet(GenotypeConcordanceScheme.NA);
        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                if (0 < getCount(truthAndCallStates)) {
                    final Set<ContingencyState> contingencyStates = getContingencyStateSet(scheme.getConcordanceStateArray(truthAndCallStates));
                    if (contingencyStates.containsAll(naContingencyStates)) {
                        throw new PicardException(String.format("Found counts for an illegal set of states: [%s, %s]", truthState.name(), callState.name()));
                    }
                }
            }
        }
    }

    private Set<ContingencyState> getContingencyStateSet(final ContingencyState[] contingencyStateArray) {
        final Set<ContingencyState> contingencyStateSet = new HashSet<ContingencyState>();
        Collections.addAll(contingencyStateSet, contingencyStateArray);
        return contingencyStateSet;
    }

    /**
     * Genotype Concordance Util calculates the Genotype Concordance and the Non-Ref Genotype Concordance, based on the includeHomRef flag.
     */
    private double calculateGenotypeConcordanceUtil(final GenotypeConcordanceScheme scheme, final boolean missingSitesFlag, final boolean includeHomRef) {
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        final TruthState[] allTruthStates = TruthState.values();
        final CallState[] allCallStates = CallState.values();

        for (final TruthState truthState : allTruthStates) {
            for (final CallState callState : allCallStates) {
                if (!missingSitesFlag && isMissing(truthState, callState)) {
                    continue;
                }
                else if (includeHomRef || isVar(truthState, callState)) {
                    final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                    final long count = getCount(truthAndCallStates);
                    if (truthState.getCode()==callState.getCode()) {
                        //If we enter this, we are 'on the diagonal'
                        numerator += count;
                    }
                    denominator += count;
                }
            }
        }
        return (denominator > 0.0 ? (numerator / denominator) : Double.NaN);
    }

    /**
     * Genotype Concordance is the number of times the truth and call states match exactly / all truth and call combinations made
     * If the GA4GH scheme is being used, any MISSING sites in truth OR call will not be included in the discordance calculations.
     */
    public double calculateGenotypeConcordance(final GenotypeConcordanceScheme scheme, final boolean missingSitesFlag) {
        return calculateGenotypeConcordanceUtil(scheme, missingSitesFlag, true);
    }

    /**
     * Non Ref Genotype Concordance is the number of times the truth and call states match exactly for *vars only* / all truth and call *var* combinations made
     * If the GA4GH scheme is being used, any MISSING sites in truth OR call will not be included in the discordance calculations.
     */
    public double calculateNonRefGenotypeConcordance(final GenotypeConcordanceScheme scheme, final boolean missingSitesFlag) {
        return calculateGenotypeConcordanceUtil(scheme, missingSitesFlag, false);
    }

    /**
     * Returns the sensitivity defined by the scheme across the subset of truth states.
     */
    public double getSensitivity(final GenotypeConcordanceScheme scheme, final TruthState[] truthStateArray) {
        /**
         * Sensitivity is the TP / P = TP / (TP + FN)
         */
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        for (final TruthState truthState : truthStateArray) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                final long count = getCount(truthAndCallStates);
                for (final ContingencyState contingencyState : scheme.getConcordanceStateArray(truthAndCallStates)) {
                    if (ContingencyState.TP == contingencyState) {
                        numerator += count;
                        denominator += count;
                    } else if (ContingencyState.FN == contingencyState) {
                        denominator += count;
                    }
                }
            }
        }

        return (numerator / denominator);
    }

    /**
     * Returns the PPV defined by the scheme across the subset of call states.
     */
    public double Ppv(final GenotypeConcordanceScheme scheme, final CallState[] callStateList) {
        /**
         * PPV is the TP / (TP + FP)
         */
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        for (final CallState callState : callStateList) {
            for (final TruthState truthState : TruthState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                final long count = getCount(truthAndCallStates);
                for (final ContingencyState contingencyState : scheme.getConcordanceStateArray(truthAndCallStates)) {
                    if (ContingencyState.TP == contingencyState) {
                        numerator += count;
                        denominator += count;
                    } else if (ContingencyState.FP == contingencyState) {
                        denominator += count;
                    }
                }
            }
        }

        return (numerator / denominator);
    }

    /**
     * Returns the specificity defined by the scheme across the subset of truth states.
     */
    public double getSpecificity(final GenotypeConcordanceScheme scheme, final TruthState[] truthStateArray) {
        /**
         * Specificity is the TN / N = TN / (FP + TN)
         */
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        for (final TruthState truthState : truthStateArray) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                final long count = getCount(truthAndCallStates);
                for (final ContingencyState contingencyState : scheme.getConcordanceStateArray(truthAndCallStates)) {
                    if (ContingencyState.TN == contingencyState) {
                        numerator += count;
                        denominator += count;
                    } else if (ContingencyState.FP == contingencyState) {
                        denominator += count;
                    }
                }
            }
        }
        return (numerator / denominator);
    }

    /**
     * Returns the count defined by the truth state set and call state set.
     */
    public long getCount(final TruthState truthState, final CallState callState) {
        return getCount(new TruthAndCallStates(truthState, callState));
    }

    /**
     * Returns the count defined by the truth state set and call state set.
     */
    public long getCount(final TruthAndCallStates truthAndCallStates) {
        final Histogram<TruthAndCallStates>.Bin bin = this.counter.get(truthAndCallStates);
        return (bin == null ? 0L : (long) bin.getValue());
    }

    /**
     * Returns true if EITHER the truth or call state is a VAR.
     * Used for calculating non ref genotype concordance.
     */
    public boolean isVar(final TruthState truthState, final CallState callState) {
        final List<TruthState> listOfTruthStates = Arrays.asList(TruthState.HOM_VAR1, TruthState.HET_REF_VAR1, TruthState.HET_VAR1_VAR2);
        final List<CallState> listOfCallStates = Arrays.asList(CallState.HET_REF_VAR1, CallState.HET_REF_VAR2, CallState.HET_REF_VAR3,
                CallState.HET_VAR1_VAR2, CallState.HET_VAR1_VAR3, CallState.HET_VAR3_VAR4,
                CallState.HOM_VAR1, CallState.HOM_VAR2, CallState.HOM_VAR3);

        return listOfTruthStates.contains(truthState) || listOfCallStates.contains(callState);
    }

    /**
     * Returns true if EITHER the truth or call state is MISSING.
     * Used for calculating genotype concordance and non-ref genotype concordance when the GA4GH scheme is used.
     */
    public boolean isMissing(final TruthState truthState, final CallState callState) {
        return truthState == TruthState.MISSING || callState == CallState.MISSING;
    }

    /**
     * Returns the sum of all pairs of tuples defined by the truth state set and call state set.
     */
    public long getSum(final Set<TruthState> truthStateSet, final Set<CallState> callStateSet) {
        long count = 0;
        for (final TruthState truthState : truthStateSet) {
            for (final CallState callState : callStateSet) {
                count += getCount(truthState, callState);
            }
        }
        return count;
    }

    /**
     * Returns the sum of all pairs of tuples defined by the truth state set and call state set.
     */
    public long getSum() {
        return getSum(new HashSet<TruthState>(Arrays.asList(TruthState.values())), new HashSet<CallState>(Arrays.asList(CallState.values())));
    }

    /**
     * Returns the total number of times each contingency state is encountered, summed across all truth/call state pairs.
     */
    public Map<ContingencyState, Long> getContingencyStateCounts(final GenotypeConcordanceScheme scheme) {
        scheme.validateScheme();

        final Map<ContingencyState, Long> counts = new HashMap<ContingencyState, Long>();
        for (final ContingencyState contingencyState : ContingencyState.values()) {
            counts.put(contingencyState, 0L);
        }

        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                final ContingencyState[] contingencyStateArray = scheme.getConcordanceStateArray(truthAndCallStates);
                for (final ContingencyState contingencyState : contingencyStateArray) {
                    final long newCount = counts.get(contingencyState) + getCount(truthAndCallStates);
                    counts.put(contingencyState, newCount);
                }
            }
        }

        return counts;
    }
}
