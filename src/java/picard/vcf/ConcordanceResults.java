package picard.vcf;

import htsjdk.samtools.util.Histogram;

import java.util.EnumSet;

/**
 * A small class to store results of GenotypeConcordance between two samples
 *
 * @author Tim Fennell
 * @author George Grant
 */
class ConcordanceResults {
    private final Histogram<ConcordanceState> counts = new Histogram<ConcordanceState>();
    private final String truthSample;
    private final String callSample;

    ConcordanceResults(final String truthSample, final String callSample) {
        this.truthSample = truthSample;
        this.callSample = callSample;
    }

    void add(final VariantCallState truthState, final VariantCallState callState, final boolean altAllelesAgree) {
        this.counts.increment(new ConcordanceState(truthState, callState, altAllelesAgree));
    }

    long getCount(final VariantCallState truthState, final VariantCallState callState, final boolean altAllelesAgree) {
        final Histogram<ConcordanceState>.Bin bin = this.counts.get(new ConcordanceState(truthState, callState, altAllelesAgree));
        if (bin == null) return 0;
        else return (long) bin.getValue();
    }

    public String getTruthSample() { return truthSample; }
    public String getCallSample() { return callSample; }

    double hetSensitivity() {
        return getSensitivity(EnumSet.of(VariantCallState.Het));
    }

    double homVarSensitivity() {
        return getSensitivity(EnumSet.of(VariantCallState.HomVar));
    }

    double varSensitivity() {
        return getSensitivity(EnumSet.of(VariantCallState.Het, VariantCallState.HomVar));
    }

    private double getSensitivity(final EnumSet<VariantCallState> variantCallStates) {
        double numerator = 0.0;
        // The numerator is where call agrees with truth.  Sum it over all possible call states (i.e. Both Het and HomVar)
        for (final VariantCallState callState : variantCallStates) {
            numerator += getCount(callState, callState, true);
        }
        return numerator / sum(variantCallStates, EnumSet.allOf(VariantCallState.class), null);
    }

    double hetSpecificity() {
        return getSpecificity(EnumSet.of(VariantCallState.Het));
    }

    double homVarSpecificity() {
        return getSpecificity(EnumSet.of(VariantCallState.HomVar));
    }

    double varSpecificity() {
        return getSpecificity(EnumSet.of(VariantCallState.Het, VariantCallState.HomVar));
    }

    private double getSpecificity(final EnumSet<VariantCallState> variantCallStates) {
        final EnumSet<VariantCallState> otherVariantCallStates = EnumSet.allOf(VariantCallState.class);
        for (final VariantCallState variantCallState : variantCallStates) {
            otherVariantCallStates.remove(variantCallState);
        }
        return (double) sum(otherVariantCallStates, otherVariantCallStates, true) / sum(otherVariantCallStates, EnumSet.allOf(VariantCallState.class), null);
    }

    double hetPpv() {
        return getPpv(EnumSet.of(VariantCallState.Het));
    }

    double homVarPpv() {
        return getPpv(EnumSet.of(VariantCallState.HomVar));
    }

    double varPpv() {
        return getPpv(EnumSet.of(VariantCallState.Het, VariantCallState.HomVar));
    }

    private double getPpv(final EnumSet<VariantCallState> variantCallStates) {
        return (double) sum(variantCallStates, variantCallStates, true) / sum(EnumSet.allOf(VariantCallState.class), variantCallStates, null);
    }

    long numFalsePositives() {
        return numHomRefFalsePositives() + numHomVarFalsePositives() + numHetFalsePositives();
    }

    /**
     * Returns the number of HomRef False positives.
     * i.e. Where truth is HomRef and Call is HomVar or Het
     * @return
     */
    long numHomRefFalsePositives() {
        return sum(EnumSet.of(VariantCallState.HomRef), EnumSet.of(VariantCallState.Het, VariantCallState.HomVar), true);
    }

    /**
     * Returns the number of HomVar False positives.
     * i.e. Where truth is HomVar and Call is HomRef or Het
     * @return
     */
    long numHomVarFalsePositives() {
        return sum(EnumSet.of(VariantCallState.HomVar), EnumSet.of(VariantCallState.Het, VariantCallState.HomRef), true);
    }

    /**
     * Returns the number of Het False positives.
     * i.e. Where truth is Het and Call is HomRef or HomVar
     * @return
     */
    long numHetFalsePositives() {
        return sum(EnumSet.of(VariantCallState.Het), EnumSet.of(VariantCallState.HomVar, VariantCallState.HomRef), true);
    }

    /** Sums the counts where the first state is contains in truthStateSet and the second state is contained in callStateSet. */
    private long sum(final EnumSet<VariantCallState> truthStateSet, final EnumSet<VariantCallState> callStateSet, final Boolean altAllelesAgree) {
        long result = 0;
        for (final VariantCallState truthState : truthStateSet) {
            for (final VariantCallState callState : callStateSet) {
                if (altAllelesAgree == null) {
                    result += getCount(truthState, callState, true);
                    result += getCount(truthState, callState, false);
                } else {
                    result += getCount(truthState, callState, altAllelesAgree);
                }
            }
        }
        return result;
    }

}
