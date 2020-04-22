package picard.fingerprint;

import picard.util.MathUtil;

public class CappedHaplotypeProbabilities extends HaplotypeProbabilitiesUsingLogLikelihoods {
    private final double cap;

    // cap should be negative to indicate that you can never be too sure of anything (since the log likelihood is the
    //    log probability of an error, it's negative)
    public CappedHaplotypeProbabilities(final HaplotypeProbabilities haplotypeProbabilities, double cap) {
        super(haplotypeProbabilities.getHaplotype());
        this.cap = cap;
        final double[] logLikelihoods = haplotypeProbabilities.getLogLikelihoods();
        final double max = MathUtil.max(logLikelihoods);
        final double[] cappedLogLikelihoods = MathUtil.sum(logLikelihoods,-max);
        this.setLogLikelihoods(MathUtil.max(cappedLogLikelihoods,cap));
    }

    @Override
    public HaplotypeProbabilities deepCopy() {
        return new CappedHaplotypeProbabilities(this,cap);
    }
}
