package picard.vcf;

import htsjdk.samtools.util.Histogram;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * A small class to store counts of co-occurring VariantContext.Types
 *
 * @author George Grant
 */
public class VariantContextTypeResults {
    private final Histogram<VariantContextTypeState> counts = new Histogram<VariantContextTypeState>();
    private final String truthSample, callSample;

    VariantContextTypeResults(final String truthSample, final String callSample) {
        this.truthSample = truthSample;
        this.callSample = callSample;
    }

    void add(final VariantContext.Type truthType, final VariantContext.Type callType) {
        this.counts.increment(new VariantContextTypeState(truthType, callType));
    }

    long getCount(final VariantContext.Type truthType, final VariantContext.Type callType) {
        final Histogram<VariantContextTypeState>.Bin bin = this.counts.get(new VariantContextTypeState(truthType, callType));
        if (bin == null) return 0;
        else return (long) bin.getValue();
    }

    public String getTruthSample() { return truthSample; }
    public String getCallSample() { return callSample; }
}

