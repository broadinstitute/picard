package picard.vcf;
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

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import picard.util.DbSnpBitSetUtil;

import java.util.List;

/**
 * An accumulator for collecting metrics about a single-sample GVCF. The main point here is to subset the
 * context of each {@link VariantContext} as it comes by to the alleles present in the genotype of the only sample.
 * Since this is a GVCF we expect a symbolic \<NON_REF\> allele to be present in each VC. If we do not subset
 * the context this symbolic allele will cause the regular {@link CallingMetricAccumulator} to return only a
 * small subset of the relevant metrics.
 *
 * @author farjoun
 */
public class GvcfMetricAccumulator extends CallingMetricAccumulator {
    String sample = null;

    public GvcfMetricAccumulator(final DbSnpBitSetUtil.DbSnpBitSets dbsnp) {
        super(dbsnp);
    }

    @Override
    public void setup(final VCFHeader vcfHeader) {
        final List<String> samples = vcfHeader.getGenotypeSamples();
        if (samples == null || samples.size() != 1) {
            throw new IllegalArgumentException("Expected to have exactly 1 sample in a GVCF, found " + ((samples == null) ? "0" : samples.size()));
        }
        sample = samples.get(0);
    }

    @Override
    public void accumulate(final VariantContext vc) {
        //since a gvcf always has a <NON_REF> allele, in order to get meaningful results we need to subset the context of
        // the variant to the alleles that actually appear in the only sample's genotype
        final VariantContext subContext = vc.subContextFromSample(sample);
        super.accumulate(subContext);
    }
}
