/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package picard.vcf.filter;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;

import java.util.List;

/**
 * Filters records based on the phred scaled p-value from the Fisher Strand test stored in
 * the FS attribute.
 *
 * @author tfennell
 */
public class FisherStrandFilter implements VariantFilter {
    private final double maxPhredScalePValue;

    public FisherStrandFilter(final double maxPhredScalePValue) {
        this.maxPhredScalePValue= maxPhredScalePValue;
    }

    @Override
    public List<VCFFilterHeaderLine> headerLines() {
        return CollectionUtil.makeList(new VCFFilterHeaderLine("StrandBias", "Site exhibits excessive allele/strand correlation."));
    }

    @Override
    public String filter(final VariantContext ctx) {
        final double fs = ctx.getAttributeAsDouble("FS", 0);
        return (fs > maxPhredScalePValue) ? "StrandBias" : null;
    }
}
