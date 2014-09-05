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
 * Filters out sites that have a QD annotation applied to them and where the QD value is lower than a
 * lower limit.
 *
 * @author Tim Fennell
 */
public class QdFilter implements VariantFilter {
    public static final String FILTER_NAME = "LowQD";
    private final double minimumQd;

    public QdFilter(final double minimumQd) {
        this.minimumQd = minimumQd;
    }

    @Override public String filter(final VariantContext ctx) {
        final double qd = ctx.getAttributeAsDouble("QD", -1d); // If QD is missing, return -1

        // QD should always be positive so a value < 0 indicates a missing value
        if (qd >= 0 && qd < minimumQd) {
            return FILTER_NAME;
        }
        else {
            return null;
        }
    }

    @Override
    public List<VCFFilterHeaderLine> headerLines() {
        return CollectionUtil.makeList(new VCFFilterHeaderLine(FILTER_NAME, "Site exhibits QD value below a hard limit."));
    }
}
