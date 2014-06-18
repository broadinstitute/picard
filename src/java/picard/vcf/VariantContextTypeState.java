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
package picard.vcf;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * A small class to store the state between two VariantContext.Types
 *
 * @author George Grant
 */

public class VariantContextTypeState implements Comparable<VariantContextTypeState> {
    final VariantContext.Type truthType, callType;

    public VariantContextTypeState(final VariantContext.Type truthType, final VariantContext.Type callType) {
        this.truthType = truthType;
        this.callType = callType;
    }

    @Override public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        return compareTo((VariantContextTypeState) o) == 0;
    }

    @Override public int hashCode() {
        int result = truthType.hashCode();
        result = 31 * result + callType.hashCode();
        return result;
    }

    @Override public int compareTo(final VariantContextTypeState that) {
        int result = this.truthType.compareTo(that.truthType);
        if (result == 0) result = this.callType.compareTo(that.callType);
        return result;
    }

}
