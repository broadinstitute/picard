/*
 * The MIT License
 *
 * Copyright (c) 2015-2016 The Broad Institute
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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;

import java.util.Iterator;
import java.util.Optional;

/**
 * An iterator that takes a pair of iterators over VariantContexts and iterates over them in tandem.
 *
 * A tuple will be returned with variant contexts for both contexts if present.  Otherwise, the missing
 * context at that site will be empty.  The contexts will be returned in coordinate order.
 *
 * */
public class PairedVariantSubContextIterator implements Iterator<PairedVariantSubContextIterator.VcfTuple> {
    private final PeekableIterator<VariantContext> leftIterator;
    private final String leftSample;
    private final PeekableIterator<VariantContext> rightIterator;
    private final String rightSample;
    private final VariantContextComparator comparator;

    public PairedVariantSubContextIterator(final Iterator<VariantContext> leftIterator, final String leftSample,
                                    final Iterator<VariantContext> rightIterator, final String rightSample,
                                    final SAMSequenceDictionary dict) {
        this.leftIterator  = new PeekableIterator<>(leftIterator);
        this.leftSample    = leftSample;
        this.rightIterator = new PeekableIterator<>(rightIterator);
        this.rightSample   = rightSample;
        this.comparator    = new VariantContextComparator(dict);
    }

    @Override
    public boolean hasNext() {
        return this.leftIterator.hasNext() || this.rightIterator.hasNext();
    }

    @Override
    public VcfTuple next() {
        if (!hasNext()) throw new IllegalStateException("next() called while hasNext() is false.");

        final Optional<VariantContext> leftVariantContext  = this.leftIterator.hasNext() ? Optional.of(this.leftIterator.peek()) : Optional.empty();
        final Optional<VariantContext> rightVariantContext = this.rightIterator.hasNext() ? Optional.of(this.rightIterator.peek()) : Optional.empty();

        // If one or the other is missing because there is no next, just return a one-sided tuple
        if (!leftVariantContext.isPresent() && !rightVariantContext.isPresent()) {
            throw new IllegalStateException("BUG: Both contexts empty.");
        }
        else if (!leftVariantContext.isPresent()) {
            return new VcfTuple(Optional.empty(), this.rightIterator.next().subContextFromSample(rightSample));
        }
        else if (!rightVariantContext.isPresent()) {
            return new VcfTuple(this.leftIterator.next().subContextFromSample(leftSample), Optional.empty());
        }
        else { // Otherwise check the ordering and do the right thing
            final int ordering = this.comparator.compare(leftVariantContext.get(), rightVariantContext.get());
            if (ordering == 0) {
                return new VcfTuple(this.leftIterator.next().subContextFromSample(leftSample), this.rightIterator.next().subContextFromSample(rightSample));
            } else if (ordering < 0) {
                return new VcfTuple(this.leftIterator.next().subContextFromSample(leftSample), Optional.empty());
            } else {
                return new VcfTuple(Optional.empty(), this.rightIterator.next().subContextFromSample(rightSample));
            }
        }
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /** Little class to hold a pair of VariantContexts that are in sync with one another. */
    public static class VcfTuple {
        public final Optional<VariantContext> leftVariantContext;
        public final Optional<VariantContext> rightVariantContext;

        private VcfTuple(final Optional<VariantContext> leftVariantContext, final Optional<VariantContext> rightVariantContext) {
            this.leftVariantContext  = leftVariantContext;
            this.rightVariantContext = rightVariantContext;
        }

        VcfTuple(final VariantContext leftVariantContext, final VariantContext rightVariantContext) {
            this(Optional.of(leftVariantContext), Optional.of(rightVariantContext));
        }

        VcfTuple(final Optional<VariantContext> leftVariantContext, final VariantContext rightVariantContext) {
            this(leftVariantContext, Optional.of(rightVariantContext));
        }

        VcfTuple(final VariantContext leftVariantContext, final Optional<VariantContext> rightVariantContext) {
            this(Optional.of(leftVariantContext), rightVariantContext);
        }
    }
}
