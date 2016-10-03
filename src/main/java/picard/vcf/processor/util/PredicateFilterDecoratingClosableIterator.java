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

package picard.vcf.processor.util;

import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.collect.Iterators;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Collection;
import java.util.Iterator;

/**
 * Performs on-the-fly filtering of the provided {@link VariantContext} {@link java.util.Iterator} such that only variants that satisfy
 * all predicates are emitted.
 *
 * This class only exists because {@link Iterators#filter(Iterator, Predicate)} won't produce a {@link CloseableIterator}, which is 
 * necessary. 
 * 
 * @author mccowan
 */
public class PredicateFilterDecoratingClosableIterator<T> implements CloseableIterator<T> {
    final CloseableIterator<T> underlyingIterator;
    final Iterator<T> filteredIterator;

    public PredicateFilterDecoratingClosableIterator(final CloseableIterator<T> underlyingIterator, final Collection<Predicate<T>> predicates) {
        Preconditions.checkArgument(!predicates.isEmpty(), "predicates must not be empty");
        Iterator<T> nestedPredicateIterator = underlyingIterator;
        for (final Predicate<T> predicate : predicates) {
           nestedPredicateIterator = Iterators.filter(nestedPredicateIterator, predicate);   
        }
        filteredIterator = nestedPredicateIterator;
        
        this.underlyingIterator = underlyingIterator;
    }
    
    @Override
    public boolean hasNext() {
        return filteredIterator.hasNext();
    }

    @Override
    public T next() {
        return filteredIterator.next();
    }

    @Override
    public void close() {
        underlyingIterator.close();
    }

    @Override
    public void remove() {
        underlyingIterator.remove();
    }
}
