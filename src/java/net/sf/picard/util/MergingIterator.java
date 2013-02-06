/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
package net.sf.picard.util;

import net.sf.samtools.util.CloseableIterator;

import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.PriorityQueue;

/**
 * An iterator over Iterators that return Ts. Calling next() returns the next T ordered according
 * to Comparator provided at construction time. Importantly, the elements in the input Iterators
 * must already be sorted according to the provided Comparator.
 */
public class MergingIterator<T> implements CloseableIterator<T> {

	/*
	 * An Iterator whose natural ordering is by the T that will be returned by the next call to
	 * next().
	 */
	private class ComparableIterator extends PeekableIterator<T> implements Comparable<ComparableIterator> {

		public ComparableIterator(final Iterator<T> iterator) {
			super(iterator);
		}

		@Override
		public int compareTo(final ComparableIterator that) {
			if (comparator.getClass() != comparator.getClass()) {
				throw new IllegalStateException("Can't compare two ComparableIterators that have different orderings.");
			}

			return comparator.compare(this.peek(), that.peek());
		}
	}

	/*
	 * The general flow is to pull the "top" (according to the ComparableIterator's compareTo())
	 * iterator off on calls to this.next(), get iterator.next() and then re-add the iterator to
	 * the queue. Readding reorders the queue so the next "top" iterator is ready.
	 */
	private final PriorityQueue<ComparableIterator> queue;

	private final Comparator<T> comparator;

	// This is the last T returned by the call to next(). It's used to make sure that the comparators
	// always return correctly ordered Ts.
	private T lastReturned;

	/**
	 * Creates a MergingIterator over the given Collection of iterators whose elements will be
	 * returned in the order defined by the given Comparator.
	 */
	public MergingIterator(final Comparator<T> comparator, final Collection<CloseableIterator<T>> iterators) {
		if (iterators.isEmpty()) throw new IllegalArgumentException("One or more CloseableIterators must be provided.");

		this.comparator = comparator;

		this.queue = new PriorityQueue<ComparableIterator>();
		for (final CloseableIterator<T> iterator : iterators) {
			this.addIfNotEmpty(new ComparableIterator(iterator));
		}

		// If there are no iterators to read from after adding them all to the prioqueue,
		// should we throw? it's prob'ly an error.
	}

	/**
	 * @see java.util.Iterator<T>.hasNext
	 */
	@Override
	public boolean hasNext() {
		return ! this.queue.isEmpty();
	}

	/**
	 * @see java.util.Iterator<T>.next
	 */
	@Override
	public T next() {
		if ( ! this.hasNext()) throw new NoSuchElementException();

		final ComparableIterator recordIterator = this.queue.poll();
		// Assumes the iterator is closed & removed from the queue before recordIterator.hasNext() == false
		final T next = recordIterator.next();
		// I don't like having to test for null here -- it's really only null before the first call
		// to next() -- but I don't see any other way
		if (this.lastReturned != null && this.comparator.compare(lastReturned, next) > 0) {
			throw new IllegalStateException(
					"The elements of the input Iterators are not sorted according to the comparator " +
							this.comparator.getClass().getName());
		}

		addIfNotEmpty(recordIterator);
		this.lastReturned = next;
		return next;
	}

	/**
	 * Unsupported.
	 */
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	/**
	 * Closes every CloseableIterator in this MergingIterator. After calling, calls to
	 * hasNext() will always return false.
	 */
	@Override
	public void close() {
		for (final ComparableIterator iterator : this.queue) {
			iterator.close();
			this.queue.remove(iterator);
		}
	}

	private void addIfNotEmpty(final ComparableIterator iterator) {
		if (iterator.hasNext()) queue.offer(iterator);
		else iterator.close();
	}
}
