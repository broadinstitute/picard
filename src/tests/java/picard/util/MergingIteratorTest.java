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
package picard.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.MergingIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;

public class MergingIteratorTest {

	private static class QueueBackedIterator<T> implements CloseableIterator<T> {

		private final Iterator<T> backing;
		QueueBackedIterator(final Queue<T> queue) {
			this.backing = queue.iterator();
		}

		@Override
		public void close() {
			// no-op
		}

		@Override
		public boolean hasNext() {
			return backing.hasNext();
		}

		@Override
		public T next() {
			return backing.next();
		}

		@Override
		public void remove() {
			backing.remove();
		}
	}

	private static final Comparator<Integer> INTEGER_COMPARATOR = new Comparator<Integer>() {
		@Override
		public int compare(Integer integer, Integer integer2) {
			return integer - integer2;
		}
	};

	@Test
	public void testOrderingAndCompleteness() {
		final Queue<Integer> queueOne = new LinkedList<Integer>();
		queueOne.add(1);
		queueOne.add(3);
		queueOne.add(5);

		final Queue<Integer> queueTwo = new LinkedList<Integer>();
		queueTwo.add(2);
		queueTwo.add(4);
		queueTwo.add(6);

		final Queue<Integer> queueThree = new LinkedList<Integer>();
		queueThree.add(0);
		queueThree.add(1);

		final Collection<CloseableIterator<Integer>> iterators = new ArrayList<CloseableIterator<Integer>>(3);
		Collections.addAll(
				iterators,
				new QueueBackedIterator<Integer>(queueOne),
				new QueueBackedIterator<Integer>(queueTwo),
				new QueueBackedIterator<Integer>(queueThree));

		final MergingIterator<Integer> mergingIterator = new MergingIterator<Integer>(
				INTEGER_COMPARATOR,
				iterators);

		int count = 0;
		int last = -1;
		while (mergingIterator.hasNext()) {
			final Integer integer = mergingIterator.next();
			count++;
			if (integer == 1) Assert.assertTrue(integer >= last);
			else Assert.assertTrue(integer > last);
			last = integer;
		}

		Assert.assertEquals(queueOne.size() + queueTwo.size() + queueThree.size(), count);
	}

	@Test
	public void testIteratorsOfUnevenLength() {
		final Queue<Integer> queueOne = new LinkedList<Integer>();
		queueOne.add(1);
		queueOne.add(3);
		queueOne.add(5);
		queueOne.add(7);
		queueOne.add(9);
		queueOne.add(11);
		queueOne.add(13);

		final Queue<Integer> queueTwo = new LinkedList<Integer>();
		queueTwo.add(2);

		final Collection<CloseableIterator<Integer>> iterators = new ArrayList<CloseableIterator<Integer>>(3);
		Collections.addAll(
				iterators,
				new QueueBackedIterator<Integer>(queueOne),
				new QueueBackedIterator<Integer>(queueTwo));

		final MergingIterator<Integer> mergingIterator = new MergingIterator<Integer>(
				INTEGER_COMPARATOR,
				iterators);

		int count = 0;
		int last = -1;
		while (mergingIterator.hasNext()) {
			final Integer integer = mergingIterator.next();
			count++;
			Assert.assertTrue(integer > last);
			last = integer;
		}

		Assert.assertEquals(queueOne.size() + queueTwo.size(), count);
	}

	@Test(expectedExceptions = IllegalStateException.class)
	public void testOutOfOrderIterators() {
		final Queue<Integer> queueOne = new LinkedList<Integer>();
		queueOne.add(1);
		queueOne.add(3);

		final Queue<Integer> queueTwo = new LinkedList<Integer>();
		queueTwo.add(4);
		queueTwo.add(2);

		final Collection<CloseableIterator<Integer>> iterators = new ArrayList<CloseableIterator<Integer>>(3);
		Collections.addAll(
				iterators,
				new QueueBackedIterator<Integer>(queueOne),
				new QueueBackedIterator<Integer>(queueTwo));

		final MergingIterator<Integer> mergingIterator = new MergingIterator<Integer>(
				INTEGER_COMPARATOR,
				iterators);

		Assert.assertEquals(mergingIterator.next().intValue(), 1);
		Assert.assertEquals(mergingIterator.next().intValue(), 3);
		Assert.assertEquals(mergingIterator.next().intValue(), 4);
		mergingIterator.next(); // fails, because the next element would be "2"
	}
}
