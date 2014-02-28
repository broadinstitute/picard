package net.sf.samtools;

/**
 * The MIT License
 * <p/>
 * Copyright (c) 2014 The Broad Institute
 * <p/>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * <p/>
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * <p/>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

import java.util.Comparator;

/**
 * Provides ordering based on SAM header records' attribute values. Provide the list of attributes to use
 * in the comparison to the constructor. Null attribute values (i.e., those attributes not present in the
 * record) sort behind those that have values.
 */
public class SAMHeaderRecordComparator<T extends AbstractSAMHeaderRecord> implements Comparator<T> {

	private final String[] attributes;

	public SAMHeaderRecordComparator(final String... attributes) {
		this.attributes = attributes;
	}

	@Override
	public int compare(final T left, final T right) {
		for (final String attribute : attributes) {
			final String leftValue = left.getAttribute(attribute);
			final String rightValue = right.getAttribute(attribute);

			if (leftValue == null) {
				// Fastest comparison possible; two empty values are
				// equivalent, so move along to the next attribute
				if (rightValue == null) continue;

					// Otherwise left < right, since right has a value
				else return -1;
			}

			// left is not null; if right is, left > right
			if (rightValue == null) return 1;

			final int compare = leftValue.compareTo(rightValue);
			if (compare != 0) return compare;
		}

		return 0;
	}
}
