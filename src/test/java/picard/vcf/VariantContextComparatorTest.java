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
package picard.vcf;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;

public class VariantContextComparatorTest {

	private static VariantContext buildVariantContext(final String source, final String contig, final long start) {
		final Collection<Allele> alleles = new ArrayList<Allele>();
		alleles.add(Allele.create("AAAA", true));
		alleles.add(Allele.create("AAGG", false));
		return new VariantContextBuilder(source, contig, start, start + 3, alleles).make();
	}

	private static List<String> getOrderedContigList(final VariantContext... variantContexts) {
		final LinkedHashSet<String> contigs = new LinkedHashSet<String>();
		for (final VariantContext context : variantContexts) {
			contigs.add(context.getChr());
		}
		return new ArrayList<String>(contigs);
	}

	@Test
	public void testIdentical() {
		final VariantContext contextOne = buildVariantContext("source", "one", 100);
		final List<String> contigs = getOrderedContigList(contextOne);
		Assert.assertEquals(0, new VariantContextComparator(contigs).compare(contextOne, contextOne));
	}

	@Test
	public void testPositions() {
		final VariantContext contextOne = buildVariantContext("source", "one", 100);
		final VariantContext contextTwo = buildVariantContext("source", "one", 150);
		final List<String> contigs = getOrderedContigList(contextOne, contextTwo);
		Assert.assertTrue(new VariantContextComparator(contigs).compare(contextOne, contextTwo) < 0);
		Assert.assertTrue(new VariantContextComparator(contigs).compare(contextTwo, contextOne) > 0);
	}

	@Test
	public void testContigs() {
		final VariantContext contextOne = buildVariantContext("source", "one", 100);
		final VariantContext contextTwo = buildVariantContext("source", "two", 100);
		final List<String> contigs = getOrderedContigList(contextOne, contextTwo);
		Assert.assertTrue(new VariantContextComparator(contigs).compare(contextOne, contextTwo) < 0);
		Assert.assertTrue(new VariantContextComparator(contigs).compare(contextTwo, contextOne) > 0);
	}

	@Test
	public void testCombinationOne() {
		final VariantContext contextOne = buildVariantContext("source", "one", 100);
		final VariantContext contextTwo = buildVariantContext("source", "two", 150);
		final List<String> contigs = getOrderedContigList(contextOne, contextTwo);
		Assert.assertTrue(new VariantContextComparator(contigs).compare(contextOne, contextTwo) < 0);
		Assert.assertTrue(new VariantContextComparator(contigs).compare(contextTwo, contextOne) > 0);
	}

	@Test
	public void testCombinationTwo() {
		final VariantContext contextOne = buildVariantContext("source", "one", 150);
		final VariantContext contextTwo = buildVariantContext("source", "two", 100);
		final List<String> contigs = getOrderedContigList(contextOne, contextTwo);
		Assert.assertTrue(new VariantContextComparator(contigs).compare(contextOne, contextTwo) < 0);
		Assert.assertTrue(new VariantContextComparator(contigs).compare(contextTwo, contextOne) > 0);
	}

	@Test (expectedExceptions = IllegalArgumentException.class)
	public void testThrowsOnDuplicateContig() {
		final List<String> contigs = new ArrayList<String>(3);
		contigs.add("one");
		contigs.add("two");
		contigs.add("one");

		new VariantContextComparator(contigs);
	}
}
