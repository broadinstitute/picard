package net.sf.picard.vcf;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A Comparator that orders VariantContexts by the ordering of the contigs/chromosomes in the List
 * provided at construction time, then by start position with each contig/chromosome.
 */
public class VariantContextComparator implements Comparator<VariantContext> {

	// For fast lookup of the contig's index in the contig list
	private final Map<String, Integer> contigIndexLookup;

	public VariantContextComparator(final List<String> contigs) {
		if (contigs.size() == 0) throw new IllegalArgumentException("One or more contigs must be in the contig list.");

		final Map<String, Integer> protoContigIndexLookup = new HashMap<String, Integer>();
		int index = 0;
		for (final String contig : contigs) {
			protoContigIndexLookup.put(contig, index++);
		}

		if (protoContigIndexLookup.size() != contigs.size()) {
			throw new IllegalArgumentException("There are duplicate contigs/chromosomes in the input contig list.");
		}

		this.contigIndexLookup = Collections.unmodifiableMap(protoContigIndexLookup);
	}

	/**
	 * Creates a VariantContextComparator from the given VCF contig header lines. The header lines'
	 * index values are used to order the contigs. Throws IllegalArgumentException if there are dupe
	 *
	 */
	public VariantContextComparator(final Collection<VCFContigHeaderLine> headerLines) {
		if (headerLines.size() == 0) throw new IllegalArgumentException("One or more header lines must be in the header line collection.");

		final Map<String, Integer> protoContigIndexLookup = new HashMap<String, Integer>();
		for (final VCFContigHeaderLine headerLine : headerLines) {
			protoContigIndexLookup.put(headerLine.getID(), headerLine.getContigIndex());
		}

		if (protoContigIndexLookup.size() != headerLines.size()) {
			throw new IllegalArgumentException("There are duplicate contigs/chromosomes in the input header line collection.");
		}

		final Set<Integer> protoIndexValues = new HashSet<Integer>(protoContigIndexLookup.values());
		if (protoIndexValues.size() != headerLines.size()) {
			throw new IllegalArgumentException("One or more contigs share the same index number.");
		}

		this.contigIndexLookup = Collections.unmodifiableMap(protoContigIndexLookup);
	}

	@Override
	public int compare(final VariantContext firstVariantContext, final VariantContext secondVariantContext) {
		// Will throw NullPointerException -- happily -- if either of the chromosomes/contigs aren't
		// present. This error checking should already have been done in the constructor but it's left
		// in as defence anyway.
		final int contigCompare =
				this.contigIndexLookup.get(firstVariantContext.getChr()) - this.contigIndexLookup.get(secondVariantContext.getChr());
		return contigCompare != 0
				? contigCompare
				: firstVariantContext.getStart() - secondVariantContext.getStart();
	}

	/**
	 * Returns true if the given header lines are from a file sorted according to this VariantContextComparator.
	 * For sorting to work properly, the contig in each header line must have the same index.
	 */
	public boolean isCompatible(final Collection<VCFContigHeaderLine> headerLines) {
		for (final VCFContigHeaderLine headerLine : headerLines) {
			final Integer existingIndex = this.contigIndexLookup.get(headerLine.getID());
			if (existingIndex == null || headerLine.getContigIndex() != existingIndex.intValue()) return false;
		}

		return true;
	}
}
