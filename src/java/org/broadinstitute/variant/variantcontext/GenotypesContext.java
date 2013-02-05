/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.variantcontext;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.*;

/**
 * Represents an ordered collection of Genotype objects
 */
public class GenotypesContext implements List<Genotype> {
    /**
     * static constant value for an empty GenotypesContext.  Useful since so many VariantContexts have no genotypes
     */
    public final static GenotypesContext NO_GENOTYPES =
            new GenotypesContext(new ArrayList<Genotype>(0), new HashMap<String, Integer>(0), Collections.<String>emptyList()).immutable();

    /**
     *sampleNamesInOrder a list of sample names, one for each genotype in genotypes, sorted in alphabetical order
     */
    List<String> sampleNamesInOrder = null;

    /**
     * a map optimized for efficient lookup.  Each genotype in genotypes must have its
     * sample name in sampleNameToOffset, with a corresponding integer value that indicates the offset of that
     * genotype in the vector of genotypes
     */
    Map<String, Integer> sampleNameToOffset = null;

    /**
     * An ArrayList of genotypes contained in this context
     *
     * WARNING: TO ENABLE THE LAZY VERSION OF THIS CLASS, NO METHODS SHOULD DIRECTLY
     * ACCESS THIS VARIABLE.  USE getGenotypes() INSTEAD.
     *
     */
    ArrayList<Genotype> notToBeDirectlyAccessedGenotypes;

    /**
     * Cached value of the maximum ploidy observed among all samples
     */
    private int maxPloidy = -1;

    /** Are we allowing users to modify the list? */
    boolean immutable = false;

    // ---------------------------------------------------------------------------
    //
    // private constructors -- you have to use static create methods to make these classes
    //
    // ---------------------------------------------------------------------------

    /**
     * Create an empty GenotypeContext
     */
    protected GenotypesContext() {
        this(10);
    }

    /**
     * Create an empty GenotypeContext, with initial capacity for n elements
     */
    @Requires("n >= 0")
    protected GenotypesContext(final int n) {
        this(new ArrayList<Genotype>(n));
    }

    /**
     * Create an GenotypeContext containing genotypes
     */
    @Requires("genotypes != null")
    protected GenotypesContext(final ArrayList<Genotype> genotypes) {
        this.notToBeDirectlyAccessedGenotypes = genotypes;
        this.sampleNameToOffset = null;
    }

    /**
     * Create a fully resolved GenotypeContext containing genotypes, sample lookup table,
     * and sorted sample names
     *
     * @param genotypes our genotypes in arbitrary
     * @param sampleNameToOffset map optimized for efficient lookup.  Each genotype in genotypes must have its
     * sample name in sampleNameToOffset, with a corresponding integer value that indicates the offset of that
     * genotype in the vector of genotypes
     * @param sampleNamesInOrder a list of sample names, one for each genotype in genotypes, sorted in alphabetical
     * order.
     */
    @Requires({"genotypes != null",
            "sampleNameToOffset != null",
            "sampleNamesInOrder != null",
            "genotypes.size() == sampleNameToOffset.size()",
            "genotypes.size() == sampleNamesInOrder.size()"})
    protected GenotypesContext(final ArrayList<Genotype> genotypes,
                             final Map<String, Integer> sampleNameToOffset,
                             final List<String> sampleNamesInOrder) {
        this.notToBeDirectlyAccessedGenotypes = genotypes;
        this.sampleNameToOffset = sampleNameToOffset;
        this.sampleNamesInOrder = sampleNamesInOrder;
    }

    // ---------------------------------------------------------------------------
    //
    // public static factory methods
    //
    // ---------------------------------------------------------------------------

    /**
     * Basic creation routine
     * @return an empty, mutable GenotypeContext
     */
    @Ensures({"result != null"})
    public static final GenotypesContext create() {
        return new GenotypesContext();
    }

    /**
     * Basic creation routine
     * @return an empty, mutable GenotypeContext with initial capacity for nGenotypes
     */
    @Requires("nGenotypes >= 0")
    @Ensures({"result != null"})
    public static final GenotypesContext create(final int nGenotypes) {
        return new GenotypesContext(nGenotypes);
    }

    /**
     * Create a fully resolved GenotypeContext containing genotypes, sample lookup table,
     * and sorted sample names
     *
     * @param genotypes our genotypes in arbitrary
     * @param sampleNameToOffset map optimized for efficient lookup.  Each genotype in genotypes must have its
     * sample name in sampleNameToOffset, with a corresponding integer value that indicates the offset of that
     * genotype in the vector of genotypes
     * @param sampleNamesInOrder a list of sample names, one for each genotype in genotypes, sorted in alphabetical
     * order.
     * @return an mutable GenotypeContext containing genotypes with already present lookup data
     */
    @Requires({"genotypes != null",
            "sampleNameToOffset != null",
            "sampleNamesInOrder != null"})
    @Ensures({"result != null"})
    public static final GenotypesContext create(final ArrayList<Genotype> genotypes,
                                                final Map<String, Integer> sampleNameToOffset,
                                                final List<String> sampleNamesInOrder) {
        return new GenotypesContext(genotypes, sampleNameToOffset, sampleNamesInOrder);
    }

    /**
     * Create a fully resolved GenotypeContext containing genotypes
     *
     * @param genotypes our genotypes in arbitrary
     * @return an mutable GenotypeContext containing genotypes
     */
    @Requires({"genotypes != null"})
    @Ensures({"result != null"})
    public static final GenotypesContext create(final ArrayList<Genotype> genotypes) {
        return genotypes == null ? NO_GENOTYPES : new GenotypesContext(genotypes);
    }

    /**
     * Create a fully resolved GenotypeContext containing genotypes
     *
     * @param genotypes our genotypes in arbitrary
     * @return an mutable GenotypeContext containing genotypes
     */
    @Requires({"genotypes != null"})
    @Ensures({"result != null"})
    public static final GenotypesContext create(final Genotype... genotypes) {
        return create(new ArrayList<Genotype>(Arrays.asList(genotypes)));
    }

    /**
     * Create a freshly allocated GenotypeContext containing the genotypes in toCopy
     *
     * @param toCopy the GenotypesContext to copy
     * @return an mutable GenotypeContext containing genotypes
     */
    @Requires({"toCopy != null"})
    @Ensures({"result != null"})
    public static final GenotypesContext copy(final GenotypesContext toCopy) {
        return create(new ArrayList<Genotype>(toCopy.getGenotypes()));
    }

    /**
     * Create a GenotypesContext containing the genotypes in iteration order contained
     * in toCopy
     *
     * @param toCopy the collection of genotypes
     * @return an mutable GenotypeContext containing genotypes
     */
    @Ensures({"result != null"})
    public static final GenotypesContext copy(final Collection<Genotype> toCopy) {
        return toCopy == null ? NO_GENOTYPES : create(new ArrayList<Genotype>(toCopy));
    }

    // ---------------------------------------------------------------------------
    //
    // Mutability methods
    //
    // ---------------------------------------------------------------------------

    public final GenotypesContext immutable() {
        immutable = true;
        return this;
    }

    public boolean isMutable() {
        return ! immutable;
    }

    public final void checkImmutability() {
        if ( immutable )
            throw new IllegalAccessError("GenotypeMap is currently immutable, but a mutator method was invoked on it");
    }

    // ---------------------------------------------------------------------------
    //
    // caches
    //
    // ---------------------------------------------------------------------------

    @Ensures({"sampleNameToOffset == null"})
    protected void invalidateSampleNameMap() {
        sampleNameToOffset = null;
    }

    @Ensures({"sampleNamesInOrder == null"})
    protected void invalidateSampleOrdering() {
        sampleNamesInOrder = null;
    }

    @Ensures({"sampleNamesInOrder != null"})
    protected void ensureSampleOrdering() {
        if ( sampleNamesInOrder == null ) {
            sampleNamesInOrder = new ArrayList<String>(size());

            for ( int i = 0; i < size(); i++ ) {
                sampleNamesInOrder.add(getGenotypes().get(i).getSampleName());
            }
            Collections.sort(sampleNamesInOrder);
        }
    }

    @Ensures({"sampleNameToOffset != null"})
    protected void ensureSampleNameMap() {
        if ( sampleNameToOffset == null ) {
            sampleNameToOffset = new HashMap<String, Integer>(size());

            for ( int i = 0; i < size(); i++ ) {
                sampleNameToOffset.put(getGenotypes().get(i).getSampleName(), i);
            }
        }
    }

    // ---------------------------------------------------------------------------
    //
    // Lazy methods
    //
    // ---------------------------------------------------------------------------

    public boolean isLazyWithData() {
        return this instanceof LazyGenotypesContext &&
                ((LazyGenotypesContext)this).getUnparsedGenotypeData() != null;
    }

    // ---------------------------------------------------------------------------
    //
    // Map methods
    //
    // ---------------------------------------------------------------------------

    protected ArrayList<Genotype> getGenotypes() {
        return notToBeDirectlyAccessedGenotypes;
    }

    @Override
    public void clear() {
        checkImmutability();
        invalidateSampleNameMap();
        invalidateSampleOrdering();
        getGenotypes().clear();
    }

    @Override
    public int size() {
        return getGenotypes().size();
    }

    @Override
    public boolean isEmpty() {
        return getGenotypes().isEmpty();
    }

    /**
     * Adds a single genotype to this context.
     *
     * There are many constraints on this input, and important
     * impacts on the performance of other functions provided by this
     * context.
     *
     * First, the sample name of genotype must be unique within this
     * context.  However, this is not enforced in the code itself, through
     * you will invalid the contract on this context if you add duplicate
     * samples and are running with CoFoJa enabled.
     *
     * Second, adding genotype also updates the sample name -> index map,
     * so add() followed by containsSample and related function is an efficient
     * series of operations.
     *
     * Third, adding the genotype invalidates the sorted list of sample names, to
     * add() followed by any of the SampleNamesInOrder operations is inefficient, as
     * each SampleNamesInOrder must rebuild the sorted list of sample names at
     * an O(n log n) cost.
     *
     * @param genotype
     * @return
     */
    @Override
    @Requires({"genotype != null", "get(genotype.getSampleName()) == null"})
    public boolean add(final Genotype genotype) {
        checkImmutability();
        invalidateSampleOrdering();

        if ( sampleNameToOffset != null ) {
            // update the name map by adding entries
            sampleNameToOffset.put(genotype.getSampleName(), size());
        }

        return getGenotypes().add(genotype);
    }

    @Override
    @Requires("! contains(genotype)")
    public void add(final int i, final Genotype genotype) {
        throw new UnsupportedOperationException();
    }

    /**
     * Adds all of the genotypes to this context
     *
     * See {@link #add(Genotype)} for important information about this functions
     * constraints and performance costs
     *
     * @param genotypes
     * @return
     */
    @Override
    @Requires("! containsAny(genotypes)")
    public boolean addAll(final Collection<? extends Genotype> genotypes) {
        checkImmutability();
        invalidateSampleOrdering();

        if ( sampleNameToOffset != null ) {
            // update the name map by adding entries
            int pos = size();
            for ( final Genotype g : genotypes ) {
                sampleNameToOffset.put(g.getSampleName(), pos++);
            }
        }

        return getGenotypes().addAll(genotypes);
    }

    @Override
    public boolean addAll(final int i, final Collection<? extends Genotype> genotypes) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean contains(final Object o) {
        return getGenotypes().contains(o);
    }

    @Override
    public boolean containsAll(final Collection<?> objects) {
        return getGenotypes().containsAll(objects);
    }

    private boolean containsAny(final Collection<? extends Genotype> genotypes) {
        for ( final Genotype g : genotypes ) {
            if ( contains(g) ) return true;
        }
        return false;
    }

    @Override
    public Genotype get(final int i) {
        return getGenotypes().get(i);
    }

    /**
     * What is the max ploidy among all samples?  Returns defaultPloidy if no genotypes are present
     *
     * @param defaultPloidy the default ploidy, if all samples are no-called
     * @return
     */
    @Ensures("result >= 0")
    public int getMaxPloidy(final int defaultPloidy) {
        if ( defaultPloidy < 0 ) throw new IllegalArgumentException("defaultPloidy must be greater than or equal to 0");

        if ( maxPloidy == -1 ) {
            maxPloidy = 0; // necessary in the case where there are no genotypes
            for ( final Genotype g : getGenotypes() ) {
                maxPloidy = Math.max(g.getPloidy(), maxPloidy);
            }

            // everything is no called so we return the default ploidy
            if ( maxPloidy == 0 ) maxPloidy = defaultPloidy;
        }

        return maxPloidy;
    }

    /**
     * Gets sample associated with this sampleName, or null if none is found
     *
     * @param sampleName
     * @return
     */
    public Genotype get(final String sampleName) {
        Integer offset = getSampleI(sampleName);
        return offset == null ? null : getGenotypes().get(offset);
    }

    private Integer getSampleI(final String sampleName) {
        ensureSampleNameMap();
        return sampleNameToOffset.get(sampleName);
    }

    @Override
    public int indexOf(final Object o) {
        return getGenotypes().indexOf(o);
    }

    @Override
    public Iterator<Genotype> iterator() {
        return getGenotypes().iterator();
    }

    @Override
    public int lastIndexOf(final Object o) {
        return getGenotypes().lastIndexOf(o);
    }

    @Override
    public ListIterator<Genotype> listIterator() {
        // todo -- must be immutable
        throw new UnsupportedOperationException();
//        return genotypes.listIterator();
    }

    @Override
    public ListIterator<Genotype> listIterator(final int i) {
        // todo -- must be immutable
        throw new UnsupportedOperationException();
//        return genotypes.listIterator(i);
    }

    /**
     * Note that remove requires us to invalidate our sample -> index
     * cache.  The loop:
     *
     * GenotypesContext gc = ...
     * for ( sample in samples )
     *   if ( gc.containsSample(sample) )
     *     gc.remove(sample)
     *
     * is extremely inefficient, as each call to remove invalidates the cache
     * and containsSample requires us to rebuild it, an O(n) operation.
     *
     * If you must remove many samples from the GC, use either removeAll or retainAll
     * to avoid this O(n * m) operation.
     *
     * @param i
     * @return
     */
    @Override
    public Genotype remove(final int i) {
        checkImmutability();
        invalidateSampleNameMap();
        invalidateSampleOrdering();
        return getGenotypes().remove(i);
    }

    /**
     * See for important warning {@link this.remove(Integer)}
     * @param o
     * @return
     */
    @Override
    public boolean remove(final Object o) {
        checkImmutability();
        invalidateSampleNameMap();
        invalidateSampleOrdering();
        return getGenotypes().remove(o);
    }

    @Override
    public boolean removeAll(final Collection<?> objects) {
        checkImmutability();
        invalidateSampleNameMap();
        invalidateSampleOrdering();
        return getGenotypes().removeAll(objects);
    }

    @Override
    public boolean retainAll(final Collection<?> objects) {
        checkImmutability();
        invalidateSampleNameMap();
        invalidateSampleOrdering();
        return getGenotypes().retainAll(objects);
    }

    @Override
    public Genotype set(final int i, final Genotype genotype) {
        checkImmutability();
        final Genotype prev = getGenotypes().set(i, genotype);

        invalidateSampleOrdering();
        if ( sampleNameToOffset != null ) {
            // update the name map by removing the old entry and replacing it with the new one
            sampleNameToOffset.remove(prev.getSampleName());
            sampleNameToOffset.put(genotype.getSampleName(), i);
        }

        return prev;
    }

    /**
     * Replaces the genotype in this context -- note for efficiency
     * reasons we do not add the genotype if it's not present.  The
     * return value will be null indicating this happened.
     *
     * Note this operation is preserves the map cache Sample -> Offset but
     * invalidates the sorted list of samples.  Using replace within a loop
     * containing any of the SampleNameInOrder operation requires an O(n log n)
     * resorting after each replace operation.
     *
     * @param genotype a non null genotype to bind in this context
     * @return null if genotype was not added, otherwise returns the previous genotype
     */
    @Requires("genotype != null")
    public Genotype replace(final Genotype genotype) {
        checkImmutability();
        Integer offset = getSampleI(genotype.getSampleName());
        if ( offset == null )
            return null;
        else
            return set(offset, genotype);
    }

    @Override
    public List<Genotype> subList(final int i, final int i1) {
        return getGenotypes().subList(i, i1);
    }

    @Override
    public Object[] toArray() {
        return getGenotypes().toArray();
    }

    @Override
    public <T> T[] toArray(final T[] ts) {
        return getGenotypes().toArray(ts);
    }

    /**
     * Iterate over the Genotypes in this context in the order specified by sampleNamesInOrder
     *
     * @param sampleNamesInOrder a Iterable of String, containing exactly one entry for each Genotype sample name in
     * this context
     * @return a Iterable over the genotypes in this context.
     */
    @Requires("sampleNamesInOrder != null")
    public Iterable<Genotype> iterateInSampleNameOrder(final Iterable<String> sampleNamesInOrder) {
        return new Iterable<Genotype>() {
            @Override
            public Iterator<Genotype> iterator() {
                return new InOrderIterator(sampleNamesInOrder.iterator());
            }
        };
    }

    /**
     * Iterate over the Genotypes in this context in their sample name order (A, B, C)
     * regardless of the underlying order in the vector of genotypes
     * @return a Iterable over the genotypes in this context.
     */
    public Iterable<Genotype> iterateInSampleNameOrder() {
        return iterateInSampleNameOrder(getSampleNamesOrderedByName());
    }

    private final class InOrderIterator implements Iterator<Genotype> {
        final Iterator<String> sampleNamesInOrder;

        private InOrderIterator(final Iterator<String> sampleNamesInOrder) {
            this.sampleNamesInOrder = sampleNamesInOrder;
        }

        @Override
        public boolean hasNext() {
            return sampleNamesInOrder.hasNext();
        }

        @Override
        public Genotype next() {
            return get(sampleNamesInOrder.next());
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    /**
     * @return The set of sample names for all genotypes in this context, in arbitrary order
     */
    @Ensures("result != null")
    public Set<String> getSampleNames() {
        ensureSampleNameMap();
        return sampleNameToOffset.keySet();
    }

    /**
     * @return The set of sample names for all genotypes in this context, in their natural ordering (A, B, C)
     */
    @Ensures("result != null")
    public List<String> getSampleNamesOrderedByName() {
        ensureSampleOrdering();
        return sampleNamesInOrder;
    }

    @Requires("sample != null")
    public boolean containsSample(final String sample) {
        ensureSampleNameMap();
        return sampleNameToOffset.containsKey(sample);
    }

    @Requires("samples != null")
    public boolean containsSamples(final Collection<String> samples) {
        return getSampleNames().containsAll(samples);
    }

    /**
     * Return a freshly allocated subcontext of this context containing only the samples
     * listed in samples.  Note that samples can contain names not in this context, they
     * will just be ignored.
     *
     * @param samples
     * @return
     */
    @Requires("samples != null")
    @Ensures("result != null")
    public GenotypesContext subsetToSamples( final Set<String> samples ) {
        final int nSamples = samples.size();

        if ( nSamples == 0 )
            return NO_GENOTYPES;
        else { // nGenotypes < nSamples
            final GenotypesContext subset = create(samples.size());
            for ( final String sample : samples ) {
                final Genotype g = get(sample);
                if ( g != null )
                    subset.add(g);
            }
            return subset;
        }
    }

    @Override
    public String toString() {
        final List<String> gS = new ArrayList<String>();
        for ( final Genotype g : this.iterateInSampleNameOrder() )
            gS.add(g.toString());
        return "[" + join(",", gS) + "]";
    }

    // copied from Utils
    private static <T> String join(final String separator, final Collection<T> objects) {
        if (objects.isEmpty()) { // fast path for empty collection
            return "";
        } else {
            final Iterator<T> iter = objects.iterator();
            final T first = iter.next();

            if ( ! iter.hasNext() ) // fast path for singleton collections
                return first.toString();
            else { // full path for 2+ collection that actually need a join
                final StringBuilder ret = new StringBuilder(first.toString());
                while(iter.hasNext()) {
                    ret.append(separator);
                    ret.append(iter.next().toString());
                }
                return ret.toString();
            }
        }
    }
}
