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

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Lazy-loading GenotypesContext.  A lazy-loading context has access to the
 * VCFParser and a unparsed string of genotype data.  If the user attempts to manipulate
 * the genotypes contained in this context, we decode the data and become a full blown
 * GenotypesContext.  However, if the user never does this we are spared a lot of expense
 * decoding the genotypes unnecessarily.
 */
public class LazyGenotypesContext extends GenotypesContext {
    /** The LazyParser we'll use to decode unparsedGenotypeData if necessary */
    final LazyParser parser;

    Object unparsedGenotypeData;

    /**
     * nUnparsedGenotypes the number of genotypes contained in the unparsedGenotypes data
     * (known already in the parser).  Useful for isEmpty and size() optimizations
     */
    final int nUnparsedGenotypes;

    /**
     * True if we've already decoded the values in unparsedGenotypeData
     */
    boolean loaded = false;

    private final static ArrayList<Genotype> EMPTY = new ArrayList<Genotype>(0);

    /**
     * Simple lazy parser interface.  Provide an object implementing this
     * interface to LazyGenotypesContext, and it's parse method will be called
     * when the use of the lazy context requires the underlying genotypes data
     * be parsed into Genotype objects.  The data argument is the data provided
     * to the LazyGenotypesContext holding encoded genotypes data
     */
    public interface LazyParser {
        @Requires("data != null")
        @Ensures("result != null")
        public LazyData parse(Object data);
    }

    /**
     * Returns the data used in the full GenotypesContext constructor
     *
     * {@link GenotypesContext#GenotypesContext(java.util.ArrayList, java.util.Map, java.util.List)}
     */
    public static class LazyData {
        final ArrayList<Genotype> genotypes;
        final Map<String, Integer> sampleNameToOffset;
        final List<String> sampleNamesInOrder;

        @Requires({"genotypes != null", "sampleNamesInOrder != null", "sampleNameToOffset != null"})
        public LazyData(final ArrayList<Genotype> genotypes,
                        final List<String> sampleNamesInOrder,
                        final Map<String, Integer> sampleNameToOffset) {
            this.genotypes = genotypes;
            this.sampleNamesInOrder = sampleNamesInOrder;
            this.sampleNameToOffset = sampleNameToOffset;
        }
    }

    /**
     * Creates a new lazy loading genotypes context using the LazyParser to create
     * genotypes data on demand.
     *
     * @param parser the parser to be used to load on-demand genotypes data
     * @param unparsedGenotypeData the encoded genotypes data that we will decode if necessary
     * @param nUnparsedGenotypes the number of genotypes that will be produced if / when we actually decode the genotypes data
     */
    @Requires({"parser != null", "unparsedGenotypeData != null", "nUnparsedGenotypes >= 0"})
    public LazyGenotypesContext(final LazyParser parser, final Object unparsedGenotypeData, final int nUnparsedGenotypes) {
        super(EMPTY);
        this.parser = parser;
        this.unparsedGenotypeData = unparsedGenotypeData;
        this.nUnparsedGenotypes = nUnparsedGenotypes;
    }

    /**
     * Overrides the genotypes accessor.  If we haven't already, decode the genotypes data
     * and store the decoded results in the appropriate variables.  Otherwise we just
     * returned the decoded result directly.  Note some care needs to be taken here as
     * the value in notToBeDirectlyAccessedGenotypes may diverge from what would be produced
     * by decode, if after the first decode the genotypes themselves are replaced
     * @return
     */
    @Override
    @Ensures("result != null")
    protected ArrayList<Genotype> getGenotypes() {
        decode();
        return notToBeDirectlyAccessedGenotypes;
    }

    /**
     * Force us to decode the genotypes, if not already done
     */
    public void decode() {
        if ( ! loaded ) {
            //System.out.printf("Loading genotypes... %s:%d%n", contig, start);
            LazyData parsed = parser.parse(unparsedGenotypeData);
            notToBeDirectlyAccessedGenotypes = parsed.genotypes;
            sampleNamesInOrder = parsed.sampleNamesInOrder;
            sampleNameToOffset = parsed.sampleNameToOffset;
            loaded = true;
            unparsedGenotypeData = null; // don't hold the unparsed data any longer

            // warning -- this path allows us to create a VariantContext that doesn't run validateGenotypes()
            // That said, it's not such an important routine -- it's just checking that the genotypes
            // are well formed w.r.t. the alleles list, but this will be enforced within the VCFCodec
        }
    }

    /**
     * Overrides the ensure* functionality.  If the data hasn't been loaded
     * yet and we want to build the cache, just decode it and we're done.  If we've
     * already decoded the data, though, go through the super class
     */
    @Override
    protected synchronized void ensureSampleNameMap() {
        if ( ! loaded ) {
            decode(); // will load up all of the necessary data
        } else {
            super.ensureSampleNameMap();
        }
    }

    @Override
    protected synchronized void ensureSampleOrdering() {
        if ( ! loaded ) {
            decode(); // will load up all of the necessary data
        } else {
            super.ensureSampleOrdering();
        }
    }

    @Override
    protected void invalidateSampleNameMap() {
        // if the cache is invalidated, and we haven't loaded our data yet, do so
        if ( ! loaded ) decode();
        super.invalidateSampleNameMap();
    }

    @Override
    protected void invalidateSampleOrdering() {
        // if the cache is invalidated, and we haven't loaded our data yet, do so
        if ( ! loaded ) decode();
        super.invalidateSampleOrdering();
    }

    @Override
    public boolean isEmpty() {
        // optimization -- we know the number of samples in the unparsed data, so use it here to
        // avoid parsing just to know if the genotypes context is empty
        return loaded ? super.isEmpty() : nUnparsedGenotypes == 0;
    }

    @Override
    public int size() {
        // optimization -- we know the number of samples in the unparsed data, so use it here to
        // avoid parsing just to know the size of the context
        return loaded ? super.size() : nUnparsedGenotypes;
    }

    public Object getUnparsedGenotypeData() {
        return unparsedGenotypeData;
    }
}
