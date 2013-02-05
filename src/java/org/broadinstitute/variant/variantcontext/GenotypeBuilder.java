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
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.vcf.VCFConstants;

import java.util.*;

/**
 * A builder class for genotypes
 *
 * Provides convenience setter methods for all of the Genotype field
 * values.  Setter methods can be used in any order, allowing you to
 * pass through states that wouldn't be allowed in the highly regulated
 * immutable Genotype class.
 *
 * All fields default to meaningful MISSING values.
 *
 * Call make() to actually create the corresponding Genotype object from
 * this builder.  Can be called multiple times to create independent copies,
 * or with intervening sets to conveniently make similar Genotypes with
 * slight modifications.
 *
 * @author Mark DePristo
 * @since 06/12
 */
@Invariant({"alleles != null"})
public final class GenotypeBuilder {
    private static final List<Allele> HAPLOID_NO_CALL = Arrays.asList(Allele.NO_CALL);
    private static final List<Allele> DIPLOID_NO_CALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    private String sampleName = null;
    private List<Allele> alleles = Collections.emptyList();

    private boolean isPhased = false;
    private int GQ = -1;
    private int DP = -1;
    private int[] AD = null;
    private int[] PL = null;
    private Map<String, Object> extendedAttributes = null;
    private String filters = null;
    private int initialAttributeMapSize = 5;

    private final static Map<String, Object> NO_ATTRIBUTES =
            Collections.unmodifiableMap(new HashMap<String, Object>(0));

    // -----------------------------------------------------------------
    //
    // Factory methods
    //
    // -----------------------------------------------------------------

    public static Genotype create(final String sampleName, final List<Allele> alleles) {
        return new GenotypeBuilder(sampleName, alleles).make();
    }

    public static Genotype create(final String sampleName,
                                        final List<Allele> alleles,
                                        final Map<String, Object> attributes) {
        return new GenotypeBuilder(sampleName, alleles).attributes(attributes).make();
    }

    protected static Genotype create(final String sampleName,
                                           final List<Allele> alleles,
                                           final double[] gls) {
        return new GenotypeBuilder(sampleName, alleles).PL(gls).make();
    }

    /**
     * Create a new Genotype object for a sample that's missing from the VC (i.e., in
     * the output header).  Defaults to a diploid no call genotype ./.
     *
     * @param sampleName the name of this sample
     * @return an initialized Genotype with sampleName that's a diploid ./. no call genotype
     */
    public static Genotype createMissing(final String sampleName, final int ploidy) {
        final GenotypeBuilder builder = new GenotypeBuilder(sampleName);
        switch ( ploidy ) {
            case 1:  builder.alleles(HAPLOID_NO_CALL); break;
            case 2:  builder.alleles(DIPLOID_NO_CALL); break;
            default: builder.alleles(Collections.nCopies(ploidy, Allele.NO_CALL)); break;
        }
        return builder.make();
    }

    /**
     * Create a empty builder.  Both a sampleName and alleles must be provided
     * before trying to make a Genotype from this builder.
     */
    public GenotypeBuilder() {}

    /**
     * Create a builder using sampleName.  Alleles must be provided
     * before trying to make a Genotype from this builder.
     * @param sampleName
     */
    public GenotypeBuilder(final String sampleName) {
        name(sampleName);
    }

    /**
     * Make a builder using sampleName and alleles for starting values
     * @param sampleName
     * @param alleles
     */
    public GenotypeBuilder(final String sampleName, final List<Allele> alleles) {
        name(sampleName);
        alleles(alleles);
    }

    /**
     * Create a new builder starting with the values in Genotype g
     * @param g
     */
    public GenotypeBuilder(final Genotype g) {
        copy(g);
    }

    /**
     * Copy all of the values for this builder from Genotype g
     * @param g
     * @return
     */
    public GenotypeBuilder copy(final Genotype g) {
        name(g.getSampleName());
        alleles(g.getAlleles());
        phased(g.isPhased());
        GQ(g.getGQ());
        DP(g.getDP());
        AD(g.getAD());
        PL(g.getPL());
        filter(g.getFilters());
        attributes(g.getExtendedAttributes());
        return this;
    }

    /**
     * Reset all of the builder attributes to their defaults.  After this
     * function you must provide sampleName and alleles before trying to
     * make more Genotypes.
     */
    public final void reset(final boolean keepSampleName) {
        if ( ! keepSampleName ) sampleName = null;
        alleles = Collections.emptyList();
        isPhased = false;
        GQ = -1;
        DP = -1;
        AD = null;
        PL = null;
        filters = null;
        extendedAttributes = null;
    }

    /**
     * Create a new Genotype object using the values set in this builder.
     *
     * After creation the values in this builder can be modified and more Genotypes
     * created, althrough the contents of array values like PL should never be modified
     * inline as they are not copied for efficiency reasons.
     *
     * @return a newly minted Genotype object with values provided from this builder
     */
    @Ensures({"result != null"})
    public Genotype make() {
        final Map<String, Object> ea = extendedAttributes == null ? NO_ATTRIBUTES : extendedAttributes;
        return new FastGenotype(sampleName, alleles, isPhased, GQ, DP, AD, PL, filters, ea);
    }

    /**
     * Set this genotype's name
     * @param sampleName
     * @return
     */
    @Requires({"sampleName != null"})
    @Ensures({"this.sampleName != null"})
    public GenotypeBuilder name(final String sampleName) {
        this.sampleName = sampleName;
        return this;
    }

    /**
     * Set this genotype's alleles
     * @param alleles
     * @return
     */
    @Ensures({"this.alleles != null"})
    public GenotypeBuilder alleles(final List<Allele> alleles) {
        if ( alleles == null )
            this.alleles = Collections.emptyList();
        else
            this.alleles = alleles;
        return this;
    }

    /**
     * Is this genotype phased?
     * @param phased
     * @return
     */
    public GenotypeBuilder phased(final boolean phased) {
        isPhased = phased;
        return this;
    }

    @Requires({"GQ >= -1"})
    @Ensures({"this.GQ == GQ", "this.GQ >= -1"})
    public GenotypeBuilder GQ(final int GQ) {
        this.GQ = GQ;
        return this;
    }

    /**
     * Adaptor interface from the pLog10Error system.
     *
     * Will be retired when
     *
     * @param pLog10Error
     * @return
     */
    @Deprecated
    public GenotypeBuilder log10PError(final double pLog10Error) {
        if ( pLog10Error == CommonInfo.NO_LOG10_PERROR )
            return GQ(-1);
        else
            return GQ((int)Math.round(pLog10Error * -10));
    }

    /**
     * This genotype has no GQ value
     * @return
     */
    public GenotypeBuilder noGQ() { GQ = -1; return this; }

    /**
     * This genotype has no AD value
     * @return
     */
    public GenotypeBuilder noAD() { AD = null; return this; }

    /**
     * This genotype has no DP value
     * @return
     */
    public GenotypeBuilder noDP() { DP = -1; return this; }

    /**
     * This genotype has no PL value
     * @return
     */
    public GenotypeBuilder noPL() { PL = null; return this; }

    /**
     * This genotype has this DP value
     * @return
     */
    @Requires({"DP >= -1"})
    @Ensures({"this.DP == DP"})
    public GenotypeBuilder DP(final int DP) {
        this.DP = DP;
        return this;
    }

    /**
     * This genotype has this AD value
     * @return
     */
    @Requires({"AD == null || AD.length > 0"})
    @Ensures({"this.AD == AD"})
    public GenotypeBuilder AD(final int[] AD) {
        this.AD = AD;
        return this;
    }

    /**
     * This genotype has this PL value, as int[].  FAST
     * @return
     */
    @Requires("PL == null || PL.length > 0")
    @Ensures({"this.PL == PL"})
    public GenotypeBuilder PL(final int[] PL) {
        this.PL = PL;
        return this;
    }

    /**
     * This genotype has this PL value, converted from double[]. SLOW
     * @return
     */
    @Requires("PL == null || PL.length > 0")
    @Ensures({"this.PL == PL"})
    public GenotypeBuilder PL(final double[] GLs) {
        this.PL = GenotypeLikelihoods.fromLog10Likelihoods(GLs).getAsPLs();
        return this;
    }

    /**
     * This genotype has these attributes.
     *
     * Cannot contain inline attributes (DP, AD, GQ, PL)
     * @return
     */
    @Requires("attributes != null")
    @Ensures("attributes.isEmpty() || extendedAttributes != null")
    public GenotypeBuilder attributes(final Map<String, Object> attributes) {
        for ( Map.Entry<String, Object> pair : attributes.entrySet() )
            attribute(pair.getKey(), pair.getValue());
        return this;
    }

    /**
     * Tells this builder to remove all extended attributes
     *
     * @return
     */
    public GenotypeBuilder noAttributes() {
        this.extendedAttributes = null;
        return this;
    }

    /**
     * This genotype has this attribute key / value pair.
     *
     * Cannot contain inline attributes (DP, AD, GQ, PL)
     * @return
     */
    @Requires({"key != null"})
    @Ensures({"extendedAttributes != null", "extendedAttributes.containsKey(key)"})
    public GenotypeBuilder attribute(final String key, final Object value) {
        if ( extendedAttributes == null )
            extendedAttributes = new HashMap<String, Object>(initialAttributeMapSize);
        extendedAttributes.put(key, value);
        return this;
    }

    /**
     * Tells this builder to make a Genotype object that has had filters applied,
     * which may be empty (passes) or have some value indicating the reasons
     * why it's been filtered.
     *
     * @param filters non-null list of filters.  empty list => PASS
     * @return this builder
     */
    @Requires("filters != null")
    public GenotypeBuilder filters(final List<String> filters) {
        if ( filters.isEmpty() )
            return filter(null);
        else if ( filters.size() == 1 )
            return filter(filters.get(0));
        else
            return filter(ParsingUtils.join(";", ParsingUtils.sortList(filters)));
    }

    /**
     * varargs version of #filters
     * @param filters
     * @return
     */
    @Requires("filters != null")
    public GenotypeBuilder filters(final String ... filters) {
        return filters(Arrays.asList(filters));
    }

    /**
     * Most efficient version of setting filters -- just set the filters string to filters
     *
     * @param filter if filters == null or filters.equals("PASS") => genotype is PASS
     * @return
     */
    public GenotypeBuilder filter(final String filter) {
        this.filters = VCFConstants.PASSES_FILTERS_v4.equals(filter) ? null : filter;
        return this;
    }

    /**
     * This genotype is unfiltered
     *
     * @return
     */
    public GenotypeBuilder unfiltered() {
        return filter(null);
    }

    /**
     * Tell's this builder that we have at most these number of attributes
     * @return
     */
    public GenotypeBuilder maxAttributes(final int i) {
        initialAttributeMapSize = i;
        return this;
    }
}