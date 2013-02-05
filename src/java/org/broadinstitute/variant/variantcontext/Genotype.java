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
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * @author Mark DePristo
 */
@Invariant({
        "getAlleles() != null",
        "getSampleName() != null",
        "getPloidy() >= 0",
        "! hasForbiddenKey(getExtendedAttributes())"})
public abstract class Genotype implements Comparable<Genotype> {
    /**
     * A list of genotype field keys corresponding to values we
     * manage inline in the Genotype object.  They must not appear in the
     * extended attributes map
     */
    public final static Collection<String> PRIMARY_KEYS = Arrays.asList(
            VCFConstants.GENOTYPE_FILTER_KEY,
            VCFConstants.GENOTYPE_KEY,
            VCFConstants.GENOTYPE_QUALITY_KEY,
            VCFConstants.DEPTH_KEY,
            VCFConstants.GENOTYPE_ALLELE_DEPTHS,
            VCFConstants.GENOTYPE_PL_KEY);

    public final static String PHASED_ALLELE_SEPARATOR = "|";
    public final static String UNPHASED_ALLELE_SEPARATOR = "/";

    private final String sampleName;
    private GenotypeType type = null;
    private final String filters;

    protected Genotype(final String sampleName, final String filters) {
        this.sampleName = sampleName;
        this.filters = filters;
    }

    /**
     * @return the alleles for this genotype.  Cannot be null.  May be empty
     */
    @Ensures("result != null")
    public abstract List<Allele> getAlleles();

    /**
     * Returns how many times allele appears in this genotype object?
     *
     * @param allele
     * @return a value >= 0 indicating how many times the allele occurred in this sample's genotype
     */
    @Requires("allele != null")
    @Ensures("result >= 0")
    public int countAllele(final Allele allele) {
        int c = 0;
        for ( final Allele a : getAlleles() )
            if ( a.equals(allele) )
                c++;

        return c;
    }

    /**
     * Get the ith allele in this genotype
     *
     * @param i the ith allele, must be < the ploidy, starting with 0
     * @return the allele at position i, which cannot be null
     */
    @Requires({"i >=0 && i < getPloidy()", "getType() != GenotypeType.UNAVAILABLE"})
    @Ensures("result != null")
    public abstract Allele getAllele(int i);

    /**
     * Are the alleles phased w.r.t. the global phasing system?
     *
     * @return true if yes
     */
    public abstract boolean isPhased();

    /**
     * What is the ploidy of this sample?
     *
     * @return the ploidy of this genotype.  0 if the site is no-called.
     */
    @Ensures("result >= 0")
    public int getPloidy() {
        return getAlleles().size();
    }

    /**
     * @return the sequencing depth of this sample, or -1 if this value is missing
     */
    @Ensures("result >= -1")
    public abstract int getDP();

    /**
     * @return the count of reads, one for each allele in the surrounding Variant context,
     *      matching the corresponding allele, or null if this value is missing.  MUST
     *      NOT BE MODIFIED!
     */
    public abstract int[] getAD();

    /**
     * Returns the name associated with this sample.
     *
     * @return a non-null String
     */
    @Ensures("result != null")
    public String getSampleName() {
        return sampleName;
    }

    /**
     * Returns a phred-scaled quality score, or -1 if none is available
     * @return
     */
    @Ensures("result >= -1")
    public abstract int getGQ();

    /**
     * Does the PL field have a value?
     * @return true if there's a PL field value
     */
    @Ensures("(result == false && getPL() == null) || (result == true && getPL() != null)")
    public boolean hasPL() {
        return getPL() != null;
    }

    /**
     * Does the AD field have a value?
     * @return true if there's a AD field value
     */
    @Ensures("(result == false && getAD() == null) || (result == true && getAD() != null)")
    public boolean hasAD() {
        return getAD() != null;
    }

    /**
     * Does the GQ field have a value?
     * @return true if there's a GQ field value
     */
    @Ensures("(result == false && getGQ() == -1) || (result == true && getGQ() >= 0)")
    public boolean hasGQ() {
        return getGQ() != -1;
    }

    /**
     * Does the DP field have a value?
     * @return true if there's a DP field value
     */
    @Ensures("(result == false && getDP() == -1) || (result == true && getDP() >= 0)")
    public boolean hasDP() {
        return getDP() != -1;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // The type of this genotype
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the high-level type of this sample's genotype
     */
    @Ensures({"type != null", "result != null"})
    public GenotypeType getType() {
        if ( type == null ) {
            type = determineType();
        }
        return type;
    }

    /**
     * Internal code to determine the type of the genotype from the alleles vector
     * @return the type
     */
    @Requires("type == null") // we should never call if already calculated
    protected GenotypeType determineType() {
        // TODO -- this code is slow and could be optimized for the diploid case
        final List<Allele> alleles = getAlleles();
        if ( alleles.isEmpty() )
            return GenotypeType.UNAVAILABLE;

        boolean sawNoCall = false, sawMultipleAlleles = false;
        Allele observedAllele = null;

        for ( final Allele allele : alleles ) {
            if ( allele.isNoCall() )
                sawNoCall = true;
            else if ( observedAllele == null )
                observedAllele = allele;
            else if ( !allele.equals(observedAllele) )
                sawMultipleAlleles = true;
        }

        if ( sawNoCall ) {
            if ( observedAllele == null )
                return GenotypeType.NO_CALL;
            return GenotypeType.MIXED;
        }

        if ( observedAllele == null )
            throw new IllegalStateException("BUG: there are no alleles present in this genotype but the alleles list is not null");

        return sawMultipleAlleles ? GenotypeType.HET : observedAllele.isReference() ? GenotypeType.HOM_REF : GenotypeType.HOM_VAR;
    }

    /**
     * @return true if all observed alleles are the same (regardless of whether they are ref or alt); if any alleles are no-calls, this method will return false.
     */
    public boolean isHom()    { return isHomRef() || isHomVar(); }

    /**
     * @return true if all observed alleles are ref; if any alleles are no-calls, this method will return false.
     */
    public boolean isHomRef() { return getType() == GenotypeType.HOM_REF; }

    /**
     * @return true if all observed alleles are alt; if any alleles are no-calls, this method will return false.
     */
    public boolean isHomVar() { return getType() == GenotypeType.HOM_VAR; }

    /**
     * @return true if we're het (observed alleles differ); if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
     */
    public boolean isHet() { return getType() == GenotypeType.HET; }

    /**
     * @return true if this genotype is not actually a genotype but a "no call" (e.g. './.' in VCF); if any alleles are not no-calls (even if some are), this method will return false.
     */
    public boolean isNoCall() { return getType() == GenotypeType.NO_CALL; }

    /**
     * @return true if this genotype is comprised of any alleles that are not no-calls (even if some are).
     */
    public boolean isCalled() { return getType() != GenotypeType.NO_CALL && getType() != GenotypeType.UNAVAILABLE; }

    /**
     * @return true if this genotype is comprised of both calls and no-calls.
     */
    public boolean isMixed() { return getType() == GenotypeType.MIXED; }

    /**
     * @return true if the type of this genotype is set.
     */
    public boolean isAvailable() { return getType() != GenotypeType.UNAVAILABLE; }

    // ------------------------------------------------------------------------------
    //
    // methods for getting genotype likelihoods for a genotype object, if present
    //
    // ------------------------------------------------------------------------------

    /**
     * @return Returns true if this Genotype has PL field values
     */
    @Ensures("(result && getLikelihoods() != null) || (! result && getLikelihoods() == null)")
    public boolean hasLikelihoods() {
        return getPL() != null;
    }

    /**
     * Convenience function that returns a string representation of the PL field of this
     * genotype, or . if none is available.
     *
     * @return a non-null String representation for the PL of this sample
     */
    @Ensures("result != null")
    public String getLikelihoodsString() {
        return hasLikelihoods() ? getLikelihoods().toString() : VCFConstants.MISSING_VALUE_v4;
    }

    /**
     * Returns the GenotypesLikelihoods data associated with this Genotype, or null if missing
     * @return null or a GenotypesLikelihood object for this sample's PL field
     */
    @Ensures("(hasLikelihoods() && result != null) || (! hasLikelihoods() && result == null)")
    public GenotypeLikelihoods getLikelihoods() {
        return hasLikelihoods() ? GenotypeLikelihoods.fromPLs(getPL()) : null;
    }

    /**
     * Are all likelihoods for this sample non-informative?
     *
     * Returns true if all PLs are 0 => 0,0,0 => true
     * 0,0,0,0,0,0 => true
     * 0,10,100 => false
     *
     * @return true if all samples PLs are equal and == 0
     */
    public boolean isNonInformative() {
        if ( getPL() == null )
            return true;
        else {
            for ( final int PL : getPL() ) {
                if ( PL != 0 )
                    return false;
            }

            return true;
        }
    }

    /**
     * Unsafe low-level accessor the PL field itself, may be null.
     *
     * @return a pointer to the underlying PL data.  MUST NOT BE MODIFIED!
     */
    public abstract int[] getPL();

    // ---------------------------------------------------------------------------------------------------------
    //
    // Many different string representations
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Return a VCF-like string representation for the alleles of this genotype.
     *
     * Does not append the reference * marker on the alleles.
     *
     * @return a string representing the genotypes, or null if the type is unavailable.
     */
    @Ensures("result != null || ! isAvailable()")
    public String getGenotypeString() {
        return getGenotypeString(true);
    }

    /**
     * Return a VCF-like string representation for the alleles of this genotype.
     *
     * If ignoreRefState is true, will not append the reference * marker on the alleles.
     *
     * @return a string representing the genotypes, or null if the type is unavailable.
     */
    @Ensures("result != null || ! isAvailable()")
    public String getGenotypeString(boolean ignoreRefState) {
        if ( getPloidy() == 0 )
            return "NA";

        // Notes:
        // 1. Make sure to use the appropriate separator depending on whether the genotype is phased
        // 2. If ignoreRefState is true, then we want just the bases of the Alleles (ignoring the '*' indicating a ref Allele)
        // 3. So that everything is deterministic with regards to integration tests, we sort Alleles (when the genotype isn't phased, of course)
        return ParsingUtils.join(isPhased() ? PHASED_ALLELE_SEPARATOR : UNPHASED_ALLELE_SEPARATOR,
                ignoreRefState ? getAlleleStrings() : (isPhased() ? getAlleles() : ParsingUtils.sortList(getAlleles())));
    }

    /**
     * Utility that returns a list of allele strings corresponding to the alleles in this sample
     * @return
     */
    protected List<String> getAlleleStrings() {
        final List<String> al = new ArrayList<String>(getPloidy());
        for ( Allele a : getAlleles() )
            al.add(a.getBaseString());

        return al;
    }

    public String toString() {
        return String.format("[%s %s%s%s%s%s%s%s]",
                getSampleName(),
                getGenotypeString(false),
                toStringIfExists(VCFConstants.GENOTYPE_QUALITY_KEY, getGQ()),
                toStringIfExists(VCFConstants.DEPTH_KEY, getDP()),
                toStringIfExists(VCFConstants.GENOTYPE_ALLELE_DEPTHS, getAD()),
                toStringIfExists(VCFConstants.GENOTYPE_PL_KEY, getPL()),
                toStringIfExists(VCFConstants.GENOTYPE_FILTER_KEY, getFilters()),
                sortedString(getExtendedAttributes()));
    }

    public String toBriefString() {
        return String.format("%s:Q%d", getGenotypeString(false), getGQ());
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Comparison operations
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * comparable genotypes -> compareTo on the sample names
     * @param genotype
     * @return
     */
    @Override
    public int compareTo(final Genotype genotype) {
        return getSampleName().compareTo(genotype.getSampleName());
    }

    public boolean sameGenotype(final Genotype other) {
        return sameGenotype(other, true);
    }

    public boolean sameGenotype(final Genotype other, boolean ignorePhase) {
        if (getPloidy() != other.getPloidy())
            return false; // gotta have the same number of allele to be equal

        // By default, compare the elements in the lists of alleles, element-by-element
        Collection<Allele> thisAlleles = this.getAlleles();
        Collection<Allele> otherAlleles = other.getAlleles();

        if (ignorePhase) { // do not care about order, only identity of Alleles
            thisAlleles = new TreeSet<Allele>(thisAlleles);   //implemented Allele.compareTo()
            otherAlleles = new TreeSet<Allele>(otherAlleles);
        }

        return thisAlleles.equals(otherAlleles);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // get routines for extended attributes
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Returns the extended attributes for this object
     * @return is never null, but is often isEmpty()
     */
    @Ensures({"result != null", "! hasForbiddenKey(result)"})
    public abstract Map<String, Object> getExtendedAttributes();

    /**
     * Is key associated with a value (even a null one) in the extended attributes?
     *
     * Note this will not return true for the inline attributes DP, GQ, AD, or PL
     *
     * @param key a non-null string key to check for an association
     * @return true if key has a value in the extendedAttributes
     */
    @Requires({"key != null", "! isForbiddenKey(key)"})
    public boolean hasExtendedAttribute(final String key) {
        return getExtendedAttributes().containsKey(key);
    }

    /**
     * Get the extended attribute value associated with key, if possible
     *
     * @param key a non-null string key to fetch a value for
     * @param defaultValue the value to return if key isn't in the extended attributes
     * @return a value (potentially) null associated with key, or defaultValue if no association exists
     */
    @Requires({"key != null", "! isForbiddenKey(key)"})
    @Ensures("hasExtendedAttribute(key) || result == defaultValue")
    public Object getExtendedAttribute(final String key, final Object defaultValue) {
        return hasExtendedAttribute(key) ? getExtendedAttributes().get(key) : defaultValue;
    }

    /**
     * Same as #getExtendedAttribute with a null default
     *
     * @param key
     * @return
     */
    public Object getExtendedAttribute(final String key) {
        return getExtendedAttribute(key, null);
    }

    /**
     * Returns the filter string associated with this Genotype.
     *
     * @return If this result == null, then the genotype is considered PASSing filters
     *   If the result != null, then the genotype has failed filtering for the reason(s)
     *   specified in result.  To be reference compliant multiple filter field
     *   string values can be encoded with a ; separator.
     */
    public final String getFilters() {
        return filters;
    }

    /**
     * Is this genotype filtered or not?
     *
     * @return returns false if getFilters() == null
     */
    @Ensures({"result != (getFilters() == null)"})
    public final boolean isFiltered() {
        return getFilters() != null;
    }

    @Deprecated public boolean hasLog10PError() { return hasGQ(); }
    @Deprecated public double getLog10PError() { return getGQ() / -10.0; }
    @Deprecated public int getPhredScaledQual() { return getGQ(); }

    @Deprecated
    public String getAttributeAsString(String key, String defaultValue) {
        Object x = getExtendedAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof String ) return (String)x;
        return String.valueOf(x); // throws an exception if this isn't a string
    }

    @Deprecated
    public int getAttributeAsInt(String key, int defaultValue) {
        Object x = getExtendedAttribute(key);
        if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
        if ( x instanceof Integer ) return (Integer)x;
        return Integer.valueOf((String)x); // throws an exception if this isn't a string
    }

    @Deprecated
    public double getAttributeAsDouble(String key, double defaultValue) {
        Object x = getExtendedAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Double ) return (Double)x;
        return Double.valueOf((String)x); // throws an exception if this isn't a string
    }

    /**
     * A totally generic getter, that allows you to specific keys that correspond
     * to even inline values (GQ, for example).  Can be very expensive.  Additionally,
     * all int[] are converted inline into List<Integer> for convenience.
     *
     * @param key
     * @return
     */
    public Object getAnyAttribute(final String key) {
        if (key.equals(VCFConstants.GENOTYPE_KEY)) {
            return getAlleles();
        } else if (key.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
            return getGQ();
        } else if (key.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS)) {
            return Arrays.asList(getAD());
        } else if (key.equals(VCFConstants.GENOTYPE_PL_KEY)) {
            return Arrays.asList(getPL());
        } else if (key.equals(VCFConstants.DEPTH_KEY)) {
            return getDP();
        } else {
            return getExtendedAttribute(key);
        }
    }

    public boolean hasAnyAttribute(final String key) {
        if (key.equals(VCFConstants.GENOTYPE_KEY)) {
            return isAvailable();
        } else if (key.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
            return hasGQ();
        } else if (key.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS)) {
            return hasAD();
        } else if (key.equals(VCFConstants.GENOTYPE_PL_KEY)) {
            return hasPL();
        } else if (key.equals(VCFConstants.DEPTH_KEY)) {
            return hasDP();
        } else {
            return hasExtendedAttribute(key);
        }
    }

    // TODO -- add getAttributesAsX interface here

    // ------------------------------------------------------------------------------
    //
    // private utilities
    //
    // ------------------------------------------------------------------------------

    /**
     * a utility method for generating sorted strings from a map key set.
     * @param c the map
     * @param <T> the key type
     * @param <V> the value type
     * @return a sting, enclosed in {}, with comma seperated key value pairs in order of the keys
     */
    @Requires("c != null")
    protected static <T extends Comparable<T>, V> String sortedString(Map<T, V> c) {

        // NOTE -- THIS IS COPIED FROM GATK UTILS TO ALLOW US TO KEEP A SEPARATION BETWEEN THE GATK AND VCF CODECS
        final List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);

        final List<String> pairs = new ArrayList<String>();
        for (final T k : t) {
            pairs.add(k + "=" + c.get(k));
        }

        return pairs.isEmpty() ? "" : " {" + ParsingUtils.join(", ", pairs.toArray(new String[pairs.size()])) + "}";
    }

    /**
     * Returns a display name for field name with value v if this isn't -1.  Otherwise returns ""
     * @param name of the field ("AD")
     * @param v the value of the field, or -1 if missing
     * @return a non-null string for display if the field is not missing
     */
    @Requires("name != null")
    @Ensures("result != null")
    protected final static String toStringIfExists(final String name, final int v) {
        return v == -1 ? "" : " " + name + " " + v;
    }

    /**
     * Returns a display name for field name with String value v if this isn't null.  Otherwise returns ""
     * @param name of the field ("FT")
     * @param v the value of the field, or null if missing
     * @return a non-null string for display if the field is not missing
     */
    protected final static String toStringIfExists(final String name, final String v) {
        return v == null ? "" : " " + name + " " + v;
    }

    /**
     * Returns a display name for field name with values vs if this isn't null.  Otherwise returns ""
     * @param name of the field ("AD")
     * @param vs the value of the field, or null if missing
     * @return a non-null string for display if the field is not missing
     */
    @Requires("name != null")
    @Ensures("result != null")
    protected final static String toStringIfExists(final String name, final int[] vs) {
        if ( vs == null )
            return "";
        else {
            StringBuilder b = new StringBuilder();
            b.append(" ").append(name).append(" ");
            for ( int i = 0; i < vs.length; i++ ) {
                if ( i != 0 ) b.append(",");
                b.append(vs[i]);
            }
            return b.toString();
        }
    }

    /**
     * Does the attribute map have a mapping involving a forbidden key (i.e.,
     * one that's managed inline by this Genotypes object?
     *
     * @param attributes the extended attributes key
     * @return
     */
    protected final static boolean hasForbiddenKey(final Map<String, Object> attributes) {
        for ( final String forbidden : PRIMARY_KEYS)
            if ( attributes.containsKey(forbidden) )
                return true;
        return false;
    }

    protected final static boolean isForbiddenKey(final String key) {
        return PRIMARY_KEYS.contains(key);
    }
}