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

import net.sf.samtools.util.StringUtil;

import java.util.Arrays;
import java.util.Collection;

/**
 * Immutable representation of an allele
 *
 * Types of alleles:
 *
 * Ref: a t C g a // C is the reference base
 *
 *    : a t G g a // C base is a G in some individuals
 *
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 *    : a t CAg a // A base is inserted w.r.t. the reference sequence
 *
 * In these cases, where are the alleles?
 *
 * SNP polymorphism of C/G  -> { C , G } -> C is the reference allele
 * 1 base deletion of C     -> { C , - } -> C is the reference allele
 * 1 base insertion of A    -> { - ; A } -> Null is the reference allele
 *
 * Suppose I see a the following in the population:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t G g a // C base is a G in some individuals
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 * How do I represent this?  There are three segregating alleles:
 *
 *  { C , G , - }
 *
 * Now suppose I have this more complex example:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t - g a
 *    : a t - - a
 *    : a t CAg a
 *
 * There are actually four segregating alleles:
 *
 *   { C g , - g, - -, and CAg } over bases 2-4
 *
 * However, the molecular equivalence explicitly listed above is usually discarded, so the actual
 * segregating alleles are:
 *
 *   { C g, g, -, C a g }
 *
 * Critically, it should be possible to apply an allele to a reference sequence to create the
 * correct haplotype sequence:
 *
 * Allele + reference => haplotype
 *
 * For convenience, we are going to create Alleles where the GenomeLoc of the allele is stored outside of the
 * Allele object itself.  So there's an idea of an A/C polymorphism independent of it's surrounding context.
 *
 * Given list of alleles it's possible to determine the "type" of the variation
 *
 *      A / C @ loc => SNP with
 *      - / A => INDEL
 *
 * If you know where allele is the reference, you can determine whether the variant is an insertion or deletion.
 *
 * Alelle also supports is concept of a NO_CALL allele.  This Allele represents a haplotype that couldn't be
 * determined. This is usually represented by a '.' allele.
 *
 * Note that Alleles store all bases as bytes, in **UPPER CASE**.  So 'atc' == 'ATC' from the perspective of an
 * Allele.

 * @author ebanks, depristo
 */
public class Allele implements Comparable<Allele> {
    private static final byte[] EMPTY_ALLELE_BASES = new byte[0];

    private boolean isRef = false;
    private boolean isNoCall = false;
    private boolean isSymbolic = false;

    private byte[] bases = null;

    public final static String NO_CALL_STRING = ".";
    /** A generic static NO_CALL allele for use */

    // no public way to create an allele
    protected Allele(byte[] bases, boolean isRef) {
        // null alleles are no longer allowed
        if ( wouldBeNullAllele(bases) ) {
            throw new IllegalArgumentException("Null alleles are not supported");
        }

        // no-calls are represented as no bases
        if ( wouldBeNoCallAllele(bases) ) {
            this.bases = EMPTY_ALLELE_BASES;
            isNoCall = true;
            if ( isRef ) throw new IllegalArgumentException("Cannot tag a NoCall allele as the reference allele");
            return;
        }

        if ( wouldBeSymbolicAllele(bases) ) {
            isSymbolic = true;
            if ( isRef ) throw new IllegalArgumentException("Cannot tag a symbolic allele as the reference allele");
        }
        else {
            StringUtil.toUpperCase(bases);
        }

        this.isRef = isRef;
        this.bases = bases;

        if ( ! acceptableAlleleBases(bases) )
            throw new IllegalArgumentException("Unexpected base in allele bases \'" + new String(bases)+"\'");
    }

    protected Allele(String bases, boolean isRef) {
        this(bases.getBytes(), isRef);
    }

    /**
     * Creates a new allele based on the provided one.  Ref state will be copied unless ignoreRefState is true
     * (in which case the returned allele will be non-Ref).
     *
     * This method is efficient because it can skip the validation of the bases (since the original allele was already validated)
     *
     * @param allele  the allele from which to copy the bases
     * @param ignoreRefState  should we ignore the reference state of the input allele and use the default ref state?
     */
    protected Allele(final Allele allele, final boolean ignoreRefState) {
        this.bases = allele.bases;
        this.isRef = ignoreRefState ? false : allele.isRef;
        this.isNoCall = allele.isNoCall;
        this.isSymbolic = allele.isSymbolic;
    }


    private final static Allele REF_A = new Allele("A", true);
    private final static Allele ALT_A = new Allele("A", false);
    private final static Allele REF_C = new Allele("C", true);
    private final static Allele ALT_C = new Allele("C", false);
    private final static Allele REF_G = new Allele("G", true);
    private final static Allele ALT_G = new Allele("G", false);
    private final static Allele REF_T = new Allele("T", true);
    private final static Allele ALT_T = new Allele("T", false);
    private final static Allele REF_N = new Allele("N", true);
    private final static Allele ALT_N = new Allele("N", false);
    public final static Allele NO_CALL = new Allele(NO_CALL_STRING, false);

    // ---------------------------------------------------------------------------------------------------------
    //
    // creation routines
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Create a new Allele that includes bases and if tagged as the reference allele if isRef == true.  If bases
     * == '-', a Null allele is created.  If bases ==  '.', a no call Allele is created.
     *
     * @param bases the DNA sequence of this variation, '-', of '.'
     * @param isRef should we make this a reference allele?
     * @throws IllegalArgumentException if bases contains illegal characters or is otherwise malformated
     */
    public static Allele create(byte[] bases, boolean isRef) {
        if ( bases == null )
            throw new IllegalArgumentException("create: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");

        if ( bases.length == 1 ) {
            // optimization to return a static constant Allele for each single base object
            switch (bases[0]) {
                case '.':
                    if ( isRef ) throw new IllegalArgumentException("Cannot tag a NoCall allele as the reference allele");
                    return NO_CALL;
                case 'A': case 'a' : return isRef ? REF_A : ALT_A;
                case 'C': case 'c' : return isRef ? REF_C : ALT_C;
                case 'G': case 'g' : return isRef ? REF_G : ALT_G;
                case 'T': case 't' : return isRef ? REF_T : ALT_T;
                case 'N': case 'n' : return isRef ? REF_N : ALT_N;
                default: throw new IllegalArgumentException("Illegal base [" + (char)bases[0] + "] seen in the allele");
            }
        } else {
            return new Allele(bases, isRef);
        }
    }

    public static Allele create(byte base, boolean isRef) {
        return create( new byte[]{ base }, isRef);
    }

    public static Allele create(byte base) {
        return create( base, false );
    }

    public static Allele extend(Allele left, byte[] right) {
        if (left.isSymbolic())
            throw new IllegalArgumentException("Cannot extend a symbolic allele");
        byte[] bases = new byte[left.length() + right.length];
        System.arraycopy(left.getBases(), 0, bases, 0, left.length());
        System.arraycopy(right, 0, bases, left.length(), right.length);

        return create(bases, left.isReference());
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the null allele
     */
    public static boolean wouldBeNullAllele(byte[] bases) {
        return (bases.length == 1 && bases[0] == '-') || bases.length == 0;
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the NO_CALL allele
     */
    public static boolean wouldBeNoCallAllele(byte[] bases) {
        return bases.length == 1 && bases[0] == '.';
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent a symbolic allele
     */
    public static boolean wouldBeSymbolicAllele(byte[] bases) {
        if ( bases.length <= 2 )
            return false;
        else {
            final String strBases = new String(bases);
            return (bases[0] == '<' && bases[bases.length-1] == '>') ||
                    (strBases.contains("[") || strBases.contains("]"));
        }
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the well formatted allele
     */
    public static boolean acceptableAlleleBases(String bases) {
        return acceptableAlleleBases(bases.getBytes(), true);
    }

    public static boolean acceptableAlleleBases(String bases, boolean allowNsAsAcceptable) {
        return acceptableAlleleBases(bases.getBytes(), allowNsAsAcceptable);
    }

    /**
     * @param bases  bases representing an allele
     * @return true if the bases represent the well formatted allele
     */
    public static boolean acceptableAlleleBases(byte[] bases) {
        return acceptableAlleleBases(bases, true); // default: N bases are acceptable
    }
    
    public static boolean acceptableAlleleBases(byte[] bases, boolean allowNsAsAcceptable) {
        if ( wouldBeNullAllele(bases) )
            return false;

        if ( wouldBeNoCallAllele(bases) || wouldBeSymbolicAllele(bases) )
            return true;

        for (byte base :  bases ) {
            switch (base) {
                case 'A': case 'C': case 'G': case 'T':  case 'a': case 'c': case 'g': case 't': 
                    break;
                case 'N' : case 'n' :
                    if (allowNsAsAcceptable)
                        break;
                    else
                        return false;
                default:
                    return false;
            }
        }

        return true;
    }

    /**
     * @see Allele(byte[], boolean)
     *
     * @param bases  bases representing an allele
     * @param isRef  is this the reference allele?
     */
    public static Allele create(String bases, boolean isRef) {
        return create(bases.getBytes(), isRef);
    }


    /**
     * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
     *
     * @param bases  bases representing an allele
     */
    public static Allele create(String bases) {
        return create(bases, false);
    }

    /**
     * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
     *
     * @param bases  bases representing an allele
     */
    public static Allele create(byte[] bases) {
        return create(bases, false);
    }

    /**
     * Creates a new allele based on the provided one.  Ref state will be copied unless ignoreRefState is true
     * (in which case the returned allele will be non-Ref).
     *
     * This method is efficient because it can skip the validation of the bases (since the original allele was already validated)
     *
     * @param allele  the allele from which to copy the bases
     * @param ignoreRefState  should we ignore the reference state of the input allele and use the default ref state?
     */
    public static Allele create(final Allele allele, final boolean ignoreRefState) {
        return new Allele(allele, ignoreRefState);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // accessor routines
    //
    // ---------------------------------------------------------------------------------------------------------

    // Returns true if this is the NO_CALL allele
    public boolean isNoCall()           { return isNoCall; }
    // Returns true if this is not the NO_CALL allele
    public boolean isCalled()           { return ! isNoCall(); }

    // Returns true if this Allele is the reference allele
    public boolean isReference()        { return isRef; }
    // Returns true if this Allele is not the reference allele
    public boolean isNonReference()     { return ! isReference(); }

    // Returns true if this Allele is symbolic (i.e. no well-defined base sequence)
    public boolean isSymbolic()         { return isSymbolic; }

    // Returns a nice string representation of this object
    public String toString() {
        return ( isNoCall() ? NO_CALL_STRING : getDisplayString() ) + (isReference() ? "*" : "");
    }

    /**
     * Return the DNA bases segregating in this allele.  Note this isn't reference polarized,
     * so the Null allele is represented by a vector of length 0
     *
     * @return the segregating bases
     */
    public byte[] getBases() { return isSymbolic ? EMPTY_ALLELE_BASES : bases; }

    /**
     * Return the DNA bases segregating in this allele in String format.
     * This is useful, because toString() adds a '*' to reference alleles and getBases() returns garbage when you call toString() on it.
     *
     * @return the segregating bases
     */
    public String getBaseString() { return isNoCall() ? NO_CALL_STRING : new String(getBases()); }

    /**
     * Return the printed representation of this allele.
     * Same as getBaseString(), except for symbolic alleles.
     * For symbolic alleles, the base string is empty while the display string contains <TAG>.
     *
     * @return the allele string representation
     */
    public String getDisplayString() { return new String(bases); }

    /**
     * Same as #getDisplayString() but returns the result as byte[].
     *
     * Slightly faster then getDisplayString()
     *
     * @return the allele string representation
     */
    public byte[] getDisplayBases() { return bases; }

    /**
     * @param other  the other allele
     *
     * @return true if these alleles are equal
     */
    public boolean equals(Object other) {
        return ( ! (other instanceof Allele) ? false : equals((Allele)other, false) );
    }

    /**
     * @return hash code
     */
    public int hashCode() {
        int hash = 1;
        for (int i = 0; i < bases.length; i++)
            hash += (i+1) * bases[i];
        return hash;
    }

    /**
     * Returns true if this and other are equal.  If ignoreRefState is true, then doesn't require both alleles has the
     * same ref tag
     *
     * @param other            allele to compare to
     * @param ignoreRefState   if true, ignore ref state in comparison
     * @return true if this and other are equal
     */
    public boolean equals(Allele other, boolean ignoreRefState) {
        return this == other || (isRef == other.isRef || ignoreRefState) && isNoCall == other.isNoCall && (bases == other.bases || Arrays.equals(bases, other.bases));
    }

    /**
     * @param test  bases to test against
     *
     * @return  true if this Allele contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(byte[] test) { return !isSymbolic && (bases == test || Arrays.equals(bases, test)); }

    /**
     * @param test  bases to test against
     *
     * @return  true if this Allele contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(String test) { return basesMatch(test.toUpperCase().getBytes()); }

    /**
     * @param test  allele to test against
     *
     * @return  true if this Allele contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
     */
    public boolean basesMatch(Allele test) { return basesMatch(test.getBases()); }

    /**
     * @return the length of this allele.  Null and NO_CALL alleles have 0 length.
     */
    public int length() {
        return isSymbolic ? 0 : bases.length;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // useful static functions
    //
    // ---------------------------------------------------------------------------------------------------------

    public static Allele getMatchingAllele(Collection<Allele> allAlleles, byte[] alleleBases) {
        for ( Allele a : allAlleles ) {
            if ( a.basesMatch(alleleBases) ) {
                return a;
            }
        }

        if ( wouldBeNoCallAllele(alleleBases) )
            return NO_CALL;
        else
            return null;    // couldn't find anything
    }

    public int compareTo(Allele other) {
        if ( isReference() && other.isNonReference() )
            return -1;
        else if ( isNonReference() && other.isReference() ) 
            return 1;
        else
            return getBaseString().compareTo(other.getBaseString()); // todo -- potential performance issue
    }

    public static boolean oneIsPrefixOfOther(Allele a1, Allele a2) {
        if ( a2.length() >= a1.length() )
            return firstIsPrefixOfSecond(a1, a2);
        else
            return firstIsPrefixOfSecond(a2, a1);
    }

    private static boolean firstIsPrefixOfSecond(Allele a1, Allele a2) {
        String a1String = a1.getBaseString();
        return a2.getBaseString().substring(0, a1String.length()).equals(a1String);
    }
}
