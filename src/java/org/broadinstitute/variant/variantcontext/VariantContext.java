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

import org.broad.tribble.Feature;
import org.broad.tribble.TribbleException;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.vcf.*;

import java.util.*;

/**
 * Class VariantContext
 *
 * == High-level overview ==
 *
 * The VariantContext object is a single general class system for representing genetic variation data composed of:
 *
 * * Allele: representing single genetic haplotypes (A, T, ATC, -)
 * * Genotype: an assignment of alleles for each chromosome of a single named sample at a particular locus
 * * VariantContext: an abstract class holding all segregating alleles at a locus as well as genotypes
 *    for multiple individuals containing alleles at that locus
 *
 * The class system works by defining segregating alleles, creating a variant context representing the segregating
 * information at a locus, and potentially creating and associating genotypes with individuals in the context.
 *
 * All of the classes are highly validating -- call validate() if you modify them -- so you can rely on the
 * self-consistency of the data once you have a VariantContext in hand.  The system has a rich set of assessor
 * and manipulator routines, as well as more complex static support routines in VariantContextUtils.
 *
 * The VariantContext (and Genotype) objects are attributed (supporting addition of arbitrary key/value pairs) and
 * filtered (can represent a variation that is viewed as suspect).
 *
 * VariantContexts are dynamically typed, so whether a VariantContext is a SNP, Indel, or NoVariant depends
 * on the properties of the alleles in the context.  See the detailed documentation on the Type parameter below.
 *
 * It's also easy to create subcontexts based on selected genotypes.
 *
 * == Working with Variant Contexts ==
 * By default, VariantContexts are immutable.  In order to access (in the rare circumstances where you need them)
 * setter routines, you need to create MutableVariantContexts and MutableGenotypes.
 *
 * === Some example data ===
 *
 * Allele A, Aref, T, Tref;
 * Allele del, delRef, ATC, ATCref;
 *
 * A [ref] / T at 10
 * GenomeLoc snpLoc = GenomeLocParser.createGenomeLoc("chr1", 10, 10);
 *
 * - / ATC [ref] from 20-23
 * GenomeLoc delLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 22);
 *
 *  // - [ref] / ATC immediately after 20
 * GenomeLoc insLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 20);
 *
 * === Alleles ===
 *
 * See the documentation in the Allele class itself
 *
 * What are they?
 *
 * Alleles can be either reference or non-reference
 *
 * Example alleles used here:
 *
 *   del = new Allele("-");
 *   A = new Allele("A");
 *   Aref = new Allele("A", true);
 *   T = new Allele("T");
 *   ATC = new Allele("ATC");
 *
 * === Creating variant contexts ===
 *
 * ==== By hand ====
 *
 * Here's an example of a A/T polymorphism with the A being reference:
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, snpLoc, Arrays.asList(Aref, T));
 * </pre>
 *
 * If you want to create a non-variant site, just put in a single reference allele
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, snpLoc, Arrays.asList(Aref));
 * </pre>
 *
 * A deletion is just as easy:
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, delLoc, Arrays.asList(ATCref, del));
 * </pre>
 *
 * The only 2 things that distinguishes between a insertion and deletion are the reference allele
 * and the location of the variation.  An insertion has a Null reference allele and at least
 * one non-reference Non-Null allele.  Additionally, the location of the insertion is immediately after
 * a 1-bp GenomeLoc (at say 20).
 *
 * <pre>
 * VariantContext vc = new VariantContext("name", insLoc, Arrays.asList(delRef, ATC));
 * </pre>
 *
 * ==== Converting rods and other data structures to VCs ====
 *
 * You can convert many common types into VariantContexts using the general function:
 *
 * <pre>
 * VariantContextAdaptors.convertToVariantContext(name, myObject)
 * </pre>
 *
 * dbSNP and VCFs, for example, can be passed in as myObject and a VariantContext corresponding to that
 * object will be returned.  A null return type indicates that the type isn't yet supported.  This is the best
 * and easiest way to create contexts using RODs.
 *
 *
 * === Working with genotypes ===
 *
 * <pre>
 * List<Allele> alleles = Arrays.asList(Aref, T);
 * Genotype g1 = new Genotype(Arrays.asList(Aref, Aref), "g1", 10);
 * Genotype g2 = new Genotype(Arrays.asList(Aref, T), "g2", 10);
 * Genotype g3 = new Genotype(Arrays.asList(T, T), "g3", 10);
 * VariantContext vc = new VariantContext(snpLoc, alleles, Arrays.asList(g1, g2, g3));
 * </pre>
 *
 * At this point we have 3 genotypes in our context, g1-g3.
 *
 * You can assess a good deal of information about the genotypes through the VariantContext:
 *
 * <pre>
 * vc.hasGenotypes()
 * vc.isMonomorphicInSamples()
 * vc.isPolymorphicInSamples()
 * vc.getSamples().size()
 *
 * vc.getGenotypes()
 * vc.getGenotypes().get("g1")
 * vc.hasGenotype("g1")
 *
 * vc.getCalledChrCount()
 * vc.getCalledChrCount(Aref)
 * vc.getCalledChrCount(T)
 * </pre>
 *
 * === NO_CALL alleles ===
 *
 * The system allows one to create Genotypes carrying special NO_CALL alleles that aren't present in the
 * set of context alleles and that represent undetermined alleles in a genotype:
 *
 * Genotype g4 = new Genotype(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), "NO_DATA_FOR_SAMPLE", 10);
 *
 *
 * === subcontexts ===
 * It's also very easy get subcontext based only the data in a subset of the genotypes:
 *
 * <pre>
 * VariantContext vc12 = vc.subContextFromGenotypes(Arrays.asList(g1,g2));
 * VariantContext vc1 = vc.subContextFromGenotypes(Arrays.asList(g1));
 * </pre>
 *
 * <s3>
 *     Fully decoding.  Currently VariantContexts support some fields, particularly those
 *     stored as generic attributes, to be of any type.  For example, a field AB might
 *     be naturally a floating point number, 0.51, but when it's read into a VC its
 *     not decoded into the Java presentation but left as a string "0.51".  A fully
 *     decoded VariantContext is one where all values have been converted to their
 *     corresponding Java object types, based on the types declared in a VCFHeader.
 *
 *     The fullyDecode() takes a header object and creates a new fully decoded VariantContext
 *     where all fields are converted to their true java representation.  The VCBuilder
 *     can be told that all fields are fully decoded, in which case no work is done when
 *     asking for a fully decoded version of the VC.
 * </s3>
 *
 * @author depristo
 */
public class VariantContext implements Feature { // to enable tribble integration
    private final static boolean WARN_ABOUT_BAD_END = true;
    private final static int MAX_ALLELE_SIZE_FOR_NON_SV = 150;
    private boolean fullyDecoded = false;
    protected CommonInfo commonInfo = null;
    public final static double NO_LOG10_PERROR = CommonInfo.NO_LOG10_PERROR;

    public final static Set<String> PASSES_FILTERS = Collections.unmodifiableSet(new LinkedHashSet<String>());

    /** The location of this VariantContext */
    final protected String contig;
    final protected long start;
    final protected long stop;
    private final String ID;

    /** The type (cached for performance reasons) of this context */
    protected Type type = null;

    /** A set of the alleles segregating in this context */
    final protected List<Allele> alleles;

    /** A mapping from sampleName -> genotype objects for all genotypes associated with this context */
    protected GenotypesContext genotypes = null;

    /** Counts for each of the possible Genotype types in this context */
    protected int[] genotypeCounts = null;

    public final static GenotypesContext NO_GENOTYPES = GenotypesContext.NO_GENOTYPES;

    // a fast cached access point to the ref / alt alleles for biallelic case
    private Allele REF = null;

    // set to the alt allele when biallelic, otherwise == null
    private Allele ALT = null;

    /* cached monomorphic value: null -> not yet computed, False, True */
    private Boolean monomorphic = null;

    // ---------------------------------------------------------------------------------------------------------
    //
    // validation mode
    //
    // ---------------------------------------------------------------------------------------------------------

    public enum Validation {
        ALLELES,
        GENOTYPES
    }

    private final static EnumSet<Validation> NO_VALIDATION = EnumSet.noneOf(Validation.class);

    // ---------------------------------------------------------------------------------------------------------
    //
    // constructors: see VariantContextBuilder
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Copy constructor
     *
     * @param other the VariantContext to copy
     */
    protected VariantContext(VariantContext other) {
        this(other.getSource(), other.getID(), other.getChr(), other.getStart(), other.getEnd(),
                other.getAlleles(), other.getGenotypes(), other.getLog10PError(),
                other.getFiltersMaybeNull(),
                other.getAttributes(),
                other.fullyDecoded, NO_VALIDATION);
    }

    /**
     * the actual constructor.  Private access only
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles         alleles
     * @param genotypes       genotypes map
     * @param log10PError  qual
     * @param filters         filters: use null for unfiltered and empty set for passes filters
     * @param attributes      attributes
     * @param validationToPerform     set of validation steps to take
     */
    protected VariantContext(final String source,
                             final String ID,
                             final String contig,
                             final long start,
                             final long stop,
                             final Collection<Allele> alleles,
                             final GenotypesContext genotypes,
                             final double log10PError,
                             final Set<String> filters,
                             final Map<String, Object> attributes,
                             final boolean fullyDecoded,
                             final EnumSet<Validation> validationToPerform ) {
        if ( contig == null ) { throw new IllegalArgumentException("Contig cannot be null"); }
        this.contig = contig;
        this.start = start;
        this.stop = stop;

        // intern for efficiency.  equals calls will generate NPE if ID is inappropriately passed in as null
        if ( ID == null || ID.equals("") ) throw new IllegalArgumentException("ID field cannot be the null or the empty string");
        this.ID = ID.equals(VCFConstants.EMPTY_ID_FIELD) ? VCFConstants.EMPTY_ID_FIELD : ID;

        this.commonInfo = new CommonInfo(source, log10PError, filters, attributes);

        if ( alleles == null ) { throw new IllegalArgumentException("Alleles cannot be null"); }

        // we need to make this a LinkedHashSet in case the user prefers a given ordering of alleles
        this.alleles = makeAlleles(alleles);

        if ( genotypes == null || genotypes == NO_GENOTYPES ) {
            this.genotypes = NO_GENOTYPES;
        } else {
            this.genotypes = genotypes.immutable();
        }

        // cache the REF and ALT alleles
        int nAlleles = alleles.size();
        for ( Allele a : alleles ) {
            if ( a.isReference() ) {
                REF = a;
            } else if ( nAlleles == 2 ) { // only cache ALT when biallelic
                ALT = a;
            }
        }

        this.fullyDecoded = fullyDecoded;

        if ( ! validationToPerform.isEmpty() ) {
            validate(validationToPerform);
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Selectors
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * This method subsets down to a set of samples.
     *
     * At the same time returns the alleles to just those in use by the samples,
     * if rederiveAllelesFromGenotypes is true, otherwise the full set of alleles
     * in this VC is returned as the set of alleles in the subContext, even if
     * some of those alleles aren't in the samples
     *
     * WARNING: BE CAREFUL WITH rederiveAllelesFromGenotypes UNLESS YOU KNOW WHAT YOU ARE DOING?
     *
     * @param sampleNames    the sample names
     * @param rederiveAllelesFromGenotypes if true, returns the alleles to just those in use by the samples, true should be default
     * @return new VariantContext subsetting to just the given samples
     */
    public VariantContext subContextFromSamples(Set<String> sampleNames, final boolean rederiveAllelesFromGenotypes ) {
        if ( sampleNames.containsAll(getSampleNames()) && ! rederiveAllelesFromGenotypes ) {
            return this; // fast path when you don't have any work to do
        } else {
            VariantContextBuilder builder = new VariantContextBuilder(this);
            GenotypesContext newGenotypes = genotypes.subsetToSamples(sampleNames);

            if ( rederiveAllelesFromGenotypes )
                builder.alleles(allelesOfGenotypes(newGenotypes));
            else {
                builder.alleles(alleles);
            }

            return builder.genotypes(newGenotypes).make();
        }
    }

    /**
     * @see #subContextFromSamples(java.util.Set, boolean) with rederiveAllelesFromGenotypes = true
     *
     * @param sampleNames
     * @return
     */
    public VariantContext subContextFromSamples(final Set<String> sampleNames) {
        return subContextFromSamples(sampleNames, true);
    }

    public VariantContext subContextFromSample(String sampleName) {
        return subContextFromSamples(Collections.singleton(sampleName));
    }

    /**
     * helper routine for subcontext
     * @param genotypes genotypes
     * @return allele set
     */
    private final Set<Allele> allelesOfGenotypes(Collection<Genotype> genotypes) {
        final Set<Allele> alleles = new HashSet<Allele>();

        boolean addedref = false;
        for ( final Genotype g : genotypes ) {
            for ( final Allele a : g.getAlleles() ) {
                addedref = addedref || a.isReference();
                if ( a.isCalled() )
                    alleles.add(a);
            }
        }
        if ( ! addedref ) alleles.add(getReference());

        return alleles;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // type operations
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * see: http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=handbook&part=ch5&rendertype=table&id=ch5.ch5_t3
     *
     * Format:
     * dbSNP variation class
     * Rules for assigning allele classes
     * Sample allele definition
     *
     * Single Nucleotide Polymorphisms (SNPs)a
     *   Strictly defined as single base substitutions involving A, T, C, or G.
     *   A/T
     *
     * Deletion/Insertion Polymorphisms (DIPs)
     *   Designated using the full sequence of the insertion as one allele, and either a fully
     *   defined string for the variant allele or a '-' character to specify the deleted allele.
     *   This class will be assigned to a variation if the variation alleles are of different lengths or
     *   if one of the alleles is deleted ('-').
     *   T/-/CCTA/G
     *
     * No-variation
     *   Reports may be submitted for segments of sequence that are assayed and determined to be invariant
     *   in the sample.
     *   (NoVariation)
     *
     * Mixed
     *   Mix of other classes
     *
     * Also supports NO_VARIATION type, used to indicate that the site isn't polymorphic in the population
     *
     *
     * Not currently supported:
     *
     * Heterozygous sequence
     * The term heterozygous is used to specify a region detected by certain methods that do not
     * resolve the polymorphism into a specific sequence motif. In these cases, a unique flanking
     * sequence must be provided to define a sequence context for the variation.
     * (heterozygous)
     *
     * Microsatellite or short tandem repeat (STR)
     * Alleles are designated by providing the repeat motif and the copy number for each allele.
     * Expansion of the allele repeat motif designated in dbSNP into full-length sequence will
     * be only an approximation of the true genomic sequence because many microsatellite markers are
     * not fully sequenced and are resolved as size variants only.
     * (CAC)8/9/10/11
     *
     * Named variant
     * Applies to insertion/deletion polymorphisms of longer sequence features, such as retroposon
     * dimorphism for Alu or line elements. These variations frequently include a deletion '-' indicator
     * for the absent allele.
     * (alu) / -
     *
     * Multi-Nucleotide Polymorphism (MNP)
     *   Assigned to variations that are multi-base variations of a single, common length
     *   GGA/AGT
     */
    public enum Type {
        NO_VARIATION,
        SNP,
        MNP,    // a multi-nucleotide polymorphism
        INDEL,
        SYMBOLIC,
        MIXED,
    }

    /**
     * Determines (if necessary) and returns the type of this variation by examining the alleles it contains.
     *
     * @return the type of this VariantContext
     **/
    public Type getType() {
        if ( type == null )
            determineType();

        return type;
    }

    /**
     * convenience method for SNPs
     *
     * @return true if this is a SNP, false otherwise
     */
    public boolean isSNP() { return getType() == Type.SNP; }


    /**
     * convenience method for variants
     *
     * @return true if this is a variant allele, false if it's reference
     */
    public boolean isVariant() { return getType() != Type.NO_VARIATION; }

    /**
     * convenience method for point events
     *
     * @return true if this is a SNP or ref site, false if it's an indel or mixed event
     */
    public boolean isPointEvent() { return isSNP() || !isVariant(); }

    /**
     * convenience method for indels
     *
     * @return true if this is an indel, false otherwise
     */
    public boolean isIndel() { return getType() == Type.INDEL; }

    /**
     * @return true if the alleles indicate a simple insertion (i.e., the reference allele is Null)
     */
    public boolean isSimpleInsertion() {
        // can't just call !isSimpleDeletion() because of complex indels
        return isSimpleIndel() && getReference().length() == 1;
    }

    /**
     * @return true if the alleles indicate a simple deletion (i.e., a single alt allele that is Null)
     */
    public boolean isSimpleDeletion() {
        // can't just call !isSimpleInsertion() because of complex indels
        return isSimpleIndel() && getAlternateAllele(0).length() == 1;
    }

    /**
     * @return true if the alleles indicate a simple indel, false otherwise.
     */
    public boolean isSimpleIndel() {
        return getType() == Type.INDEL                   // allelic lengths differ
                && isBiallelic()                         // exactly 2 alleles
                && getReference().length() > 0           // ref is not null or symbolic
                && getAlternateAllele(0).length() > 0    // alt is not null or symbolic
                && getReference().getBases()[0] == getAlternateAllele(0).getBases()[0]    // leading bases match for both alleles
                && (getReference().length() == 1 || getAlternateAllele(0).length() == 1);
    }

    /**
     * @return true if the alleles indicate neither a simple deletion nor a simple insertion
     */
    public boolean isComplexIndel() {
        return isIndel() && !isSimpleDeletion() && !isSimpleInsertion();
    }

    public boolean isSymbolic() {
        return getType() == Type.SYMBOLIC;
    }

    public boolean isStructuralIndel() {
        if ( getType() == Type.INDEL ) {
            List<Integer> sizes = getIndelLengths();
            if ( sizes != null ) {
                for ( Integer length : sizes ) {
                    if ( length > MAX_ALLELE_SIZE_FOR_NON_SV ) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     *
     * @return true if the variant is symbolic or a large indel
     */
    public boolean isSymbolicOrSV() {
        return isSymbolic() || isStructuralIndel();
    }

    public boolean isMNP() {
        return getType() == Type.MNP;
    }

    /**
     * convenience method for indels
     *
     * @return true if this is an mixed variation, false otherwise
     */
    public boolean isMixed() { return getType() == Type.MIXED; }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Generic accessors
    //
    // ---------------------------------------------------------------------------------------------------------

    public boolean hasID() {
        return getID() != VCFConstants.EMPTY_ID_FIELD;
    }

    public boolean emptyID() {
        return ! hasID();
    }

    public String getID() {
        return ID;
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // get routines to access context info fields
    //
    // ---------------------------------------------------------------------------------------------------------
    public String getSource()                   { return commonInfo.getName(); }
    public Set<String> getFiltersMaybeNull()    { return commonInfo.getFiltersMaybeNull(); }
    public Set<String> getFilters()             { return commonInfo.getFilters(); }
    public boolean isFiltered()                 { return commonInfo.isFiltered(); }
    public boolean isNotFiltered()              { return commonInfo.isNotFiltered(); }
    public boolean filtersWereApplied()         { return commonInfo.filtersWereApplied(); }
    public boolean hasLog10PError()             { return commonInfo.hasLog10PError(); }
    public double getLog10PError()              { return commonInfo.getLog10PError(); }
    public double getPhredScaledQual()          { return commonInfo.getPhredScaledQual(); }

    public Map<String, Object>  getAttributes() { return commonInfo.getAttributes(); }
    public boolean hasAttribute(String key)     { return commonInfo.hasAttribute(key); }
    public Object getAttribute(String key)      { return commonInfo.getAttribute(key); }

    public Object getAttribute(String key, Object defaultValue) {
        return commonInfo.getAttribute(key, defaultValue);
    }

    public String getAttributeAsString(String key, String defaultValue)   { return commonInfo.getAttributeAsString(key, defaultValue); }
    public int getAttributeAsInt(String key, int defaultValue)            { return commonInfo.getAttributeAsInt(key, defaultValue); }
    public double getAttributeAsDouble(String key, double  defaultValue)  { return commonInfo.getAttributeAsDouble(key, defaultValue); }
    public boolean getAttributeAsBoolean(String key, boolean  defaultValue)  { return commonInfo.getAttributeAsBoolean(key, defaultValue); }

    public CommonInfo getCommonInfo() {
        return commonInfo;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with alleles
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the reference allele for this context
     */
    public Allele getReference() {
        Allele ref = REF;
        if ( ref == null )
            throw new IllegalStateException("BUG: no reference allele found at " + this);
        return ref;
    }


    /**
     * @return true if the context is strictly bi-allelic
     */
    public boolean isBiallelic() {
        return getNAlleles() == 2;
    }

    /**
     * @return The number of segregating alleles in this context
     */
    public int getNAlleles() {
        return alleles.size();
    }

    /**
     * Returns the maximum ploidy of all samples in this VC, or default if there are no genotypes
     *
     * This function is caching, so it's only expensive on the first call
     *
     * @param defaultPloidy the default ploidy, if all samples are no-called
     * @return default, or the max ploidy
     */
    public int getMaxPloidy(final int defaultPloidy) {
        return genotypes.getMaxPloidy(defaultPloidy);
    }

    /**
     * @return The allele sharing the same bases as this String.  A convenience method; better to use byte[]
     */
    public Allele getAllele(String allele) {
        return getAllele(allele.getBytes());
    }

    /**
     * @return The allele sharing the same bases as this byte[], or null if no such allele is present.
     */
    public Allele getAllele(byte[] allele) {
        return Allele.getMatchingAllele(getAlleles(), allele);
    }

    /**
     * @return True if this context contains Allele allele, or false otherwise
     */
    public boolean hasAllele(final Allele allele) {
        return hasAllele(allele, false, true);
    }

    public boolean hasAllele(final Allele allele, final boolean ignoreRefState) {
        return hasAllele(allele, ignoreRefState, true);
    }

    public boolean hasAlternateAllele(final Allele allele) {
        return hasAllele(allele, false, false);
    }

    public boolean hasAlternateAllele(final Allele allele, final boolean ignoreRefState) {
        return hasAllele(allele, ignoreRefState, false);
    }

    private boolean hasAllele(final Allele allele, final boolean ignoreRefState, final boolean considerRefAllele) {
        if ( (considerRefAllele && allele == REF) || allele == ALT ) // optimization for cached cases
            return true;

        final List<Allele> allelesToConsider = considerRefAllele ? getAlleles() : getAlternateAlleles();
        for ( Allele a : allelesToConsider ) {
            if ( a.equals(allele, ignoreRefState) )
                return true;
        }

        return false;
    }


    /**
     * Gets the alleles.  This method should return all of the alleles present at the location,
     * including the reference allele.  There are no constraints imposed on the ordering of alleles
     * in the set. If the reference is not an allele in this context it will not be included.
     *
     * @return the set of alleles
     */
    public List<Allele> getAlleles() { return alleles; }

    /**
     * Gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference allele.  There are no constraints imposed on the ordering of alleles
     * in the set.
     *
     * @return the set of alternate alleles
     */
    public List<Allele> getAlternateAlleles() {
        return alleles.subList(1, alleles.size());
    }

    /**
     * Gets the sizes of the alternate alleles if they are insertion/deletion events, and returns a list of their sizes
     *
     * @return a list of indel lengths ( null if not of type indel or mixed )
     */
    public List<Integer> getIndelLengths() {
        if ( getType() != Type.INDEL && getType() != Type.MIXED ) {
            return null;
        }

        List<Integer> lengths = new ArrayList<Integer>();
        for ( Allele a : getAlternateAlleles() ) {
            lengths.add(a.length() - getReference().length());
        }

        return lengths;
    }

    /**
     * @param i -- the ith allele (from 0 to n - 2 for a context with n alleles including a reference allele)
     * @return the ith non-reference allele in this context
     * @throws IllegalArgumentException if i is invalid
     */
    public Allele getAlternateAllele(int i) {
        return alleles.get(i+1);
    }

    /**
     * @param  other  VariantContext whose alleles to compare against
     * @return true if this VariantContext has the same alleles (both ref and alts) as other,
     *         regardless of ordering. Otherwise returns false.
     */
    public boolean hasSameAllelesAs ( final VariantContext other ) {
        return hasSameAlternateAllelesAs(other) && other.getReference().equals(getReference(), false);
    }

    /**
     * @param  other  VariantContext whose alternate alleles to compare against
     * @return true if this VariantContext has the same alternate alleles as other,
     *         regardless of ordering. Otherwise returns false.
     */
    public boolean hasSameAlternateAllelesAs ( final VariantContext other ) {
        List<Allele> thisAlternateAlleles = getAlternateAlleles();
        List<Allele> otherAlternateAlleles = other.getAlternateAlleles();

        if ( thisAlternateAlleles.size() != otherAlternateAlleles.size() ) {
            return false;
        }

        for ( Allele allele : thisAlternateAlleles ) {
            if ( ! otherAlternateAlleles.contains(allele) ) {
                return false;
            }
        }

        return true;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with genotypes
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the number of samples in the context
     */
    public int getNSamples() {
        return genotypes.size();
    }

    /**
     * @return true if the context has associated genotypes
     */
    public boolean hasGenotypes() {
        return ! genotypes.isEmpty();
    }

    public boolean hasGenotypes(Collection<String> sampleNames) {
        return genotypes.containsSamples(sampleNames);
    }

    /**
     * @return set of all Genotypes associated with this context
     */
    public GenotypesContext getGenotypes() {
        return genotypes;
    }

    public Iterable<Genotype> getGenotypesOrderedByName() {
        return genotypes.iterateInSampleNameOrder();
    }

    public Iterable<Genotype> getGenotypesOrderedBy(Iterable<String> sampleOrdering) {
        return genotypes.iterateInSampleNameOrder(sampleOrdering);
    }

    /**
     * Returns a map from sampleName -> Genotype for the genotype associated with sampleName.  Returns a map
     * for consistency with the multi-get function.
     *
     * @param sampleName   the sample name
     * @return mapping from sample name to genotype
     * @throws IllegalArgumentException if sampleName isn't bound to a genotype
     */
    public GenotypesContext getGenotypes(String sampleName) {
        return getGenotypes(Collections.singleton(sampleName));
    }

    /**
     * Returns a map from sampleName -> Genotype for each sampleName in sampleNames.  Returns a map
     * for consistency with the multi-get function.
     *
     * For testing convenience only
     *
     * @param sampleNames a unique list of sample names
     * @return subsetting genotypes context
     * @throws IllegalArgumentException if sampleName isn't bound to a genotype
     */
    protected GenotypesContext getGenotypes(Collection<String> sampleNames) {
        return getGenotypes().subsetToSamples(new HashSet<String>(sampleNames));
    }

    public GenotypesContext getGenotypes(Set<String> sampleNames) {
        return getGenotypes().subsetToSamples(sampleNames);
    }


    /**
     * @return the set of all sample names in this context, not ordered
     */
    public Set<String> getSampleNames() {
        return getGenotypes().getSampleNames();
    }

    public List<String> getSampleNamesOrderedByName() {
        return getGenotypes().getSampleNamesOrderedByName();
    }

    /**
     * @param sample  the sample name
     *
     * @return the Genotype associated with the given sample in this context or null if the sample is not in this context
     */
    public Genotype getGenotype(String sample) {
        return getGenotypes().get(sample);
    }

    public boolean hasGenotype(String sample) {
        return getGenotypes().containsSample(sample);
    }

    public Genotype getGenotype(int ith) {
        return genotypes.get(ith);
    }


    /**
     * Returns the number of chromosomes carrying any allele in the genotypes (i.e., excluding NO_CALLS)
     *
     * @return chromosome count
     */
    public int getCalledChrCount() {
        final Set<String> noSamples = Collections.emptySet();
        return  getCalledChrCount(noSamples);
    }

    /**
     * Returns the number of chromosomes carrying any allele in the genotypes (i.e., excluding NO_CALLS)
     *
     * @param sampleIds IDs of samples to take into account. If empty then all samples are included.
     * @return chromosome count
     */
    public int getCalledChrCount(Set<String> sampleIds) {
        int n = 0;
        GenotypesContext genotypes = sampleIds.isEmpty() ? getGenotypes() : getGenotypes(sampleIds);

        for ( final Genotype g : genotypes) {
            for ( final Allele a : g.getAlleles() )
                n += a.isNoCall() ? 0 : 1;
        }

        return n;
    }

    /**
     * Returns the number of chromosomes carrying allele A in the genotypes
     *
     * @param a allele
     * @return chromosome count
     */
    public int getCalledChrCount(Allele a) {
        return getCalledChrCount(a,new HashSet<String>(0));
    }

    /**
     * Returns the number of chromosomes carrying allele A in the genotypes
     *
     * @param a allele
     * @param sampleIds - IDs of samples to take into account. If empty then all samples are included.
     * @return chromosome count
     */
    public int getCalledChrCount(Allele a, Set<String> sampleIds) {
        int n = 0;
        GenotypesContext genotypes = sampleIds.isEmpty() ? getGenotypes() : getGenotypes(sampleIds);

        for ( final Genotype g : genotypes ) {
            n += g.countAllele(a);
        }

        return n;
    }

    /**
     * Genotype-specific functions -- are the genotypes monomorphic w.r.t. to the alleles segregating at this
     * site?  That is, is the number of alternate alleles among all fo the genotype == 0?
     *
     * @return true if it's monomorphic
     */
    public boolean isMonomorphicInSamples() {
        if ( monomorphic == null )
            monomorphic = ! isVariant() || (hasGenotypes() && getCalledChrCount(getReference()) == getCalledChrCount());
        return monomorphic;
    }

    /**
     * Genotype-specific functions -- are the genotypes polymorphic w.r.t. to the alleles segregating at this
     * site?  That is, is the number of alternate alleles among all fo the genotype > 0?
     *
     * @return true if it's polymorphic
     */
    public boolean isPolymorphicInSamples() {
        return ! isMonomorphicInSamples();
    }

    private void calculateGenotypeCounts() {
        if ( genotypeCounts == null ) {
            genotypeCounts = new int[GenotypeType.values().length];

            for ( final Genotype g : getGenotypes() ) {
                genotypeCounts[g.getType().ordinal()]++;
            }
        }
    }

    /**
     * Genotype-specific functions -- how many no-calls are there in the genotypes?
     *
     * @return number of no calls
     */
    public int getNoCallCount() {
        calculateGenotypeCounts();
        return genotypeCounts[GenotypeType.NO_CALL.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many hom ref calls are there in the genotypes?
     *
     * @return number of hom ref calls
     */
    public int getHomRefCount() {
        calculateGenotypeCounts();
        return genotypeCounts[GenotypeType.HOM_REF.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many het calls are there in the genotypes?
     *
     * @return number of het calls
     */
    public int getHetCount() {
        calculateGenotypeCounts();
        return genotypeCounts[GenotypeType.HET.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many hom var calls are there in the genotypes?
     *
     * @return number of hom var calls
     */
    public int getHomVarCount() {
        calculateGenotypeCounts();
        return genotypeCounts[GenotypeType.HOM_VAR.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many mixed calls are there in the genotypes?
     *
     * @return number of mixed calls
     */
    public int getMixedCount() {
        calculateGenotypeCounts();
        return genotypeCounts[GenotypeType.MIXED.ordinal()];
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // validation: extra-strict validation routines for paranoid users
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Run all extra-strict validation tests on a Variant Context object
     *
     * @param reportedReference   the reported reference allele
     * @param observedReference   the actual reference allele
     * @param rsIDs               the true dbSNP IDs
     */
    public void extraStrictValidation(final Allele reportedReference, final Allele observedReference, final Set<String> rsIDs) {
        // validate the reference
        validateReferenceBases(reportedReference, observedReference);

        // validate the RS IDs
        validateRSIDs(rsIDs);

        // validate the altenate alleles
        validateAlternateAlleles();

        // validate the AN and AC fields
        validateChromosomeCounts();

        // TODO: implement me
        //checkReferenceTrack();
    }

    public void validateReferenceBases(final Allele reportedReference, final Allele observedReference) {
        if ( reportedReference != null && !reportedReference.basesMatch(observedReference) ) {
            throw new TribbleException.InternalCodecException(String.format("the REF allele is incorrect for the record at position %s:%d, fasta says %s vs. VCF says %s", getChr(), getStart(), observedReference.getBaseString(), reportedReference.getBaseString()));
        }
    }

    public void validateRSIDs(Set<String> rsIDs) {
        if ( rsIDs != null && hasID() ) {
            for ( String id : getID().split(VCFConstants.ID_FIELD_SEPARATOR) ) {
                if ( id.startsWith("rs") && !rsIDs.contains(id) )
                    throw new TribbleException.InternalCodecException(String.format("the rsID %s for the record at position %s:%d is not in dbSNP", id, getChr(), getStart()));
            }
        }
    }

    public void validateAlternateAlleles() {
        if ( !hasGenotypes() )
            return;

        List<Allele> reportedAlleles = getAlleles();
        Set<Allele> observedAlleles = new HashSet<Allele>();
        observedAlleles.add(getReference());
        for ( final Genotype g : getGenotypes() ) {
            if ( g.isCalled() )
                observedAlleles.addAll(g.getAlleles());
        }
        if ( observedAlleles.contains(Allele.NO_CALL) )
            observedAlleles.remove(Allele.NO_CALL);

        if ( reportedAlleles.size() != observedAlleles.size() )
            throw new TribbleException.InternalCodecException(String.format("one or more of the ALT allele(s) for the record at position %s:%d are not observed at all in the sample genotypes", getChr(), getStart()));

        int originalSize = reportedAlleles.size();
        // take the intersection and see if things change
        observedAlleles.retainAll(reportedAlleles);
        if ( observedAlleles.size() != originalSize )
            throw new TribbleException.InternalCodecException(String.format("one or more of the ALT allele(s) for the record at position %s:%d are not observed at all in the sample genotypes", getChr(), getStart()));
    }

    public void validateChromosomeCounts() {
        if ( !hasGenotypes() )
            return;

        // AN
        if ( hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) ) {
            int reportedAN = Integer.valueOf(getAttribute(VCFConstants.ALLELE_NUMBER_KEY).toString());
            int observedAN = getCalledChrCount();
            if ( reportedAN != observedAN )
                throw new TribbleException.InternalCodecException(String.format("the Allele Number (AN) tag is incorrect for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedAN, observedAN));
        }

        // AC
        if ( hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
            ArrayList<Integer> observedACs = new ArrayList<Integer>();

            // if there are alternate alleles, record the relevant tags
            if ( getAlternateAlleles().size() > 0 ) {
                for ( Allele allele : getAlternateAlleles() ) {
                    observedACs.add(getCalledChrCount(allele));
                }
            }
            else { // otherwise, set them to 0
                observedACs.add(0);
            }

            if ( getAttribute(VCFConstants.ALLELE_COUNT_KEY) instanceof List ) {
                Collections.sort(observedACs);
                List reportedACs = (List)getAttribute(VCFConstants.ALLELE_COUNT_KEY);
                Collections.sort(reportedACs);
                if ( observedACs.size() != reportedACs.size() )
                    throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag doesn't have the correct number of values for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedACs.size(), observedACs.size()));
                for (int i = 0; i < observedACs.size(); i++) {
                    if ( Integer.valueOf(reportedACs.get(i).toString()) != observedACs.get(i) )
                        throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag is incorrect for the record at position %s:%d, %s vs. %d", getChr(), getStart(), reportedACs.get(i), observedACs.get(i)));
                }
            } else {
                if ( observedACs.size() != 1 )
                    throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag doesn't have enough values for the record at position %s:%d", getChr(), getStart()));
                int reportedAC = Integer.valueOf(getAttribute(VCFConstants.ALLELE_COUNT_KEY).toString());
                if ( reportedAC != observedACs.get(0) )
                    throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag is incorrect for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedAC, observedACs.get(0)));
            }
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // validation: the normal validation routines are called automatically upon creation of the VC
    //
    // ---------------------------------------------------------------------------------------------------------

    private boolean validate(final EnumSet<Validation> validationToPerform) {
        validateStop();
        for (final Validation val : validationToPerform ) {
            switch (val) {
                case ALLELES: validateAlleles(); break;
                case GENOTYPES: validateGenotypes(); break;
                default: throw new IllegalArgumentException("Unexpected validation mode " + val);
            }
        }

        return true;
    }

    /**
     * Check that getEnd() == END from the info field, if it's present
     */
    private void validateStop() {
        if ( hasAttribute(VCFConstants.END_KEY) ) {
            final int end = getAttributeAsInt(VCFConstants.END_KEY, -1);
            assert end != -1;
            if ( end != getEnd() ) {
                final String message = "Badly formed variant context at location " + getChr() + ":"
                        + getStart() + "; getEnd() was " + getEnd()
                        + " but this VariantContext contains an END key with value " + end;
                if ( GeneralUtils.DEBUG_MODE_ENABLED && WARN_ABOUT_BAD_END ) {
                    System.err.println(message);
                }
                else {
                    throw new TribbleException(message);
                }
            }
        } else {
            final long length = (stop - start) + 1;
            if ( ! hasSymbolicAlleles() && length != getReference().length() ) {
                throw new IllegalStateException("BUG: GenomeLoc " + contig + ":" + start + "-" + stop + " has a size == " + length + " but the variation reference allele has length " + getReference().length() + " this = " + this);
            }
        }
    }

    private void validateAlleles() {

        boolean alreadySeenRef = false;

        for ( final Allele allele : alleles ) {
            // make sure there's only one reference allele
            if ( allele.isReference() ) {
                if ( alreadySeenRef ) throw new IllegalArgumentException("BUG: Received two reference tagged alleles in VariantContext " + alleles + " this=" + this);
                alreadySeenRef = true;
            }

            if ( allele.isNoCall() ) {
                throw new IllegalArgumentException("BUG: Cannot add a no call allele to a variant context " + alleles + " this=" + this);
            }
        }

        // make sure there's one reference allele
        if ( ! alreadySeenRef )
            throw new IllegalArgumentException("No reference allele found in VariantContext");
    }

    private void validateGenotypes() {
        if ( this.genotypes == null ) throw new IllegalStateException("Genotypes is null");

        for ( final Genotype g : this.genotypes ) {
            if ( g.isAvailable() ) {
                for ( Allele gAllele : g.getAlleles() ) {
                    if ( ! hasAllele(gAllele) && gAllele.isCalled() )
                        throw new IllegalStateException("Allele in genotype " + gAllele + " not in the variant context " + alleles);
                }
            }
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // utility routines
    //
    // ---------------------------------------------------------------------------------------------------------

    private void determineType() {
        if ( type == null ) {
            switch ( getNAlleles() ) {
                case 0:
                    throw new IllegalStateException("Unexpected error: requested type of VariantContext with no alleles!" + this);
                case 1:
                    // note that this doesn't require a reference allele.  You can be monomorphic independent of having a
                    // reference allele
                    type = Type.NO_VARIATION;
                    break;
                default:
                    determinePolymorphicType();
            }
        }
    }

    private void determinePolymorphicType() {
        type = null;

        // do a pairwise comparison of all alleles against the reference allele
        for ( Allele allele : alleles ) {
            if ( allele == REF )
                continue;

            // find the type of this allele relative to the reference
            Type biallelicType = typeOfBiallelicVariant(REF, allele);

            // for the first alternate allele, set the type to be that one
            if ( type == null ) {
                type = biallelicType;
            }
            // if the type of this allele is different from that of a previous one, assign it the MIXED type and quit
            else if ( biallelicType != type ) {
                type = Type.MIXED;
                return;
            }
        }
    }

    private static Type typeOfBiallelicVariant(Allele ref, Allele allele) {
        if ( ref.isSymbolic() )
            throw new IllegalStateException("Unexpected error: encountered a record with a symbolic reference allele");

        if ( allele.isSymbolic() )
            return Type.SYMBOLIC;

        if ( ref.length() == allele.length() ) {
            if ( allele.length() == 1 )
                return Type.SNP;
            else
                return Type.MNP;
        }

        // Important note: previously we were checking that one allele is the prefix of the other.  However, that's not an
        // appropriate check as can be seen from the following example:
        // REF = CTTA and ALT = C,CT,CA
        // This should be assigned the INDEL type but was being marked as a MIXED type because of the prefix check.
        // In truth, it should be absolutely impossible to return a MIXED type from this method because it simply
        // performs a pairwise comparison of a single alternate allele against the reference allele (whereas the MIXED type
        // is reserved for cases of multiple alternate alleles of different types).  Therefore, if we've reached this point
        // in the code (so we're not a SNP, MNP, or symbolic allele), we absolutely must be an INDEL.

        return Type.INDEL;

        // old incorrect logic:
        // if (oneIsPrefixOfOther(ref, allele))
        //     return Type.INDEL;
        // else
        //     return Type.MIXED;
    }

    public String toString() {
        return String.format("[VC %s @ %s Q%s of type=%s alleles=%s attr=%s GT=%s",
                getSource(), contig + ":" + (start - stop == 0 ? start : start + "-" + stop),
                hasLog10PError() ? String.format("%.2f", getPhredScaledQual()) : ".",
                this.getType(),
                ParsingUtils.sortList(this.getAlleles()),
                ParsingUtils.sortedString(this.getAttributes()),
                this.getGenotypes());
    }

    public String toStringWithoutGenotypes() {
        return String.format("[VC %s @ %s Q%s of type=%s alleles=%s attr=%s",
                getSource(), contig + ":" + (start - stop == 0 ? start : start + "-" + stop),
                hasLog10PError() ? String.format("%.2f", getPhredScaledQual()) : ".",
                this.getType(),
                ParsingUtils.sortList(this.getAlleles()),
                ParsingUtils.sortedString(this.getAttributes()));
    }

    // protected basic manipulation routines
    private static List<Allele> makeAlleles(Collection<Allele> alleles) {
        final List<Allele> alleleList = new ArrayList<Allele>(alleles.size());

        boolean sawRef = false;
        for ( final Allele a : alleles ) {
            for ( final Allele b : alleleList ) {
                if ( a.equals(b, true) )
                    throw new IllegalArgumentException("Duplicate allele added to VariantContext: " + a);
            }

            // deal with the case where the first allele isn't the reference
            if ( a.isReference() ) {
                if ( sawRef )
                    throw new IllegalArgumentException("Alleles for a VariantContext must contain at most one reference allele: " + alleles);
                alleleList.add(0, a);
                sawRef = true;
            }
            else
                alleleList.add(a);
        }

        if ( alleleList.isEmpty() )
            throw new IllegalArgumentException("Cannot create a VariantContext with an empty allele list");

        if ( alleleList.get(0).isNonReference() )
            throw new IllegalArgumentException("Alleles for a VariantContext must contain at least one reference allele: " + alleles);

        return alleleList;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Fully decode
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Return a VC equivalent to this one but where all fields are fully decoded
     *
     * See VariantContext document about fully decoded
     *
     * @param header containing types about all fields in this VC
     * @return a fully decoded version of this VC
     */
    public VariantContext fullyDecode(final VCFHeader header, final boolean lenientDecoding) {
        if ( isFullyDecoded() )
            return this;
        else {
            // TODO -- warning this is potentially very expensive as it creates copies over and over
            final VariantContextBuilder builder = new VariantContextBuilder(this);
            fullyDecodeInfo(builder, header, lenientDecoding);
            fullyDecodeGenotypes(builder, header);
            builder.fullyDecoded(true);
            return builder.make();
        }
    }

    /**
     * See VariantContext document about fully decoded
     * @return true if this is a fully decoded VC
     */
    public boolean isFullyDecoded() {
        return fullyDecoded;
    }

    private final void fullyDecodeInfo(final VariantContextBuilder builder, final VCFHeader header, final boolean lenientDecoding) {
        builder.attributes(fullyDecodeAttributes(getAttributes(), header, lenientDecoding));
    }

    private final Map<String, Object> fullyDecodeAttributes(final Map<String, Object> attributes,
                                                            final VCFHeader header,
                                                            final boolean lenientDecoding) {
        final Map<String, Object> newAttributes = new HashMap<String, Object>(10);

        for ( final Map.Entry<String, Object> attr : attributes.entrySet() ) {
            final String field = attr.getKey();

            if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY) )
                continue; // gross, FT is part of the extended attributes

            final VCFCompoundHeaderLine format = VariantContextUtils.getMetaDataForField(header, field);
            final Object decoded = decodeValue(field, attr.getValue(), format);

            if ( decoded != null &&
                    ! lenientDecoding
                    && format.getCountType() != VCFHeaderLineCount.UNBOUNDED
                    && format.getType() != VCFHeaderLineType.Flag ) { // we expect exactly the right number of elements
                final int obsSize = decoded instanceof List ? ((List) decoded).size() : 1;
                final int expSize = format.getCount(this);
                if ( obsSize != expSize ) {
                    throw new TribbleException.InvalidHeader("Discordant field size detected for field " +
                            field + " at " + getChr() + ":" + getStart() + ".  Field had " + obsSize + " values " +
                            "but the header says this should have " + expSize + " values based on header record " +
                            format);
                }
            }
            newAttributes.put(field, decoded);
        }

        return newAttributes;
    }

    private final Object decodeValue(final String field, final Object value, final VCFCompoundHeaderLine format) {
        if ( value instanceof String ) {
            if ( field.equals(VCFConstants.GENOTYPE_PL_KEY) )
                return GenotypeLikelihoods.fromPLField((String)value);

            final String string = (String)value;
            if ( string.indexOf(",") != -1 ) {
                final String[] splits = string.split(",");
                final List<Object> values = new ArrayList<Object>(splits.length);
                for ( int i = 0; i < splits.length; i++ )
                    values.add(decodeOne(field, splits[i], format));
                return values;
            } else {
                return decodeOne(field, string, format);
            }
        } else if ( value instanceof List && (((List) value).get(0)) instanceof String ) {
            final List<String> asList = (List<String>)value;
            final List<Object> values = new ArrayList<Object>(asList.size());
            for ( final String s : asList )
                values.add(decodeOne(field, s, format));
            return values;
        } else {
            return value;
        }

        // allowMissingValuesComparedToHeader
    }

    private final Object decodeOne(final String field, final String string, final VCFCompoundHeaderLine format) {
        try {
            if ( string.equals(VCFConstants.MISSING_VALUE_v4) )
                return null;
            else {
                switch ( format.getType() ) {
                    case Character: return string;
                    case Flag:
                        final boolean b = Boolean.valueOf(string) || string.equals("1");
                        if ( b == false )
                            throw new TribbleException("VariantContext FLAG fields " + field + " cannot contain false values"
                             + " as seen at " + getChr() + ":" + getStart());
                        return b;
                    case String:    return string;
                    case Integer:   return Integer.valueOf(string);
                    case Float:     return Double.valueOf(string);
                    default: throw new TribbleException("Unexpected type for field" + field);
                }
            }
        } catch (NumberFormatException e) {
            throw new TribbleException("Could not decode field " + field + " with value " + string + " of declared type " + format.getType());
        }
    }

    private final void fullyDecodeGenotypes(final VariantContextBuilder builder, final VCFHeader header) {
        final GenotypesContext gc = new GenotypesContext();
        for ( final Genotype g : getGenotypes() ) {
            gc.add(fullyDecodeGenotypes(g, header));
        }
        builder.genotypesNoValidation(gc);
    }

    private final Genotype fullyDecodeGenotypes(final Genotype g, final VCFHeader header) {
        final Map<String, Object> map = fullyDecodeAttributes(g.getExtendedAttributes(), header, true);
        return new GenotypeBuilder(g).attributes(map).make();
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // tribble integration routines -- not for public consumption
    //
    // ---------------------------------------------------------------------------------------------------------
    public String getChr() {
        return contig;
    }

    public int getStart() {
        return (int)start;
    }

    public int getEnd() {
        return (int)stop;
    }

    public boolean hasSymbolicAlleles() {
        return hasSymbolicAlleles(getAlleles());
    }

    public static boolean hasSymbolicAlleles( final List<Allele> alleles ) {
        for ( final Allele a: alleles ) {
            if (a.isSymbolic()) {
                return true;
            }
        }
        return false;
    }

    public Allele getAltAlleleWithHighestAlleleCount() {
        // optimization: for bi-allelic sites, just return the 1only alt allele
        if ( isBiallelic() )
            return getAlternateAllele(0);

        Allele best = null;
        int maxAC1 = 0;
        for ( Allele a : getAlternateAlleles() ) {
            final int ac = getCalledChrCount(a);
            if ( ac >= maxAC1 ) {
                maxAC1 = ac;
                best = a;
            }

        }
        return best;
    }

    /**
     * Lookup the index of allele in this variant context
     *
     * @param allele the allele whose index we want to get
     * @return the index of the allele into getAlleles(), or -1 if it cannot be found
     */
    public int getAlleleIndex(final Allele allele) {
        return getAlleles().indexOf(allele);
    }

    /**
     * Return the allele index #getAlleleIndex for each allele in alleles
     *
     * @param alleles the alleles we want to look up
     * @return a list of indices for each allele, in order
     */
    public List<Integer> getAlleleIndices(final Collection<Allele> alleles) {
        final List<Integer> indices = new LinkedList<Integer>();
        for ( final Allele allele : alleles )
            indices.add(getAlleleIndex(allele));
        return indices;
    }

    public int[] getGLIndecesOfAlternateAllele(Allele targetAllele) {
        final int index = getAlleleIndex(targetAllele);
        if ( index == -1 ) throw new IllegalArgumentException("Allele " + targetAllele + " not in this VariantContex " + this);
        return GenotypeLikelihoods.getPLIndecesOfAlleles(0, index);
    }
}
