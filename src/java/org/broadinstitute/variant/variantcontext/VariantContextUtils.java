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
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlEngine;
import org.broad.tribble.TribbleException;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.vcf.*;

import java.util.*;

public class VariantContextUtils {

    private static Set<String> MISSING_KEYS_WARNED_ABOUT = new HashSet<String>();

    final public static JexlEngine engine = new JexlEngine();
    private final static boolean ASSUME_MISSING_FIELDS_ARE_STRINGS = false;

    static {
        engine.setSilent(false); // will throw errors now for selects that don't evaluate properly
        engine.setLenient(false);
        engine.setDebug(false);
    }

    /**
     * Update the attributes of the attributes map given the VariantContext to reflect the
     * proper chromosome-based VCF tags
     *
     * @param vc          the VariantContext
     * @param attributes  the attributes map to populate; must not be null; may contain old values
     * @param removeStaleValues should we remove stale values from the mapping?
     * @return the attributes map provided as input, returned for programming convenience
     */
    public static Map<String, Object> calculateChromosomeCounts(VariantContext vc, Map<String, Object> attributes, boolean removeStaleValues) {
        return calculateChromosomeCounts(vc, attributes,  removeStaleValues, new HashSet<String>(0));
    }

    /**
     * Update the attributes of the attributes map given the VariantContext to reflect the
     * proper chromosome-based VCF tags
     *
     * @param vc          the VariantContext
     * @param attributes  the attributes map to populate; must not be null; may contain old values
     * @param removeStaleValues should we remove stale values from the mapping?
     * @param founderIds - Set of founders Ids to take into account. AF and FC will be calculated over the founders.
     *                  If empty or null, counts are generated for all samples as unrelated individuals
     * @return the attributes map provided as input, returned for programming convenience
     */
    public static Map<String, Object> calculateChromosomeCounts(VariantContext vc, Map<String, Object> attributes, boolean removeStaleValues, final Set<String> founderIds) {
        final int AN = vc.getCalledChrCount();

        // if everyone is a no-call, remove the old attributes if requested
        if ( AN == 0 && removeStaleValues ) {
            if ( attributes.containsKey(VCFConstants.ALLELE_COUNT_KEY) )
                attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
            if ( attributes.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY) )
                attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
            if ( attributes.containsKey(VCFConstants.ALLELE_NUMBER_KEY) )
                attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
            return attributes;
        }

        if ( vc.hasGenotypes() ) {
            attributes.put(VCFConstants.ALLELE_NUMBER_KEY, AN);

            // if there are alternate alleles, record the relevant tags
            if ( vc.getAlternateAlleles().size() > 0 ) {
                ArrayList<Double> alleleFreqs = new ArrayList<Double>();
                ArrayList<Integer> alleleCounts = new ArrayList<Integer>();
                ArrayList<Integer> foundersAlleleCounts = new ArrayList<Integer>();
                double totalFoundersChromosomes = (double)vc.getCalledChrCount(founderIds);
                int foundersAltChromosomes;
                for ( Allele allele : vc.getAlternateAlleles() ) {
                    foundersAltChromosomes = vc.getCalledChrCount(allele,founderIds);
                    alleleCounts.add(vc.getCalledChrCount(allele));
                    foundersAlleleCounts.add(foundersAltChromosomes);
                    if ( AN == 0 ) {
                        alleleFreqs.add(0.0);
                    } else {
                        final Double freq = (double)foundersAltChromosomes / totalFoundersChromosomes;
                        alleleFreqs.add(freq);
                    }
                }

                attributes.put(VCFConstants.ALLELE_COUNT_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
                attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqs.size() == 1 ? alleleFreqs.get(0) : alleleFreqs);
            } else {
                // if there's no alt AC and AF shouldn't be present
                attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
                attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
            }
        }

        return attributes;
    }

    /**
     * Update the attributes of the attributes map in the VariantContextBuilder to reflect the proper
     * chromosome-based VCF tags based on the current VC produced by builder.make()
     *
     * @param builder     the VariantContextBuilder we are updating
     * @param removeStaleValues should we remove stale values from the mapping?
     */
    public static void calculateChromosomeCounts(VariantContextBuilder builder, boolean removeStaleValues) {
        VariantContext vc = builder.make();
        builder.attributes(calculateChromosomeCounts(vc, new HashMap<String, Object>(vc.getAttributes()), removeStaleValues, new HashSet<String>(0)));
    }

    /**
     * Update the attributes of the attributes map in the VariantContextBuilder to reflect the proper
     * chromosome-based VCF tags based on the current VC produced by builder.make()
     *
     * @param builder     the VariantContextBuilder we are updating
     * @param founderIds - Set of founders to take into account. AF and FC will be calculated over the founders only.
     *                   If empty or null, counts are generated for all samples as unrelated individuals
     * @param removeStaleValues should we remove stale values from the mapping?
     */
    public static void calculateChromosomeCounts(VariantContextBuilder builder, boolean removeStaleValues, final Set<String> founderIds) {
        VariantContext vc = builder.make();
        builder.attributes(calculateChromosomeCounts(vc, new HashMap<String, Object>(vc.getAttributes()), removeStaleValues, founderIds));
    }

    public final static VCFCompoundHeaderLine getMetaDataForField(final VCFHeader header, final String field) {
        VCFCompoundHeaderLine metaData = header.getFormatHeaderLine(field);
        if ( metaData == null ) metaData = header.getInfoHeaderLine(field);
        if ( metaData == null ) {
            if ( ASSUME_MISSING_FIELDS_ARE_STRINGS ) {
                if ( ! MISSING_KEYS_WARNED_ABOUT.contains(field) ) {
                    MISSING_KEYS_WARNED_ABOUT.add(field);
                    if ( GeneralUtils.DEBUG_MODE_ENABLED )
                        System.err.println("Field " + field + " missing from VCF header, assuming it is an unbounded string type");
                }
                return new VCFInfoHeaderLine(field, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Auto-generated string header for " + field);
            }
            else
                throw new TribbleException("Fully decoding VariantContext requires header line for all fields, but none was found for " + field);
        }
        return metaData;
    }

    /**
     * A simple but common wrapper for matching VariantContext objects using JEXL expressions
     */
    public static class JexlVCMatchExp {
        public String name;
        public Expression exp;

        /**
         * Create a new matcher expression with name and JEXL expression exp
         * @param name name
         * @param exp  expression
         */
        public JexlVCMatchExp(String name, Expression exp) {
            this.name = name;
            this.exp = exp;
        }
    }

    /**
     * Method for creating JexlVCMatchExp from input walker arguments names and exps.  These two arrays contain
     * the name associated with each JEXL expression. initializeMatchExps will parse each expression and return
     * a list of JexlVCMatchExp, in order, that correspond to the names and exps.  These are suitable input to
     * match() below.
     *
     * @param names names
     * @param exps  expressions
     * @return list of matches
     */
    public static List<JexlVCMatchExp> initializeMatchExps(String[] names, String[] exps) {
        if ( names == null || exps == null )
            throw new IllegalArgumentException("BUG: neither names nor exps can be null: names " + Arrays.toString(names) + " exps=" + Arrays.toString(exps) );

        if ( names.length != exps.length )
            throw new IllegalArgumentException("Inconsistent number of provided filter names and expressions: names=" + Arrays.toString(names) + " exps=" + Arrays.toString(exps));

        Map<String, String> map = new HashMap<String, String>();
        for ( int i = 0; i < names.length; i++ ) { map.put(names[i], exps[i]); }

        return VariantContextUtils.initializeMatchExps(map);
    }

    public static List<JexlVCMatchExp> initializeMatchExps(ArrayList<String> names, ArrayList<String> exps) {
        String[] nameArray = new String[names.size()];
        String[] expArray = new String[exps.size()];
        return initializeMatchExps(names.toArray(nameArray), exps.toArray(expArray));
    }


    /**
     * Method for creating JexlVCMatchExp from input walker arguments mapping from names to exps.  These two arrays contain
     * the name associated with each JEXL expression. initializeMatchExps will parse each expression and return
     * a list of JexlVCMatchExp, in order, that correspond to the names and exps.  These are suitable input to
     * match() below.
     *
     * @param names_and_exps mapping of names to expressions
     * @return list of matches
     */
    public static List<JexlVCMatchExp> initializeMatchExps(Map<String, String> names_and_exps) {
        List<JexlVCMatchExp> exps = new ArrayList<JexlVCMatchExp>();

        for ( Map.Entry<String, String> elt : names_and_exps.entrySet() ) {
            String name = elt.getKey();
            String expStr = elt.getValue();

            if ( name == null || expStr == null ) throw new IllegalArgumentException("Cannot create null expressions : " + name +  " " + expStr);
            try {
                Expression exp = engine.createExpression(expStr);
                exps.add(new JexlVCMatchExp(name, exp));
            } catch (Exception e) {
                throw new IllegalArgumentException("Argument " + name + "has a bad value. Invalid expression used (" + expStr + "). Please see the JEXL docs for correct syntax.") ;
            }
        }

        return exps;
    }

    /**
     * Returns true if exp match VC.  See collection<> version for full docs.
     * @param vc    variant context
     * @param exp   expression
     * @return true if there is a match
     */
    public static boolean match(VariantContext vc, JexlVCMatchExp exp) {
        return match(vc,Arrays.asList(exp)).get(exp);
    }

    /**
     * Matches each JexlVCMatchExp exp against the data contained in vc, and returns a map from these
     * expressions to true (if they matched) or false (if they didn't).  This the best way to apply JEXL
     * expressions to VariantContext records.  Use initializeMatchExps() to create the list of JexlVCMatchExp
     * expressions.
     *
     * @param vc   variant context
     * @param exps expressions
     * @return true if there is a match
     */
    public static Map<JexlVCMatchExp, Boolean> match(VariantContext vc, Collection<JexlVCMatchExp> exps) {
        return new JEXLMap(exps,vc);

    }

    /**
     * Returns true if exp match VC/g.  See collection<> version for full docs.
     * @param vc   variant context
     * @param g    genotype
     * @param exp   expression
     * @return true if there is a match
     */
    public static boolean match(VariantContext vc, Genotype g, JexlVCMatchExp exp) {
        return match(vc,g,Arrays.asList(exp)).get(exp);
    }

    /**
     * Matches each JexlVCMatchExp exp against the data contained in vc/g, and returns a map from these
     * expressions to true (if they matched) or false (if they didn't).  This the best way to apply JEXL
     * expressions to VariantContext records/genotypes.  Use initializeMatchExps() to create the list of JexlVCMatchExp
     * expressions.
     *
     * @param vc   variant context
     * @param g    genotype
     * @param exps expressions
     * @return true if there is a match
     */
    public static Map<JexlVCMatchExp, Boolean> match(VariantContext vc, Genotype g, Collection<JexlVCMatchExp> exps) {
        return new JEXLMap(exps,vc,g);
    }

    /**
     * Returns a newly allocated VC that is the same as VC, but without genotypes
     * @param vc  variant context
     * @return  new VC without genotypes
     */
    @Requires("vc != null")
    @Ensures("result != null")
    public static VariantContext sitesOnlyVariantContext(VariantContext vc) {
        return new VariantContextBuilder(vc).noGenotypes().make();
    }

    /**
     * Returns a newly allocated list of VC, where each VC is the same as the input VCs, but without genotypes
     * @param vcs  collection of VCs
     * @return new VCs without genotypes
     */
    @Requires("vcs != null")
    @Ensures("result != null")
    public static Collection<VariantContext> sitesOnlyVariantContexts(Collection<VariantContext> vcs) {
        List<VariantContext> r = new ArrayList<VariantContext>();
        for ( VariantContext vc : vcs )
            r.add(sitesOnlyVariantContext(vc));
        return r;
    }

    //    TODO: remove that after testing
//    static private void verifyUniqueSampleNames(Collection<VariantContext> unsortedVCs) {
//        Set<String> names = new HashSet<String>();
//        for ( VariantContext vc : unsortedVCs ) {
//            for ( String name : vc.getSampleNames() ) {
//                //System.out.printf("Checking %s %b%n", name, names.contains(name));
//                if ( names.contains(name) )
//                    throw new IllegalStateException("REQUIRE_UNIQUE sample names is true but duplicate names were discovered " + name);
//            }
//
//            names.addAll(vc.getSampleNames());
//        }
//    }


    public static int getSize( VariantContext vc ) {
        return vc.getEnd() - vc.getStart() + 1;
    }

    public static final Set<String> genotypeNames(final Collection<Genotype> genotypes) {
        final Set<String> names = new HashSet<String>(genotypes.size());
        for ( final Genotype g : genotypes )
            names.add(g.getSampleName());
        return names;
    }

    /**
     * Compute the end position for this VariantContext from the alleles themselves
     *
     * In the trivial case this is a single BP event and end = start (open intervals)
     * In general the end is start + ref length - 1, handling the case where ref length == 0
     * However, if alleles contains a symbolic allele then we use endForSymbolicAllele in all cases
     *
     * @param alleles the list of alleles to consider.  The reference allele must be the first one
     * @param start the known start position of this event
     * @param endForSymbolicAlleles the end position to use if any of the alleles is symbolic.  Can be -1
     *                              if no is expected but will throw an error if one is found
     * @return this builder
     */
    @Requires({"! alleles.isEmpty()", "start > 0", "endForSymbolicAlleles == -1 || endForSymbolicAlleles > 0" })
    public static int computeEndFromAlleles(final List<Allele> alleles, final int start, final int endForSymbolicAlleles) {
        final Allele ref = alleles.get(0);

        if ( ref.isNonReference() )
            throw new IllegalStateException("computeEndFromAlleles requires first allele to be reference");

        if ( VariantContext.hasSymbolicAlleles(alleles) ) {
            if ( endForSymbolicAlleles == -1 )
                throw new IllegalStateException("computeEndFromAlleles found a symbolic allele but endForSymbolicAlleles was provided");
            return endForSymbolicAlleles;
        } else {
            return start + Math.max(ref.length() - 1, 0);
        }
    }

}
