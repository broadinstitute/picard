package org.broadinstitute.variant.variantcontext;

import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import org.broadinstitute.variant.vcf.VCFConstants;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * this is an implementation of a Map of JexlVCMatchExp to true or false values.  It lazy initializes each value
 * as requested to save as much processing time as possible.
 *
 * Compatible with JEXL 1.1 (this code will be easier if we move to 2.0, all of the functionality can go into the
 * JexlContext's get()
 *
 */

class JEXLMap implements Map<JexlVCMatchExp, Boolean> {
    // our variant context and/or Genotype
    private final VariantContext vc;
    private final Genotype g;

    // our context
    private JexlContext jContext = null;

    // our mapping from JEXLVCMatchExp to Booleans, which will be set to NULL for previously uncached JexlVCMatchExp
    private Map<JexlVCMatchExp,Boolean> jexl;


    public JEXLMap(Collection<JexlVCMatchExp> jexlCollection, VariantContext vc, Genotype g) {
        this.vc = vc;
        this.g = g;
        initialize(jexlCollection);
    }

    public JEXLMap(Collection<JexlVCMatchExp> jexlCollection, VariantContext vc) {
        this(jexlCollection, vc, null);
    }

    private void initialize(Collection<JexlVCMatchExp> jexlCollection) {
        jexl = new HashMap<JexlVCMatchExp,Boolean>();
        for (JexlVCMatchExp exp: jexlCollection) {
            jexl.put(exp, null);
        }
    }

    /**
     * create the internal JexlContext, only when required.  This code is where new JEXL context variables
     * should get added.
     *
     */
    private void createContext() {
        if ( g == null ) {
            // todo -- remove dependancy on g to the entire system
            jContext = new VariantJEXLContext(vc);
        } else {
            //
            // this whole branch is here just to support G jexl operations
            //
            Map<String, Object> infoMap = new HashMap<String, Object>();

            if ( vc != null ) {
                // create a mapping of what we know about the variant context, its Chromosome, positions, etc.
                infoMap.put("CHROM", vc.getChr());
                infoMap.put("POS", vc.getStart());
                infoMap.put("TYPE", vc.getType().toString());
                infoMap.put("QUAL", String.valueOf(vc.getPhredScaledQual()));

                // add alleles
                infoMap.put("ALLELES", GeneralUtils.join(";", vc.getAlleles()));
                infoMap.put("N_ALLELES", String.valueOf(vc.getNAlleles()));

                // add attributes
                addAttributesToMap(infoMap, vc.getAttributes());

                // add filter fields
                infoMap.put("FILTER", vc.isFiltered() ? "1" : "0");
                for ( Object filterCode : vc.getFilters() ) {
                    infoMap.put(String.valueOf(filterCode), "1");
                }

                // add genotype-specific fields
                // TODO -- implement me when we figure out a good way to represent this
                //    for ( Genotype g : vc.getGenotypes().values() ) {
                //        String prefix = g.getSampleName() + ".";
                //        addAttributesToMap(infoMap, g.getAttributes(), prefix);
                //        infoMap.put(prefix + "GT", g.getGenotypeString());
                //    }

                // add specific genotype if one is provided
                infoMap.put(VCFConstants.GENOTYPE_KEY, g.getGenotypeString());
                infoMap.put("isHomRef", g.isHomRef() ? "1" : "0");
                infoMap.put("isHet", g.isHet() ? "1" : "0");
                infoMap.put("isHomVar", g.isHomVar() ? "1" : "0");
                infoMap.put(VCFConstants.GENOTYPE_QUALITY_KEY, g.getGQ());
                if ( g.hasDP() )
                    infoMap.put(VCFConstants.DEPTH_KEY, g.getDP());
                for ( Entry<String, Object> e : g.getExtendedAttributes().entrySet() ) {
                    if ( e.getValue() != null && !e.getValue().equals(VCFConstants.MISSING_VALUE_v4) )
                        infoMap.put(e.getKey(), e.getValue());
                }
            }

            // create the internal context that we can evaluate expressions against

            jContext = new MapContext(infoMap);
        }
    }

    /**
     * @return the size of the internal data structure
     */
    public int size() {
        return jexl.size();
    }

    /**
     * @return true if we're empty
     */
    public boolean isEmpty() { return this.jexl.isEmpty(); }

    /**
     * do we contain the specified key
     * @param o the key
     * @return true if we have a value for that key
     */
    public boolean containsKey(Object o) { return jexl.containsKey(o); }

    public Boolean get(Object o) {
        // if we've already determined the value, return it
        if (jexl.containsKey(o) && jexl.get(o) != null) return jexl.get(o);

        // try and cast the expression
        JexlVCMatchExp e = (JexlVCMatchExp) o;
        evaluateExpression(e);
        return jexl.get(e);
    }

    /**
     * get the keyset of map
     * @return a set of keys of type JexlVCMatchExp
     */
    public Set<JexlVCMatchExp> keySet() {
        return jexl.keySet();
    }

    /**
     * get all the values of the map.  This is an expensive call, since it evaluates all keys that haven't
     * been evaluated yet.  This is fine if you truely want all the keys, but if you only want a portion, or  know
     * the keys you want, you would be better off using get() to get them by name.
     * @return a collection of boolean values, representing the results of all the variants evaluated
     */
    public Collection<Boolean> values() {
        // this is an expensive call
        for (JexlVCMatchExp exp : jexl.keySet())
            if (jexl.get(exp) == null)
                evaluateExpression(exp);
        return jexl.values();
    }

    /**
     * evaulate a JexlVCMatchExp's expression, given the current context (and setup the context if it's null)
     * @param exp the JexlVCMatchExp to evaluate
     */
    private void evaluateExpression(JexlVCMatchExp exp) {
        // if the context is null, we need to create it to evaluate the JEXL expression
        if (this.jContext == null) createContext();
        try {
            final Boolean value = (Boolean) exp.exp.evaluate(jContext);
            // treat errors as no match
            jexl.put(exp, value == null ? false : value);
        } catch (Exception e) {
            // if exception happens because variable is undefined (i.e. field in expression is not present), evaluate to FALSE
            // todo - might be safer if we explicitly checked for an exception type, but Apache's API doesn't seem to have that ability
            if (e.getMessage().contains("undefined variable"))
                jexl.put(exp,false);
            else
                throw new IllegalArgumentException(String.format("Invalid JEXL expression detected for %s with message %s", exp.name, e.getMessage()));
        }
    }

    /**
     * helper function: adds the list of attributes to the information map we're building
     * @param infoMap the map
     * @param attributes the attributes
     */
    private static void addAttributesToMap(Map<String, Object> infoMap, Map<String, ?> attributes ) {
        for (Entry<String, ?> e : attributes.entrySet()) {
            infoMap.put(e.getKey(), String.valueOf(e.getValue()));
        }
    }

    public Boolean put(JexlVCMatchExp jexlVCMatchExp, Boolean aBoolean) {
        return jexl.put(jexlVCMatchExp,aBoolean);
    }

    public void putAll(Map<? extends JexlVCMatchExp, ? extends Boolean> map) {
        jexl.putAll(map);
    }

    // //////////////////////////////////////////////////////////////////////////////////////
    // The Following are unsupported at the moment
    // //////////////////////////////////////////////////////////////////////////////////////

    // this doesn't make much sense to implement, boolean doesn't offer too much variety to deal
    // with evaluating every key in the internal map.
    public boolean containsValue(Object o) {
        throw new UnsupportedOperationException("containsValue() not supported on a JEXLMap");
    }

    // this doesn't make much sense
    public Boolean remove(Object o) {
        throw new UnsupportedOperationException("remove() not supported on a JEXLMap");
    }


    public Set<Entry<JexlVCMatchExp, Boolean>> entrySet() {
        throw new UnsupportedOperationException("clear() not supported on a JEXLMap");
    }

    // nope
    public void clear() {
        throw new UnsupportedOperationException("clear() not supported on a JEXLMap");
    }
}
