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

import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.vcf.VCFConstants;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author aaron
 * @author depristo
 *
 * Class VariantJEXLContext
 *
 * implements the JEXML context for VariantContext; this saves us from
 * having to generate a JEXML context lookup map everytime we want to evaluate an expression.
 *
 * This is package protected, only classes in variantcontext should have access to it.
 *
 * // todo -- clean up to remove or better support genotype filtering 
 */

class VariantJEXLContext implements JexlContext {
    // our stored variant context
    private VariantContext vc;

    private interface AttributeGetter {
        public Object get(VariantContext vc);
    }

    private static Map<String, AttributeGetter> x = new HashMap<String, AttributeGetter>();

    static {
        x.put("vc",   new AttributeGetter() { public Object get(VariantContext vc) { return vc; }});
        x.put("CHROM",   new AttributeGetter() { public Object get(VariantContext vc) { return vc.getChr(); }});
        x.put("POS",     new AttributeGetter() { public Object get(VariantContext vc) { return vc.getStart(); }});
        x.put("TYPE",    new AttributeGetter() { public Object get(VariantContext vc) { return vc.getType().toString(); }});
        x.put("QUAL",    new AttributeGetter() { public Object get(VariantContext vc) { return -10 * vc.getLog10PError(); }});
        x.put("ALLELES", new AttributeGetter() { public Object get(VariantContext vc) { return vc.getAlleles(); }});
        x.put("N_ALLELES", new AttributeGetter() { public Object get(VariantContext vc) { return vc.getNAlleles(); }});
        x.put("FILTER",    new AttributeGetter() { public Object get(VariantContext vc) { return vc.isFiltered() ? "1" : "0"; }});

//        x.put("GT",        new AttributeGetter() { public Object get(VariantContext vc) { return g.getGenotypeString(); }});
        x.put("homRefCount",  new AttributeGetter() { public Object get(VariantContext vc) { return vc.getHomRefCount(); }});
        x.put("hetCount",     new AttributeGetter() { public Object get(VariantContext vc) { return vc.getHetCount(); }});
        x.put("homVarCount",  new AttributeGetter() { public Object get(VariantContext vc) { return vc.getHomVarCount(); }});
    }

    public VariantJEXLContext(VariantContext vc) {
        this.vc = vc;
    }

    public Object get(String name) {
        Object result = null;
        if ( x.containsKey(name) ) { // dynamic resolution of name -> value via map
            result = x.get(name).get(vc);
        } else if ( vc.hasAttribute(name)) {
            result = vc.getAttribute(name);
        } else if ( vc.getFilters().contains(name) ) {
            result = "1";
        }

        //System.out.printf("dynamic lookup %s => %s%n", name, result);

        return result;
    }

    public boolean has(String name) {
        return get(name) != null;
    }

    public void	set(String name, Object value) {
        throw new UnsupportedOperationException("remove() not supported on a VariantJEXLContext");
    }
}




