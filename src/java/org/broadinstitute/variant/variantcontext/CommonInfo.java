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


import org.broadinstitute.variant.vcf.VCFConstants;

import java.util.*;


/**
 * Common utility routines for VariantContext and Genotype
 *
 * @author depristo
 */
public final class CommonInfo {
    public static final double NO_LOG10_PERROR = 1.0;

    private static Set<String> NO_FILTERS = Collections.emptySet();
    private static Map<String, Object> NO_ATTRIBUTES = Collections.unmodifiableMap(new HashMap<String, Object>());

    private double log10PError = NO_LOG10_PERROR;
    private String name = null;
    private Set<String> filters = null;
    private Map<String, Object> attributes = NO_ATTRIBUTES;

    public CommonInfo(String name, double log10PError, Set<String> filters, Map<String, Object> attributes) {
        this.name = name;
        setLog10PError(log10PError);
        this.filters = filters;
        if ( attributes != null && ! attributes.isEmpty() ) {
            this.attributes = attributes;
        }
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * Sets the name
     *
     * @param name    the name associated with this information
     */
    public void setName(String name) {
        if ( name == null ) throw new IllegalArgumentException("Name cannot be null " + this);
        this.name = name;
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Filter
    //
    // ---------------------------------------------------------------------------------------------------------

    public Set<String> getFiltersMaybeNull() {
        return filters;
    }

    public Set<String> getFilters() {
        return filters == null ? NO_FILTERS : Collections.unmodifiableSet(filters);
    }

    public boolean filtersWereApplied() {
        return filters != null;
    }

    public boolean isFiltered() {
        return filters == null ? false : filters.size() > 0;
    }

    public boolean isNotFiltered() {
        return ! isFiltered();
    }

    public void addFilter(String filter) {
        if ( filters == null ) // immutable -> mutable
            filters = new HashSet<String>();

        if ( filter == null ) throw new IllegalArgumentException("BUG: Attempting to add null filter " + this);
        if ( getFilters().contains(filter) ) throw new IllegalArgumentException("BUG: Attempting to add duplicate filter " + filter + " at " + this);
        filters.add(filter);
    }

    public void addFilters(Collection<String> filters) {
        if ( filters == null ) throw new IllegalArgumentException("BUG: Attempting to add null filters at" + this);
        for ( String f : filters )
            addFilter(f);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with log error rates
    //
    // ---------------------------------------------------------------------------------------------------------

    public boolean hasLog10PError() {
        return getLog10PError() != NO_LOG10_PERROR;
    }

    /**
     * @return the -1 * log10-based error estimate
     */
    public double getLog10PError() { return log10PError; }
    public double getPhredScaledQual() { return getLog10PError() * -10; }

    public void setLog10PError(double log10PError) {
        if ( log10PError > 0 && log10PError != NO_LOG10_PERROR)
            throw new IllegalArgumentException("BUG: log10PError cannot be > 0 : " + this.log10PError);
        if ( Double.isInfinite(this.log10PError) )
            throw new IllegalArgumentException("BUG: log10PError should not be Infinity");
        if ( Double.isNaN(this.log10PError) )
            throw new IllegalArgumentException("BUG: log10PError should not be NaN");
        this.log10PError = log10PError;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with attributes
    //
    // ---------------------------------------------------------------------------------------------------------
    public void clearAttributes() {
        attributes = new HashMap<String, Object>();
    }

    /**
     * @return the attribute map
     */
    public Map<String, Object> getAttributes() {
        return Collections.unmodifiableMap(attributes);
    }

    // todo -- define common attributes as enum

    public void setAttributes(Map<String, ?> map) {
        clearAttributes();
        putAttributes(map);
    }

    public void putAttribute(String key, Object value) {
        putAttribute(key, value, false);
    }

    public void putAttribute(String key, Object value, boolean allowOverwrites) {
        if ( ! allowOverwrites && hasAttribute(key) )
            throw new IllegalStateException("Attempting to overwrite key->value binding: key = " + key + " this = " + this);

        if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
            attributes = new HashMap<String, Object>();

        attributes.put(key, value);
    }

    public void removeAttribute(String key) {
        if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
            attributes = new HashMap<String, Object>();
        attributes.remove(key);
    }

    public void putAttributes(Map<String, ?> map) {
        if ( map != null ) {
            // for efficiency, we can skip the validation if the map is empty
            if ( attributes.size() == 0 ) {
                if ( attributes == NO_ATTRIBUTES ) // immutable -> mutable
                    attributes = new HashMap<String, Object>();
                attributes.putAll(map);
            } else {
                for ( Map.Entry<String, ?> elt : map.entrySet() ) {
                    putAttribute(elt.getKey(), elt.getValue(), false);
                }
            }
        }
    }

    public boolean hasAttribute(String key) {
        return attributes.containsKey(key);
    }

    public int getNumAttributes() {
        return attributes.size();
    }

    /**
     * @param key    the attribute key
     *
     * @return the attribute value for the given key (or null if not set)
     */
    public Object getAttribute(String key) {
        return attributes.get(key);
    }

    public Object getAttribute(String key, Object defaultValue) {
        if ( hasAttribute(key) )
            return attributes.get(key);
        else
            return defaultValue;
    }

    public String getAttributeAsString(String key, String defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof String ) return (String)x;
        return String.valueOf(x); // throws an exception if this isn't a string
    }

    public int getAttributeAsInt(String key, int defaultValue) {
        Object x = getAttribute(key);
        if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
        if ( x instanceof Integer ) return (Integer)x;
        return Integer.valueOf((String)x); // throws an exception if this isn't a string
    }

    public double getAttributeAsDouble(String key, double defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Double ) return (Double)x;
        if ( x instanceof Integer ) return (Integer)x;
        return Double.valueOf((String)x); // throws an exception if this isn't a string
    }

    public boolean getAttributeAsBoolean(String key, boolean defaultValue) {
        Object x = getAttribute(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Boolean ) return (Boolean)x;
        return Boolean.valueOf((String)x); // throws an exception if this isn't a string
    }

//    public String getAttributeAsString(String key)      { return (String.valueOf(getExtendedAttribute(key))); } // **NOTE**: will turn a null Object into the String "null"
//    public int getAttributeAsInt(String key)            { Object x = getExtendedAttribute(key); return x instanceof Integer ? (Integer)x : Integer.valueOf((String)x); }
//    public double getAttributeAsDouble(String key)      { Object x = getExtendedAttribute(key); return x instanceof Double ? (Double)x : Double.valueOf((String)x); }
//    public boolean getAttributeAsBoolean(String key)      { Object x = getExtendedAttribute(key); return x instanceof Boolean ? (Boolean)x : Boolean.valueOf((String)x); }
//    public Integer getAttributeAsIntegerNoException(String key)  { try {return getAttributeAsInt(key);} catch (Exception e) {return null;} }
//    public Double getAttributeAsDoubleNoException(String key)    { try {return getAttributeAsDouble(key);} catch (Exception e) {return null;} }
//    public String getAttributeAsStringNoException(String key)    { if (getExtendedAttribute(key) == null) return null; return getAttributeAsString(key); }
//    public Boolean getAttributeAsBooleanNoException(String key)  { try {return getAttributeAsBoolean(key);} catch (Exception e) {return null;} }
}