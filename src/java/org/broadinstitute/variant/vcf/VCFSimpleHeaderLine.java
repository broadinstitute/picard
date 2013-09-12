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

package org.broadinstitute.variant.vcf;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


/**
 * @author ebanks
 * A class representing a key=value entry for simple VCF header types
 */
public class VCFSimpleHeaderLine extends VCFHeaderLine implements VCFIDHeaderLine {

    private String name;
    private Map<String, String> genericFields = new LinkedHashMap<String, String>();

    /**
     * create a VCF filter header line
     *
     * @param key            the key for this header line
     * @param name           the name for this header line
     * @param description    description for this header line
     */
    public VCFSimpleHeaderLine(String key, String name, String description) {
        super(key, "");
        Map<String, String> map = new LinkedHashMap<String, String>(1);
        map.put("Description", description);
        initialize(name, map);
    }

    /**
     * create a VCF info header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     * @param key            the key for this header line
     * @param expectedTagOrdering the tag ordering expected for this header line
     */
    public VCFSimpleHeaderLine(final String line, final VCFHeaderVersion version, final String key, final List<String> expectedTagOrdering) {
        this(key, VCFHeaderLineTranslator.parseLine(version, line, expectedTagOrdering));
    }

    public VCFSimpleHeaderLine(final String key, final Map<String, String> mapping) {
        super(key, "");
        name = mapping.get("ID");
        initialize(name, mapping);
    }

	/**
	 * Returns the String value associated with the given key. Returns null if there is no value. Key
	 * must not be null.
	 */
	String getGenericFieldValue(final String key) {
		return this.genericFields.get(key);
	}

    protected void initialize(String name, Map<String, String> genericFields) {
        if ( name == null || genericFields == null || genericFields.isEmpty() )
            throw new IllegalArgumentException(String.format("Invalid VCFSimpleHeaderLine: key=%s name=%s", super.getKey(), name));
        if ( name.contains("<") || name.contains(">") )
            throw new IllegalArgumentException("VCFHeaderLine: ID cannot contain angle brackets");
        if ( name.contains("=") )
            throw new IllegalArgumentException("VCFHeaderLine: ID cannot contain an equals sign");

        this.name = name;
        this.genericFields.putAll(genericFields);
    }

    protected String toStringEncoding() {
        Map<String, Object> map = new LinkedHashMap<String, Object>();
        map.put("ID", name);
        map.putAll(genericFields);
        return getKey() + "=" + VCFHeaderLine.toStringEncoding(map);
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFSimpleHeaderLine) )
            return false;
        VCFSimpleHeaderLine other = (VCFSimpleHeaderLine)o;
        if ( !name.equals(other.name) || genericFields.size() != other.genericFields.size() )
            return false;
        for ( Map.Entry<String, String> entry : genericFields.entrySet() ) {
            if ( !entry.getValue().equals(other.genericFields.get(entry.getKey())) )
                return false;
        }
        
        return true;       
    }

    public String getID() {
        return name;
    }
}