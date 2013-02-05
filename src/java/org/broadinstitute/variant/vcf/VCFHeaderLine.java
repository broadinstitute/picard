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

import org.broad.tribble.TribbleException;

import java.util.Map;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFHeaderLine
 *         <p/>
 *         A class representing a key=value entry in the VCF header
 */
public class VCFHeaderLine implements Comparable {
    protected static final boolean ALLOW_UNBOUND_DESCRIPTIONS = true;
    protected static final String UNBOUND_DESCRIPTION = "Not provided in original VCF header";

    private String mKey = null;
    private String mValue = null;


    /**
     * create a VCF header line
     *
     * @param key     the key for this header line
     * @param value   the value for this header line
     */
    public VCFHeaderLine(String key, String value) {
        if ( key == null )
            throw new IllegalArgumentException("VCFHeaderLine: key cannot be null");
        mKey = key;
        mValue = value;
    }

    /**
     * Get the key
     *
     * @return the key
     */
    public String getKey() {
        return mKey;
    }

    /**
     * Get the value
     *
     * @return the value
     */
    public String getValue() {
        return mValue;
    }

    public String toString() {
        return toStringEncoding();
    }

    /**
     * Should be overloaded in sub classes to do subclass specific
     *
     * @return the string encoding
     */
    protected String toStringEncoding() {
        return mKey + "=" + mValue;
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFHeaderLine) )
            return false;
        return mKey.equals(((VCFHeaderLine)o).getKey()) && mValue.equals(((VCFHeaderLine)o).getValue());
    }

    public int compareTo(Object other) {
        return toString().compareTo(other.toString());
    }

    /**
     * @param line    the line
     * @return true if the line is a VCF meta data line, or false if it is not
     */
    public static boolean isHeaderLine(String line) {
        return line != null && line.length() > 0 && VCFHeader.HEADER_INDICATOR.equals(line.substring(0,1));
    }

    /**
     * create a string of a mapping pair for the target VCF version
     * @param keyValues a mapping of the key->value pairs to output
     * @return a string, correctly formatted
     */
    public static String toStringEncoding(Map<String, ? extends Object> keyValues) {
        StringBuilder builder = new StringBuilder();
        builder.append("<");
        boolean start = true;
        for (Map.Entry<String,?> entry : keyValues.entrySet()) {
            if (start) start = false;
            else builder.append(",");

            if ( entry.getValue() == null ) throw new TribbleException.InternalCodecException("Header problem: unbound value at " + entry + " from " + keyValues);

            builder.append(entry.getKey());
            builder.append("=");
            builder.append(entry.getValue().toString().contains(",") ||
                           entry.getValue().toString().contains(" ") ||
                           entry.getKey().equals("Description") ? "\""+ entry.getValue() + "\"" : entry.getValue());
        }
        builder.append(">");
        return builder.toString();
    }
}