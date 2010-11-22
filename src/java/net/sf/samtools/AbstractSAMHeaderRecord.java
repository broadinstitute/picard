/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.samtools;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Base class for the various concrete records in a SAM header, providing uniform
 * access to the attributes.
 */
public abstract class AbstractSAMHeaderRecord {
    private final Map<String,String> mAttributes = new HashMap<String, String>();

    public String getAttribute(final String key) {
        return mAttributes.get(key);
    }

    /**
     * Set the given value for the attribute named 'key'.  Replaces an existing value, if any.
     * If value is null, the attribute is removed.
     * Otherwise, the value will be converted to a String with toString.
     * @param key attribute name
     * @param value attribute value
     * @deprecated Use the version that takes a String value instead
     */
    public void setAttribute(final String key, final Object value) {
        setAttribute(key, value == null? null: value.toString());
    }

    /**
     * Set the given value for the attribute named 'key'.  Replaces an existing value, if any.
     * If value is null, the attribute is removed.
     * Supported types are Character, Integer, Float and String.  Byte and Short may also be
     * passed in but they will be converted to Integer.
     * @param key attribute name
     * @param value attribute value
     */
    public void setAttribute(final String key, final String value) {
        if (value == null) {
            mAttributes.remove(key);
        } else {
            mAttributes.put(key, value);
        }
    }
    /**
     * Returns the Set of attributes.
     */
    public Set<Map.Entry<String,String>> getAttributes() {
        return mAttributes.entrySet();
    }


    /**
     * Returns the ID tag (or equivalent) for this header record. The
     * default implementation throws a PicardException to indicate "not implemented".
     */
    public String getId() {
        throw new UnsupportedOperationException("Method not implemented for: " + this.getClass());
    }

    /**
     * For use in the equals() method of the concrete class.
     */
    protected boolean attributesEqual(final AbstractSAMHeaderRecord that) {
        return mAttributes.equals(that.mAttributes);
    }

    /**
     * For use in the hashCode() method of the concrete class.
     */
    protected int attributesHashCode() {
        return (mAttributes != null ? mAttributes.hashCode() : 0);
    }

    /**
     * Standard tags are the tags defined in SAM spec.  These do not have type information in the test
     * representation, because the type information is predefined for each tag.
     * @return list of predefined tags for the concrete SAMHeader record type.
     */
    abstract Set<String> getStandardTags();
}
