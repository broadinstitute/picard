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

import java.util.Map;
import java.util.Set;
import java.util.HashMap;

/**
 * Base class for the various concrete records in a SAM header, providing uniform
 * access to the attributes.
 */
public abstract class AbstractSAMHeaderRecord {
    private final Map<String,Object> mAttributes = new HashMap<String, Object>();

    public Object getAttribute(final String key) {
        return mAttributes.get(key);
    }

    /**
     * Set the given value for the attribute named 'key'.  Replaces an existing value, if any.
     * If value is null, the attribute is removed.
     * @param key attribute name
     * @param value attribute value
     */
    public void setAttribute(final String key, final Object value) {
        if (value == null) {
            mAttributes.remove(key);
        } else {
            mAttributes.put(key, value);
        }
    }

    public Set<Map.Entry<String,Object>> getAttributes() {
        return mAttributes.entrySet();
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
