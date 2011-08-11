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

/**
 * Holds a SAMRecord attribute and the tagname (in binary form) for that attribute.
 * SAMRecord stores tag name and value in this form, because much String creation is avoided this way.
 * See SAMTagUtil to convert the tag to String form.
 *
 * @author alecw@broadinstitute.org
 */
class SAMBinaryTagAndValue {
    public final short tag;
    public final Object value;
    protected SAMBinaryTagAndValue next = null;

    public SAMBinaryTagAndValue(final short tag, final Object value) {
        this.tag = tag;
        this.value = value;
    }

    @Override public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        return typeSafeEquals((SAMBinaryTagAndValue) o);
    }

    /** Type safe equals method that recurses down the list looking for equality. */
    private boolean typeSafeEquals(final SAMBinaryTagAndValue that) {
        if (this.tag != that.tag) return false;
        if ((this.value == null) ? that.value == null : this.value.equals(that.value)) {
            if (this.next == null) return that.next == null;
            else return this.next.equals(that.next);
        }
        else {
            return false;
        }
    }

    @Override
    public int hashCode() {
        int result = (int) tag;
        result = 31 * result + value.hashCode();
        return result;
    }

    /** Creates and returns a deep copy of the list of tag/values. */
    public SAMBinaryTagAndValue copy() {
        final SAMBinaryTagAndValue retval = new SAMBinaryTagAndValue(this.tag, this.value);
        if (next != null) retval.next = next.copy();
        return retval;
    }

    // The methods below are for implementing a light-weight, single-direction linked list

    public SAMBinaryTagAndValue getNext() { return this.next; }

    /** Inserts at item into the ordered list of attributes and returns the head of the list/sub-list */
    public SAMBinaryTagAndValue insert(final SAMBinaryTagAndValue attr) {
        if (attr == null) return this;
        if (attr.next != null) throw new IllegalStateException("Can only insert single tag/value combinations.");

        if (attr.tag < this.tag) {
            // attr joins the list ahead of this element
            attr.next = this;
            return attr;
        }
        else if (this.tag == attr.tag) {
            // attr replaces this in the list
            attr.next = this.next;
            return attr;
        }
        else if (this.next == null) {
            // attr gets stuck on the end
            this.next = attr;
            return this;
        }
        else {
            // attr gets inserted somewhere in the tail
            this.next = this.next.insert(attr);
            return this;
        }
    }

    /** Removes a tag from the list and returns the new head of the list/sub-list. */
    public SAMBinaryTagAndValue remove(final short tag) {
        if (this.tag == tag) return this.next;
        else {
            if (this.next != null) this.next = this.next.remove(tag);
            return this;
        }
    }

    /** Returns the SAMBinaryTagAndValue that contains the required tag, or null if not contained. */
    public SAMBinaryTagAndValue find(final short tag) {
        if (this.tag == tag) return this;
        else if (this.tag > tag || this.next == null) return null;
        else return this.next.find(tag); 
    }

    public boolean isUnsignedArray() {
        return false;
    }
}
