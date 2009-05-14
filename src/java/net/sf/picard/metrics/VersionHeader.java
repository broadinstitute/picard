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

package net.sf.picard.metrics;

import net.sf.samtools.util.StringUtil;

/**
 * Header that stores information about the version of some piece of software or
 * data used to create the metrics file.  Payload consists of a name or description
 * of the versioned item and a version string.
 *
 * @author Tim Fennell
 */
public class VersionHeader implements Header {
    private String versionedItem;
    private String versionString;

    public void parse(String in) {
        String[] fields = in.split("\t");
        this.versionedItem = fields[0];
        this.versionString = fields[1];
    }

    public String toString() {
        return this.versionedItem + "\t" + this.versionString;
    }

    public String getVersionedItem() { return versionedItem; }
    public void setVersionedItem(String versionedItem) {
        this.versionedItem = StringUtil.assertCharactersNotInString(versionedItem, '\t', '\n');
    }

    public String getVersionString() { return versionString; }
    public void setVersionString(String versionString) {
        this.versionString = StringUtil.assertCharactersNotInString(versionString, '\t', '\n');
    }

    /** Equals method that checks that both the item and version string are equal. */
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VersionHeader that = (VersionHeader) o;

        if (versionString != null ? !versionString.equals(that.versionString) : that.versionString != null)
            return false;
        if (versionedItem != null ? !versionedItem.equals(that.versionedItem) : that.versionedItem != null)
            return false;

        return true;
    }
}
