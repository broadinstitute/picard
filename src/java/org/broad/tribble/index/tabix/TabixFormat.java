/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package org.broad.tribble.index.tabix;

import org.broad.tribble.TribbleException;

/**
 * The values in a Tabix header that define the format of the file being indexed, e.g. gff, bed, vcf
 */
public class TabixFormat implements Cloneable {
    public static final int ZERO_BASED    = 0x10000;
    public static final int GENERIC_FLAGS = 0;
    public static final int SAM_FLAGS     = 1;
    public static final int VCF_FLAGS     = 2;
    public static final int UCSC_FLAGS    = GENERIC_FLAGS | ZERO_BASED;

    /** Predefined headers for known formats */
    public static TabixFormat GFF = new TabixFormat(GENERIC_FLAGS, 1, 4, 5, '#', 0);
    public static TabixFormat BED = new TabixFormat(UCSC_FLAGS, 1, 2, 3, '#', 0);
    public static TabixFormat PSLTBL = new TabixFormat(UCSC_FLAGS, 15, 17, 18, '#', 0);
    public static TabixFormat SAM = new TabixFormat(SAM_FLAGS, 3, 4, 0, '@', 0);
    public static TabixFormat VCF = new TabixFormat(VCF_FLAGS, 1, 2, 0, '#', 0);

    /** Describes interpretation of file being indexed.  See FLAGS constants above. */
    public int flags;
    /** One-based index of the column in the file being indexed containing the sequence name */
    public int sequenceColumn;
    /** One-based index of the column in the file being indexed containing the start position. */
    public int startPositionColumn;
    /**
     * One-based index of the column in the file being indexed containing the end position. Zero implies
     * there is no end position column.
     */
    public int endPositionColumn;
    /** Lines in the file being indexed that start with this character are ignored. */
    public char metaCharacter;
    /** This is part of the index header, but does not appear to be used. */
    public int numHeaderLinesToSkip;

    public TabixFormat() {
    }

    public TabixFormat(final int flags, final int sequenceColumn, final int startPositionColumn, final int endPositionColumn, final char metaCharacter, final int numHeaderLinesToSkip) {
        this.flags = flags;
        this.sequenceColumn = sequenceColumn;
        this.startPositionColumn = startPositionColumn;
        this.endPositionColumn = endPositionColumn;
        this.metaCharacter = metaCharacter;
        this.numHeaderLinesToSkip = numHeaderLinesToSkip;
    }

    @Override
    public TabixFormat clone() {
        try {
            return (TabixFormat)super.clone();
        } catch (final CloneNotSupportedException e) {
            throw new TribbleException("unpossible!");
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final TabixFormat that = (TabixFormat) o;

        if (endPositionColumn != that.endPositionColumn) return false;
        if (flags != that.flags) return false;
        if (metaCharacter != that.metaCharacter) return false;
        if (numHeaderLinesToSkip != that.numHeaderLinesToSkip) return false;
        if (sequenceColumn != that.sequenceColumn) return false;
        if (startPositionColumn != that.startPositionColumn) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = flags;
        result = 31 * result + sequenceColumn;
        result = 31 * result + startPositionColumn;
        result = 31 * result + endPositionColumn;
        result = 31 * result + (int) metaCharacter;
        result = 31 * result + numHeaderLinesToSkip;
        return result;
    }
}
