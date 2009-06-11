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


import java.util.*;
import java.math.BigInteger;

/**
 * Header information about a reference sequence.  Corresponds to @SQ header record in SAM text header.
 */
public class SAMSequenceRecord extends AbstractSAMHeaderRecord implements Cloneable
{
    private String mSequenceName = null;
    private int mSequenceIndex = -1;
    private int mSequenceLength = 0;
    public static final String SEQUENCE_NAME_TAG = "SN";
    public static final String SEQUENCE_LENGTH_TAG = "LN";
    public static final String MD5_TAG = "M5";
    public static final String ASSEMBLY_TAG = "AS";
    public static final String URI_TAG = "UR";
    public static final String SPECIES_TAG = "SP";

    /**
     * The standard tags are stored in text header without type information, because the type of these tags is known.
     */
    public static final Set<String> STANDARD_TAGS =
            new HashSet<String>(Arrays.asList(SEQUENCE_NAME_TAG, SEQUENCE_LENGTH_TAG, ASSEMBLY_TAG, MD5_TAG, URI_TAG,
                                                SPECIES_TAG));

    /**
     * @deprecated Use SAMSequenceRecord(final String name, final int sequenceLength) instead.
     * sequenceLength is required for the object to be considered valid.
     */
    public SAMSequenceRecord(final String name) {
        mSequenceName = name;
    }

    public SAMSequenceRecord(final String name, final int sequenceLength) {
        mSequenceName = name;
        mSequenceLength = sequenceLength;
    }

    public String getSequenceName() {
        return mSequenceName;
    }

    public int getSequenceLength() {
        return mSequenceLength;
    }

    public void setSequenceLength(final int value) {
        mSequenceLength = value;
    }

    public String getAssembly() {
        return (String) getAttribute("AS");
    }

    public void setAssembly(final String value) {
        setAttribute("AS", value);
    }

    public String getSpecies() {
        return (String) getAttribute("SP");
    }

    public void setSpecies(final String value) {
        setAttribute("SP", value);
    }


    /**
     * @return Index of this record in the sequence dictionary it lives in. 
     */
    public int getSequenceIndex() {
        return mSequenceIndex;
    }

    // Private state used only by SAM implementation.
    void setSequenceIndex(final int value) {
        mSequenceIndex = value;
    }

    /**
     * Looser comparison than equals().  If one SAMSequenceRecord has an attribute that the other does not
     * have, that is not considered inequality.  However, if they both have an attribute, but have different
     * values for that atttribute, then they are considered unequal.  This results in an intransitive equality test,
     * i.e. a.isSameSequence(b) && b.isSameSequence(c) does not necessarily imply a.isSameSequence(c)
     */
    public boolean isSameSequence(final SAMSequenceRecord that) {
        if (this == that) return true;
        if (that == null) return false;

        if (mSequenceIndex != that.mSequenceIndex) return false;
        if (mSequenceLength != that.mSequenceLength) return false;
        if (mSequenceName != null ? !mSequenceName.equals(that.mSequenceName) : that.mSequenceName != null)
            return false;
        // If one record has an optional attribute and the other does not, that is not considered inequality.
        
            for (final Map.Entry<String, Object> entry: getAttributes()) {
                final Object thatAttribute = that.getAttribute(entry.getKey());
                if (thatAttribute != null) {
                    if (entry.getKey().equals(SAMSequenceRecord.MD5_TAG)) {
                        // handle the fact that sometimes the MD5 had leading zeroes trimmed, by converting
                        // to BigInteger before comparing.
                        BigInteger thisMd5 = new BigInteger((String)entry.getValue(), 16);
                        BigInteger thatMd5 = new BigInteger((String)thatAttribute, 16);
                        if (!thisMd5.equals(thatMd5)) {
                            return false;
                        }
                    } else if (!entry.getValue().equals(thatAttribute)) {
                        return false;
                    }
                }
        }

        return true;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof SAMSequenceRecord)) return false;

        final SAMSequenceRecord that = (SAMSequenceRecord) o;

        if (mSequenceIndex != that.mSequenceIndex) return false;
        if (mSequenceLength != that.mSequenceLength) return false;
        if (!attributesEqual(that)) return false;
        if (mSequenceName != null ? !mSequenceName.equals(that.mSequenceName) : that.mSequenceName != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = mSequenceName != null ? mSequenceName.hashCode() : 0;
        result = 31 * result + mSequenceIndex;
        result = 31 * result + mSequenceLength;
        result = 31 * result + attributesHashCode();
        return result;
    }

    Set<String> getStandardTags() {
        return STANDARD_TAGS;
    }

    public final SAMSequenceRecord clone() {
        SAMSequenceRecord ret = new SAMSequenceRecord(this.mSequenceName, this.mSequenceLength);
        ret.mSequenceIndex = this.mSequenceIndex;
        for (Map.Entry<String, Object> entry : this.getAttributes()) {
            ret.setAttribute(entry.getKey(), entry.getValue());
        }
        return ret;
    }
}

