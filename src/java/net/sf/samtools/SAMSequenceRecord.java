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
import java.net.URI;
import java.net.URISyntaxException;
import java.util.regex.Pattern;

/**
 * Header information about a reference sequence.  Corresponds to @SQ header record in SAM text header.
 */
public class SAMSequenceRecord extends AbstractSAMHeaderRecord implements Cloneable
{
    private String mSequenceName = null; // Value must be interned() if it's ever set/modified
    private int mSequenceIndex = -1;
    private int mSequenceLength = 0;
    public static final String SEQUENCE_NAME_TAG = "SN";
    public static final String SEQUENCE_LENGTH_TAG = "LN";
    public static final String MD5_TAG = "M5";
    public static final String ASSEMBLY_TAG = "AS";
    public static final String URI_TAG = "UR";
    public static final String SPECIES_TAG = "SP";

    /** If one sequence has this length, and another sequence had a different length, isSameSequence will
     * not complain that they are different sequences. */
    public static final int UNKNOWN_SEQUENCE_LENGTH = 0;


    /**
     * This is not a valid sequence name, because it is reserved in the MRNM field of SAM text format
     * to mean "same reference as RNAME field."
     */
    public static final String RESERVED_MRNM_SEQUENCE_NAME = "=";

    /**
     * The standard tags are stored in text header without type information, because the type of these tags is known.
     */
    public static final Set<String> STANDARD_TAGS =
            new HashSet<String>(Arrays.asList(SEQUENCE_NAME_TAG, SEQUENCE_LENGTH_TAG, ASSEMBLY_TAG, MD5_TAG, URI_TAG,
                                                SPECIES_TAG));

    // Split on any whitespace
    private static Pattern SEQUENCE_NAME_SPLITTER = Pattern.compile("\\s");

    /**
     * @deprecated Use SAMSequenceRecord(final String name, final int sequenceLength) instead.
     * sequenceLength is required for the object to be considered valid.
     */
    public SAMSequenceRecord(final String name) {
        this(name, UNKNOWN_SEQUENCE_LENGTH);
    }

    public SAMSequenceRecord(final String name, final int sequenceLength) {
        if (name != null) {
            if (SEQUENCE_NAME_SPLITTER.matcher(name).find()) {
                throw new SAMException("Sequence name contains invalid character: " + name);
            }
            validateSequenceName(name);
            mSequenceName = name.intern();
        }
        mSequenceLength = sequenceLength;
    }

    public String getSequenceName() { return mSequenceName; }
    // We don't think this method should ever really be used, but we left it here
    // in case we forget and go to implement it later!
    private void setSequenceName(final String name) {
        if (name != null) {
            mSequenceName = name.intern();
        }
        else {
            mSequenceName = null;
        }
    }
    
    public int getSequenceLength() { return mSequenceLength; }
    public void setSequenceLength(final int value) { mSequenceLength = value; }

    public String getAssembly() { return (String) getAttribute("AS"); }
    public void setAssembly(final String value) { setAttribute("AS", value); }

    public String getSpecies() { return (String) getAttribute("SP"); }
    public void setSpecies(final String value) { setAttribute("SP", value); }


    /**
     * @return Index of this record in the sequence dictionary it lives in. 
     */
    public int getSequenceIndex() { return mSequenceIndex; }

    // Private state used only by SAM implementation.
    void setSequenceIndex(final int value) { mSequenceIndex = value; }

    /**
     * Looser comparison than equals().  We look only at sequence index, sequence length, and MD5 tag value
     * (or sequence names, if there is no MD5 tag in either record.
     */
    public boolean isSameSequence(final SAMSequenceRecord that) {
        if (this == that) return true;
        if (that == null) return false;

        if (mSequenceIndex != that.mSequenceIndex) return false;
        // PIC-439.  Allow undefined length.
        if (mSequenceLength != UNKNOWN_SEQUENCE_LENGTH && that.mSequenceLength != UNKNOWN_SEQUENCE_LENGTH && mSequenceLength != that.mSequenceLength)
            return false;
        if (this.getAttribute(SAMSequenceRecord.MD5_TAG) != null && that.getAttribute(SAMSequenceRecord.MD5_TAG) != null) {
            final BigInteger thisMd5 = new BigInteger((String)this.getAttribute(SAMSequenceRecord.MD5_TAG), 16);
            final BigInteger thatMd5 = new BigInteger((String)that.getAttribute(SAMSequenceRecord.MD5_TAG), 16);
            if (!thisMd5.equals(thatMd5)) {
                return false;
            }
        }
        else {
            if (mSequenceName != that.mSequenceName) return false; // Compare using == since we intern() the Strings
        }

        return true;
    }

    private URI makeURI(final String s) throws URISyntaxException {
        URI uri = new URI(s);
        if (uri.getScheme() == null) {
            uri = new URI("file", uri.getUserInfo(), uri.getHost(), uri.getPort(), uri.getPath(), uri.getQuery(), uri.getFragment());
        }
        return uri;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof SAMSequenceRecord)) return false;

        final SAMSequenceRecord that = (SAMSequenceRecord) o;

        if (mSequenceIndex != that.mSequenceIndex) return false;
        if (mSequenceLength != that.mSequenceLength) return false;
        if (!attributesEqual(that)) return false;
        if (mSequenceName != that.mSequenceName) return false; // Compare using == since we intern() the name

        return true;
    }

    @Override
    public int hashCode() {
        return mSequenceName != null ? mSequenceName.hashCode() : 0;
    }

    Set<String> getStandardTags() {
        return STANDARD_TAGS;
    }

    public final SAMSequenceRecord clone() {
        final SAMSequenceRecord ret = new SAMSequenceRecord(this.mSequenceName, this.mSequenceLength);
        ret.mSequenceIndex = this.mSequenceIndex;
        for (final Map.Entry<String, Object> entry : this.getAttributes()) {
            ret.setAttribute(entry.getKey(), entry.getValue());
        }
        return ret;
    }

    /**
     * Truncate sequence name at first space, @ or =.
     */
    public static String truncateSequenceName(final String sequenceName) {
        return SEQUENCE_NAME_SPLITTER.split(sequenceName, 2)[0];
    }

    /**
     * Throw an exception if the sequence name is not valid.
     */
    public static void validateSequenceName(final String name) {
        if (RESERVED_MRNM_SEQUENCE_NAME.equals(name)) {
            throw new SAMException("'" + RESERVED_MRNM_SEQUENCE_NAME + "' is not a valid sequence name");
        }
    }

}

