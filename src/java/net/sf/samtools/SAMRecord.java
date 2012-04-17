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


import net.sf.samtools.util.CoordMath;
import net.sf.samtools.util.StringUtil;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


/**
 * Java binding for a SAM file record.  c.f. http://samtools.sourceforge.net/SAM1.pdf
 *
 * The presence of reference name/reference index and alignment start
 * do not necessarily mean that a read is aligned.  Those values may merely be set to force a SAMRecord
 * to appear in a certain place in the sort order.  The readUnmappedFlag must be checked to determine whether
 * or not a read is mapped.  Only if the readUnmappedFlag is false can the reference name/index and alignment start
 * be interpreted as indicating an actual alignment position.
 *
 * Likewise, presence of mate reference name/index and mate alignment start do not necessarily mean that the
 * mate is aligned.  These may be set for an unaligned mate if the mate has been forced into a particular place
 * in the sort order per the above paragraph.  Only if the mateUnmappedFlag is false can the mate reference name/index
 * and mate alignment start be interpreted as indicating the actual alignment position of the mate.
 *
 * Note also that there are a number of getters & setters that are linked, i.e. they present different representations
 * of the same underlying data.  In these cases there is typically a representation that is preferred because it
 * ought to be faster than some other representation.  The following are the preferred representations:
 *
 * getReadNameLength() is preferred to getReadName().length()
 * get/setReadBases() is preferred to get/setReadString()
 * get/setBaseQualities() is preferred to get/setBaseQualityString()
 * get/setReferenceIndex() is preferred to get/setReferenceName()
 * get/setMateReferenceIndex() is preferred to get/setMateReferenceName()
 * getCigarLength() is preferred to getCigar().getNumElements()
 * get/setCigar() is preferred to get/setCigarString()
 *
 * Note that setIndexingBin() need not be called when writing SAMRecords.  It will be computed as necessary.  It is only
 * present as an optimization in the event that the value is already known and need not be computed.
 *
 * setHeader() need not be called when writing SAMRecords.  It may be convenient to call it, however, because
 * get/setReferenceIndex() and get/setMateReferenceIndex() must have access to the SAM header, either as an argument
 * or previously passed to setHeader().
 *
 * setHeader() is called by the SAM reading code, so the get/setReferenceIndex() and get/setMateReferenceIndex()
 * methods will have access to the sequence dictionary.
 *
 * Some of the get() methods return values that are mutable, due to the limitations of Java.  A caller should
 * never change the value returned by a get() method.  If you want to change the value of some attribute of a
 * SAMRecord, create a new value object and call the appropriate set() method.
 *
 * By default, extensive validation of SAMRecords is done when they are read.  Very limited validation is done when
 * values are set onto SAMRecords.
 */
public class SAMRecord implements Cloneable
{
    /**
     * Alignment score for a good alignment, but where computing a Phred-score is not feasible. 
     */
    public static final int UNKNOWN_MAPPING_QUALITY = 255;

    /**
     * Alignment score for an unaligned read.
     */
    public static final int NO_MAPPING_QUALITY = 0;

    /**
     * If a read has this reference name, it is unaligned, but not all unaligned reads have
     * this reference name (see above).
     */
    public static final String NO_ALIGNMENT_REFERENCE_NAME = "*";

    /**
     * If a read has this reference index, it is unaligned, but not all unaligned reads have
     * this reference index (see above).
     */
    public static final int NO_ALIGNMENT_REFERENCE_INDEX = -1;

    /**
     * Cigar string for an unaligned read.
     */
    public static final String NO_ALIGNMENT_CIGAR = "*";

    /**
     * If a read has reference name "*", it will have this value for position.
     */
    public static final int NO_ALIGNMENT_START = 0;

    /**
     * This should rarely be used, since a read with no sequence doesn't make much sense.
     */
    public static final byte[] NULL_SEQUENCE = new byte[0];

    public static final String NULL_SEQUENCE_STRING = "*";

    /**
     * This should rarely be used, since all reads should have quality scores.
     */
    public static final byte[] NULL_QUALS = new byte[0];
    public static final String NULL_QUALS_STRING = "*";

    /**
     * abs(insertSize) must be <= this
     */
    public static final int MAX_INSERT_SIZE = 1<<29;

    /**
     * It is not necessary in general to use the flag constants, because there are getters
     * & setters that handles these symbolically.
     */
    private static final int READ_PAIRED_FLAG = 0x1;
    private static final int PROPER_PAIR_FLAG = 0x2;
    private static final int READ_UNMAPPED_FLAG = 0x4;
    private static final int MATE_UNMAPPED_FLAG = 0x8;
    private static final int READ_STRAND_FLAG = 0x10;
    private static final int MATE_STRAND_FLAG = 0x20;
    private static final int FIRST_OF_PAIR_FLAG = 0x40;
    private static final int SECOND_OF_PAIR_FLAG = 0x80;
    private static final int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    private static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    private static final int DUPLICATE_READ_FLAG = 0x400;


    private String mReadName = null;
    private byte[] mReadBases = NULL_SEQUENCE;
    private byte[] mBaseQualities = NULL_QUALS;
    private String mReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
    private int mAlignmentStart = NO_ALIGNMENT_START;
    private transient int mAlignmentEnd = NO_ALIGNMENT_START;
    private int mMappingQuality = NO_MAPPING_QUALITY;
    private String mCigarString = NO_ALIGNMENT_CIGAR;
    private Cigar mCigar = null;
    private List<AlignmentBlock> mAlignmentBlocks = null;
    private int mFlags = 0;
    private String mMateReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
    private int mMateAlignmentStart = 0;
    private int mInferredInsertSize = 0;
    private SAMBinaryTagAndValue mAttributes = null;
    protected Integer mReferenceIndex = null;
    protected Integer mMateReferenceIndex = null;
    private Integer mIndexingBin = null;

    /**
     * Some attributes (e.g. CIGAR) are not decoded immediately.  Use this to decide how to validate when decoded.
     */
    private SAMFileReader.ValidationStringency mValidationStringency = SAMFileReader.ValidationStringency.SILENT;

    private SAMFileSource mFileSource;
    private SAMFileHeader mHeader = null;

    public SAMRecord(final SAMFileHeader header) {
        mHeader = header;
    }

    public String getReadName() {
        return mReadName;
    }

    /**
     * This method is preferred over getReadName().length(), because for BAMRecord
     * it may be faster.
     * @return length not including a null terminator.
     */
    public int getReadNameLength() {
        return mReadName.length();
    }

    public void setReadName(final String value) {
        mReadName = value;
    }

    /**
     * @return read sequence as a string of ACGTN=.
     */
    public String getReadString() {
        final byte[] readBases = getReadBases();
        if (readBases.length == 0) {
            return NULL_SEQUENCE_STRING;
        }
        return StringUtil.bytesToString(readBases);
    }

    public void setReadString(final String value) {
        if (NULL_SEQUENCE_STRING.equals(value)) {
            mReadBases = NULL_SEQUENCE;
        } else {
            final byte[] bases = StringUtil.stringToBytes(value);
            SAMUtils.normalizeBases(bases);
            setReadBases(bases);
        }
    }


    /**
     * Do not modify the value returned by this method.  If you want to change the bases, create a new
     * byte[] and call setReadBases() or call setReadString().
     * @return read sequence as ASCII bytes ACGTN=.
     */
    public byte[] getReadBases() {
        return mReadBases;
    }

    public void setReadBases(final byte[] value) {
        mReadBases = value;
    }

    /**
     * This method is preferred over getReadBases().length, because for BAMRecord it may be faster.
     * @return number of bases in the read.
     */
    public int getReadLength() {
        return getReadBases().length;
    }

    /**
     * @return Base qualities, encoded as a FASTQ string.
     */
    public String getBaseQualityString() {
        if (Arrays.equals(NULL_QUALS, getBaseQualities())) {
            return NULL_QUALS_STRING;
        }
        return SAMUtils.phredToFastq(getBaseQualities());
    }

    public void setBaseQualityString(final String value) {
        if (NULL_QUALS_STRING.equals(value)) {
            setBaseQualities(NULL_QUALS);
        } else {
            setBaseQualities(SAMUtils.fastqToPhred(value));
        }
    }

    /**
     * Do not modify the value returned by this method.  If you want to change the qualities, create a new
     * byte[] and call setBaseQualities() or call setBaseQualityString().
     * @return Base qualities, as binary phred scores (not ASCII).
     */
    public byte[] getBaseQualities() {
        return mBaseQualities;
    }

    public void setBaseQualities(final byte[] value) {
        mBaseQualities = value;
    }

    /**
     * If the original base quality scores have been store in the "OQ" tag will return the numeric
     * score as a byte[]
     */
    public byte[] getOriginalBaseQualities() {
        final String oqString = (String) getAttribute("OQ");
        if (oqString != null && oqString.length() > 0) {
            return SAMUtils.fastqToPhred(oqString);
        }
        else {
            return null;
        }
    }

    /**
     * Sets the original base quality scores into the "OQ" tag as a String.  Supplied value should be
     * as phred-scaled numeric qualities.
     */
    public void setOriginalBaseQualities(final byte[] oq) {
        setAttribute("OQ", SAMUtils.phredToFastq(oq));
    }

    private static boolean hasReferenceName(final Integer referenceIndex, final String referenceName) {
        return (referenceIndex != null && referenceIndex != NO_ALIGNMENT_REFERENCE_INDEX) ||
                !NO_ALIGNMENT_REFERENCE_NAME.equals(referenceName);
    }

    /**
     * @return true if this SAMRecord has a reference, either as a String or index (or both).
     */
    private boolean hasReferenceName() {
        return hasReferenceName(mReferenceIndex, mReferenceName);
    }

    /**
     * @return true if this SAMRecord has a mate reference, either as a String or index (or both).
     */
    private boolean hasMateReferenceName() {
        return hasReferenceName(mMateReferenceIndex, mMateReferenceName);
    }

    /**
     * @return Reference name, or null if record has no reference.
     */
    public String getReferenceName() {
        return mReferenceName;
    }

    public void setReferenceName(final String value) {
        /* String.intern() is surprisingly expensive, so avoid it by looking up in sequence dictionary if possible */
        if (NO_ALIGNMENT_REFERENCE_NAME.equals(value)) {
            mReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
            mReferenceIndex = NO_ALIGNMENT_REFERENCE_INDEX;
            return;
        } else if (mHeader != null) {
            int referenceIndex = mHeader.getSequenceIndex(value);
            if (referenceIndex != -1) {
                setReferenceIndex(referenceIndex);
                return;
            }
        }
        // Drop through from above if nothing done.
        mReferenceName = value.intern();
        mReferenceIndex = null;
    }

    /**
     * @return index of the reference sequence for this read in the sequence dictionary, or -1
     * if read has no reference sequence set, or if a String reference name is not found in the sequence index..
     */
    public Integer getReferenceIndex() {
        if (mReferenceIndex == null) {
            if (mReferenceName == null) {
                mReferenceIndex = NO_ALIGNMENT_REFERENCE_INDEX;
            } else if (NO_ALIGNMENT_REFERENCE_NAME.equals(mReferenceName)) {
                mReferenceIndex = NO_ALIGNMENT_REFERENCE_INDEX;
            } else {
                mReferenceIndex = mHeader.getSequenceIndex(mReferenceName);
            }
        }
        return mReferenceIndex;
    }

    /**
     * @param referenceIndex Must either equal -1 (indicating no reference), or exist in the sequence dictionary
     * in the header associated with this record.
     */
    public void setReferenceIndex(final int referenceIndex) {
        mReferenceIndex = referenceIndex;
        if (mReferenceIndex == NO_ALIGNMENT_REFERENCE_INDEX) {
            mReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
        } else {
            try {
                mReferenceName = mHeader.getSequence(referenceIndex).getSequenceName();
            } catch (NullPointerException e) {
                throw new IllegalArgumentException("Reference index " + referenceIndex + " not found in sequence dictionary.", e);
            }
        }
    }

    /**
     * @return Mate reference name, or null if one is not assigned.
     */
    public String getMateReferenceName() {
        return mMateReferenceName;
    }

    public void setMateReferenceName(final String mateReferenceName) {
        /* String.intern() is surprisingly expensive, so avoid it by looking up in sequence dictionary if possible */
        if (NO_ALIGNMENT_REFERENCE_NAME.equals(mateReferenceName)) {
            mMateReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
            mMateReferenceIndex = NO_ALIGNMENT_REFERENCE_INDEX;
            return;
        } else if (mHeader != null) {
            int referenceIndex = mHeader.getSequenceIndex(mateReferenceName);
            if (referenceIndex != -1) {
                setMateReferenceIndex(referenceIndex);
                return;
            }
        }
        // Drop through from above if nothing done.
        this.mMateReferenceName = mateReferenceName.intern();
        mMateReferenceIndex = null;
    }

    /**
     * @return index of the reference sequence for this read's mate in the sequence dictionary, or -1
     * if mate has no reference sequence set.
     */
    public Integer getMateReferenceIndex() {
        if (mMateReferenceIndex == null) {
            if (mMateReferenceName == null) {
                mMateReferenceIndex = NO_ALIGNMENT_REFERENCE_INDEX;
            } else if (NO_ALIGNMENT_REFERENCE_NAME.equals(mMateReferenceName)){
                mMateReferenceIndex = NO_ALIGNMENT_REFERENCE_INDEX;
            } else {
                mMateReferenceIndex = mHeader.getSequenceIndex(mMateReferenceName);
            }
        }
        return mMateReferenceIndex;
    }

    /**
     * @param referenceIndex Must either equal -1 (indicating no reference), or exist in the sequence dictionary
     * in the header associated with this record.
     */
    public void setMateReferenceIndex(final int referenceIndex) {
        mMateReferenceIndex = referenceIndex;
        if (mMateReferenceIndex == NO_ALIGNMENT_REFERENCE_INDEX) {
            mMateReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
        } else {
            try {
                mMateReferenceName = mHeader.getSequence(referenceIndex).getSequenceName();
            } catch (NullPointerException e) {
                throw new IllegalArgumentException("Reference index " + referenceIndex + " not found in sequence dictionary.", e);
            }
        }
    }

    /**
     * @return 1-based inclusive leftmost position of the clippped sequence, or 0 if there is no position.
     */
    public int getAlignmentStart() {
        return mAlignmentStart;
    }

    /**
     * @param value 1-based inclusive leftmost position of the clippped sequence, or 0 if there is no position.
     */
    public void setAlignmentStart(final int value) {
        mAlignmentStart = value;
        // Clear cached alignment end
        mAlignmentEnd = NO_ALIGNMENT_START;
        // Change to alignmentStart could change indexing bin
        setIndexingBin(null);
    }

    /**
     * @return 1-based inclusive rightmost position of the clippped sequence, or 0 read if unmapped.
     */
    public int getAlignmentEnd() {
        if (getReadUnmappedFlag()) {
            return NO_ALIGNMENT_START;
        }
        else if (this.mAlignmentEnd == NO_ALIGNMENT_START) {
            this.mAlignmentEnd = mAlignmentStart + getCigar().getReferenceLength() - 1;
        }

        return this.mAlignmentEnd;
    }

    /**
     * @return the alignment start (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment start of 100 but the first 4 bases were clipped (hard or soft clipped)
     * then this method will return 96.
     *
     * Invalid to call on an unmapped read.
     */
    public int getUnclippedStart() {
        int pos = getAlignmentStart();

        for (final CigarElement cig : getCigar().getCigarElements()) {
            final CigarOperator op = cig.getOperator();
            if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
                pos -= cig.getLength();
            }
            else {
                break;
            }
        }

        return pos;
    }

    /**
     * @return the alignment end (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment end of 100 but the last 7 bases were clipped (hard or soft clipped)
     * then this method will return 107.
     *
     * Invalid to call on an unmapped read.
     */
    public int getUnclippedEnd() {
        int pos = getAlignmentEnd();
        final List<CigarElement> cigs = getCigar().getCigarElements();
        for (int i=cigs.size() - 1; i>=0; --i) {
            final CigarElement cig = cigs.get(i);
            final CigarOperator op = cig.getOperator();

            if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
                pos += cig.getLength();
            }
            else {
                break;
            }
        }

        return pos;               
    }

    /**
     * @return 1-based inclusive reference position of the unclippped sequence at a given offset,
     *         or 0 if there is no position.
     *         For example, given the sequence NNNAAACCCGGG, cigar 3S9M, and an alignment start of 1,
     *         and a (1-based)offset 10 (start of GGG) it returns 7 (1-based offset starting after the soft clip.
     *         For example: given the sequence AAACCCGGGTTT, cigar 4M1D6M, an alignment start of 1,
     *         an offset of 4 returns reference position 4, an offset of 5 returns reference position 6.
     *         Another example: given the sequence AAACCCGGGTTT, cigar 4M1I6M, an alignment start of 1,
     *         an offset of 4 returns reference position 4, an offset of 5 returns 0.
     * @offset 1-based location within the unclipped sequence
     */
    public int getReferencePositionAtReadPosition(final int offset) {

        if (offset == 0) return 0;

        for (final AlignmentBlock alignmentBlock : getAlignmentBlocks()) {
            if (CoordMath.getEnd(alignmentBlock.getReadStart(), alignmentBlock.getLength()) < offset) {
                continue;
            } else if (offset < alignmentBlock.getReadStart()) {
                return 0;
            } else {
                return alignmentBlock.getReferenceStart() + offset - alignmentBlock.getReadStart();
            }
        }
        return 0; // offset not located in an alignment block
    }

    /**
     * Unsupported.  This property is derived from alignment start and CIGAR. 
     */
    public void setAlignmentEnd(final int value) {
        throw new UnsupportedOperationException("Not supported: setAlignmentEnd");
    }

    /**
     * @return 1-based inclusive leftmost position of the clippped mate sequence, or 0 if there is no position.
     */
    public int getMateAlignmentStart() {
        return mMateAlignmentStart;
    }

    public void setMateAlignmentStart(final int mateAlignmentStart) {
        this.mMateAlignmentStart = mateAlignmentStart;
    }

    /**
     * @return insert size (difference btw 5' end of read & 5' end of mate), if possible, else 0.
     * Negative if mate maps to lower position than read.
     */
    public int getInferredInsertSize() {
        return mInferredInsertSize;
    }

    public void setInferredInsertSize(final int inferredInsertSize) {
        this.mInferredInsertSize = inferredInsertSize;
    }

    /**
     * @return phred scaled mapping quality.  255 implies valid mapping but quality is hard to compute.
     */
    public int getMappingQuality() {
        return mMappingQuality;
    }

    public void setMappingQuality(final int value) {
        mMappingQuality = value;
    }

    public String getCigarString() {
        if (mCigarString == null && getCigar() != null) {
            mCigarString = TextCigarCodec.getSingleton().encode(getCigar());
        }
        return mCigarString;
    }

    public void setCigarString(final String value) {
        mCigarString = value;
        mCigar = null;
        mAlignmentBlocks = null;
        // Clear cached alignment end
        mAlignmentEnd = NO_ALIGNMENT_START;
        // Change to cigar could change alignmentEnd, and thus indexing bin
        setIndexingBin(null);
    }

    /**
     * Do not modify the value returned by this method.  If you want to change the Cigar, create a new
     * Cigar and call setCigar() or call setCigarString()
     * @return Cigar object for the read, or null if there is none.
     */
    public Cigar getCigar() {
        if (mCigar == null && mCigarString != null) {
            mCigar = TextCigarCodec.getSingleton().decode(mCigarString);
            if (getValidationStringency() != SAMFileReader.ValidationStringency.SILENT && !this.getReadUnmappedFlag()) {
                // Don't know line number, and don't want to force read name to be decoded.
                SAMUtils.processValidationErrors(validateCigar(-1L), -1L, getValidationStringency());
            }
        }
        return mCigar;
    }

    /**
     * This method is preferred over getCigar().getNumElements(), because for BAMRecord it may be faster.
     * @return number of cigar elements (number + operator) in the cigar string.
     */
    public int getCigarLength() {
        return getCigar().numCigarElements();
    }

    public void setCigar(final Cigar cigar) {
        initializeCigar(cigar);
        // Change to cigar could change alignmentEnd, and thus indexing bin
        setIndexingBin(null);
    }

    /**
     * For setting the Cigar string when BAMRecord has decoded it.  Use this rather than setCigar()
     * so that indexing bin doesn't get clobbered.
     */
    protected void initializeCigar(final Cigar cigar) {
        this.mCigar = cigar;
        mCigarString = null;
        mAlignmentBlocks = null;
        // Clear cached alignment end
        mAlignmentEnd = NO_ALIGNMENT_START;
    }

    /**
     * Get the SAMReadGroupRecord for this SAMRecord.
     * @return The SAMReadGroupRecord from the SAMFileHeader for this SAMRecord, or null if
     * 1) this record has no RG tag, or 2) the header doesn't contain the read group with
     * the given ID.
     * @throws NullPointerException if this.getHeader() returns null.
     * @throws ClassCastException if RG tag does not have a String value.
     */
    public SAMReadGroupRecord getReadGroup() {
        final String rgId = (String)getAttribute(SAMTagUtil.getSingleton().RG);
        if (rgId == null) {
            return null;
        }
        return getHeader().getReadGroup(rgId);
    }

    /**
     * It is preferrable to use the get*Flag() methods that handle the flag word symbolically.
     */
    public int getFlags() {
        return mFlags;
    }

    public void setFlags(final int value) {
        mFlags = value;
        // Could imply change to readUnmapped flag, which could change indexing bin
        setIndexingBin(null);
    }

    /**
     * the read is paired in sequencing, no matter whether it is mapped in a pair.
     */
    public boolean getReadPairedFlag() {
        return (mFlags & READ_PAIRED_FLAG) != 0;
    }

    private void requireReadPaired() {
        if (!getReadPairedFlag()) {
            throw new IllegalStateException("Inappropriate call if not paired read");
        }
    }

    /**
     * the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment).
     */
    public boolean getProperPairFlag() {
        requireReadPaired();
        return getProperPairFlagUnchecked();
    }

    private boolean getProperPairFlagUnchecked() {
        return (mFlags & PROPER_PAIR_FLAG) != 0;
    }

    /**
     * the query sequence itself is unmapped.
     */
    public boolean getReadUnmappedFlag() {
        return (mFlags & READ_UNMAPPED_FLAG) != 0;
    }

    /**
     * the mate is unmapped.
     */
    public boolean getMateUnmappedFlag() {
        requireReadPaired();
        return getMateUnmappedFlagUnchecked();
    }

    private boolean getMateUnmappedFlagUnchecked() {
        return (mFlags & MATE_UNMAPPED_FLAG) != 0;
    }

    /**
     * strand of the query (false for forward; true for reverse strand).
     */
    public boolean getReadNegativeStrandFlag() {
        return (mFlags & READ_STRAND_FLAG) != 0;
    }

    /**
     * strand of the mate (false for forward; true for reverse strand).
     */
    public boolean getMateNegativeStrandFlag() {
        requireReadPaired();
        return getMateNegativeStrandFlagUnchecked();
    }

    private boolean getMateNegativeStrandFlagUnchecked() {
        return (mFlags & MATE_STRAND_FLAG) != 0;
    }

    /**
     * the read is the first read in a pair.
     */
    public boolean getFirstOfPairFlag() {
        requireReadPaired();
        return getFirstOfPairFlagUnchecked();
    }

    private boolean getFirstOfPairFlagUnchecked() {
        return (mFlags & FIRST_OF_PAIR_FLAG) != 0;
    }

    /**
     * the read is the second read in a pair.
     */
    public boolean getSecondOfPairFlag() {
        requireReadPaired();
        return getSecondOfPairFlagUnchecked();
    }

    private boolean getSecondOfPairFlagUnchecked() {
        return (mFlags & SECOND_OF_PAIR_FLAG) != 0;
    }

    /**
     * the alignment is not primary (a read having split hits may have multiple primary alignment records).
     */
    public boolean getNotPrimaryAlignmentFlag() {
        return (mFlags & NOT_PRIMARY_ALIGNMENT_FLAG) != 0;
    }

    /**
     * the read fails platform/vendor quality checks.
     */
    public boolean getReadFailsVendorQualityCheckFlag() {
        return (mFlags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }

    /**
     * the read is either a PCR duplicate or an optical duplicate.
     */
    public boolean getDuplicateReadFlag() {
        return (mFlags & DUPLICATE_READ_FLAG) != 0;
    }

    /**
     * the read is paired in sequencing, no matter whether it is mapped in a pair.
     */
    public void setReadPairedFlag(final boolean flag) {
        setFlag(flag, READ_PAIRED_FLAG);
    }

    /**
     * the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment).
     */
    public void setProperPairFlag(final boolean flag) {
        setFlag(flag, PROPER_PAIR_FLAG);
    }

    /**
     * the query sequence itself is unmapped.  This method name is misspelled.
     * Use setReadUnmappedFlag instead.
     * @deprecated
     */
    public void setReadUmappedFlag(final boolean flag) {
        setReadUnmappedFlag(flag);
    }

    /**
     * the query sequence itself is unmapped.
     */
    public void setReadUnmappedFlag(final boolean flag) {
        setFlag(flag, READ_UNMAPPED_FLAG);
        // Change to readUnmapped could change indexing bin
        setIndexingBin(null);
    }

    /**
     * the mate is unmapped.
     */
    public void setMateUnmappedFlag(final boolean flag) {
        setFlag(flag, MATE_UNMAPPED_FLAG);
    }

    /**
     * strand of the query (false for forward; true for reverse strand).
     */
    public void setReadNegativeStrandFlag(final boolean flag) {
        setFlag(flag, READ_STRAND_FLAG);
    }

    /**
     * strand of the mate (false for forward; true for reverse strand).
     */
    public void setMateNegativeStrandFlag(final boolean flag) {
        setFlag(flag, MATE_STRAND_FLAG);
    }

    /**
     * the read is the first read in a pair.
     */
    public void setFirstOfPairFlag(final boolean flag) {
        setFlag(flag, FIRST_OF_PAIR_FLAG);
    }

    /**
     * the read is the second read in a pair.
     */
    public void setSecondOfPairFlag(final boolean flag) {
        setFlag(flag, SECOND_OF_PAIR_FLAG);
    }

    /**
     * the alignment is not primary (a read having split hits may have multiple primary alignment records).
     */
    public void setNotPrimaryAlignmentFlag(final boolean flag) {
        setFlag(flag, NOT_PRIMARY_ALIGNMENT_FLAG);
    }

    /**
     * the read fails platform/vendor quality checks.
     */
    public void setReadFailsVendorQualityCheckFlag(final boolean flag) {
        setFlag(flag, READ_FAILS_VENDOR_QUALITY_CHECK_FLAG);
    }

    /**
     * the read is either a PCR duplicate or an optical duplicate.
     */
    public void setDuplicateReadFlag(final boolean flag) {
        setFlag(flag, DUPLICATE_READ_FLAG);
    }

    private void setFlag(final boolean flag, final int bit) {
        if (flag) {
            mFlags |= bit;
        } else {
            mFlags &= ~bit;
        }
    }

    public SAMFileReader.ValidationStringency getValidationStringency() {
        return mValidationStringency;
    }

    /**
     * Control validation of lazily-decoded elements.
     */
    public void setValidationStringency(final SAMFileReader.ValidationStringency validationStringency) {
        this.mValidationStringency = validationStringency;
    }

    /**
     * Get the value for a SAM tag.
     * WARNING: Some value types (e.g. byte[]) are mutable.  It is dangerous to change one of these values in
     * place, because some SAMRecord implementations keep track of when attributes have been changed.  If you
     * want to change an attribute value, call setAttribute() to replace the value.
     *
     * @param tag Two-character tag name.
     * @return Appropriately typed tag value, or null if the requested tag is not present.
     */
    public Object getAttribute(final String tag) {
        return getAttribute(SAMTagUtil.getSingleton().makeBinaryTag(tag));
    }

    /**
     * Get the tag value and attempt to coerce it into the requested type.
     * @param tag The requested tag.
     * @return The value of a tag, converted into an Integer if possible.
     * @throws RuntimeException If the value is not an integer type, or will not fit in an Integer.
     */
    public Integer getIntegerAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof Integer) {
            return (Integer)val;
        }
        if (!(val instanceof Number)) {
            throw new RuntimeException("Value for tag " + tag + " is not Number: " + val.getClass());
        }
        final long longVal = ((Number)val).longValue();
        if (longVal < Integer.MIN_VALUE || longVal > Integer.MAX_VALUE) {
            throw new RuntimeException("Value for tag " + tag + " is not in Integer range: " + longVal);
        }
        return (int)longVal;
    }

    /**
     * Get the tag value and attempt to coerce it into the requested type.
     * @param tag The requested tag.
     * @return The value of a tag, converted into a Short if possible.
     * @throws RuntimeException If the value is not an integer type, or will not fit in a Short.
     */
    public Short getShortAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof Short) {
            return (Short)val;
        }
        if (!(val instanceof Number)) {
            throw new RuntimeException("Value for tag " + tag + " is not Number: " + val.getClass());
        }
        final long longVal = ((Number)val).longValue();
        if (longVal < Short.MIN_VALUE || longVal > Short.MAX_VALUE) {
            throw new RuntimeException("Value for tag " + tag + " is not in Short range: " + longVal);
        }
        return (short)longVal;
    }

    /**
     * Get the tag value and attempt to coerce it into the requested type.
     * @param tag The requested tag.
     * @return The value of a tag, converted into a Byte if possible.
     * @throws RuntimeException If the value is not an integer type, or will not fit in a Byte.
     */
    public Byte getByteAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof Byte) {
            return (Byte)val;
        }
        if (!(val instanceof Number)) {
            throw new RuntimeException("Value for tag " + tag + " is not Number: " + val.getClass());
        }
        final long longVal = ((Number)val).longValue();
        if (longVal < Byte.MIN_VALUE || longVal > Byte.MAX_VALUE) {
            throw new RuntimeException("Value for tag " + tag + " is not in Short range: " + longVal);
        }
        return (byte)longVal;
    }

    public String getStringAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof String) {
            return (String)val;
        }
        throw new SAMException("Value for tag " + tag + " is not a String: " + val.getClass());
    }

    public Character getCharacterAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof Character) {
            return (Character)val;
        }
        throw new SAMException("Value for tag " + tag + " is not a Character: " + val.getClass());
    }

    public Float getFloatAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof Float) {
            return (Float)val;
        }
        throw new SAMException("Value for tag " + tag + " is not a Float: " + val.getClass());
    }

    /** Will work for signed byte array, unsigned byte array, or old-style hex array */
    public byte[] getByteArrayAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof byte[]) {
            return (byte[])val;
        }
        throw new SAMException("Value for tag " + tag + " is not a byte[]: " + val.getClass());
    }

    public byte[] getUnsignedByteArrayAttribute(final String tag) {
        final byte[] ret = getByteArrayAttribute(tag);
        if (ret != null) requireUnsigned(tag);
        return ret;
    }

    /** Will work for signed byte array or old-style hex array */
    public byte[] getSignedByteArrayAttribute(final String tag) {
        final byte[] ret = getByteArrayAttribute(tag);
        if (ret != null) requireSigned(tag);
        return ret;
    }

    public short[] getUnsignedShortArrayAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof short[]) {
            requireUnsigned(tag);
            return (short[]) val;
        }
        throw new SAMException("Value for tag " + tag + " is not a short[]: " + val.getClass());
    }

    public short[] getSignedShortArrayAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof short[]) {
            requireSigned(tag);
            return (short[]) val;
        }
        throw new SAMException("Value for tag " + tag + " is not a short[]: " + val.getClass());
    }

    public int[] getUnsignedIntArrayAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof int[]) {
            requireUnsigned(tag);
            return (int[]) val;
        }
        throw new SAMException("Value for tag " + tag + " is not a int[]: " + val.getClass());
    }

    public int[] getSignedIntArrayAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val == null) return null;
        if (val instanceof int[]) {
            requireSigned(tag);
            return (int[]) val;
        }
        throw new SAMException("Value for tag " + tag + " is not a int[]: " + val.getClass());
    }

    public float[] getFloatArrayAttribute(final String tag) {
        final Object val = getAttribute(tag);
        if (val != null && !(val instanceof float[])) {
            throw new SAMException("Value for tag " + tag + " is not a float[]: " + val.getClass());
        }
        return (float[]) val;
    }

    /**
     * @return True if this tag is an unsigned array, else false.
     * @throws SAMException if the tag is not present.
     */
    public boolean isUnsignedArrayAttribute(final String tag) {
        final SAMBinaryTagAndValue tmp = this.mAttributes.find(SAMTagUtil.getSingleton().makeBinaryTag(tag));
        if (tmp != null) return tmp.isUnsignedArray();
        throw new SAMException("Tag " + tag + " is not present in this SAMRecord");
    }

    private void requireSigned(final String tag) {
        if (isUnsignedArrayAttribute(tag))  throw new SAMException("Value for tag " + tag + " is not signed");
    }

    private void requireUnsigned(final String tag) {
        if (!isUnsignedArrayAttribute(tag))  throw new SAMException("Value for tag " + tag + " is not unsigned");
    }

    /**
     * @see SAMRecord#getAttribute(java.lang.String)
     * @param tag Binary representation of a 2-char String tag as created by SAMTagUtil.
     */
    public Object getAttribute(final short tag) {
        if (this.mAttributes == null) return null;
        else {
            final SAMBinaryTagAndValue tmp = this.mAttributes.find(tag);
            if (tmp != null) return tmp.value;
            else return null;
        }
    }

    /**
     * Set a named attribute onto the SAMRecord.  Passing a null value causes the attribute to be cleared.
     * @param tag two-character tag name.  See http://samtools.sourceforge.net/SAM1.pdf for standard and user-defined tags.
     * @param value Supported types are String, Char, Integer, Float, byte[], short[]. int[], float[].
     * If value == null, tag is cleared.
     *
     * Byte and Short are allowed but discouraged.  If written to a SAM file, these will be converted to Integer,
     * whereas if written to BAM, getAttribute() will return as Byte or Short, respectively.
     *
     * Long with value between 0 and MAX_UINT is allowed for BAM but discouraged.  Attempting to write such a value
     * to SAM will cause an exception to be thrown.
     *
     * To set unsigned byte[], unsigned short[] or unsigned int[] (which is discouraged because of poor Java language
     * support), setUnsignedArrayAttribute() must be used instead of this method.
     *
     * String values are not validated to ensure that they conform to SAM spec.
     */
    public void setAttribute(final String tag, final Object value) {
        if (value != null && value.getClass().isArray() && Array.getLength(value) == 0) {
            throw new IllegalArgumentException("Empty value passed for tag " + tag);
        }
        setAttribute(SAMTagUtil.getSingleton().makeBinaryTag(tag), value);
    }

    /**
     * Because Java does not support unsigned integer types, we think it is a bad idea to encode them in SAM
     * files.  If you must do so, however, you must call this method rather than setAttribute, because calling
     * this method is the way to indicate that, e.g. a short array should be interpreted as unsigned shorts.
     * @param value must be one of byte[], short[], int[]
     */
    public void setUnsignedArrayAttribute(final String tag, final Object value) {
        if (!value.getClass().isArray()) {
            throw new IllegalArgumentException("Non-array passed to setUnsignedArrayAttribute for tag " + tag);
        }
        if (Array.getLength(value) == 0) {
            throw new IllegalArgumentException("Empty array passed to setUnsignedArrayAttribute for tag " + tag);
        }
        setAttribute(SAMTagUtil.getSingleton().makeBinaryTag(tag), value, true);
    }

    /**
     * @see net.sf.samtools.SAMRecord#setAttribute(java.lang.String, java.lang.Object)
     * @param tag Binary representation of a 2-char String tag as created by SAMTagUtil.
     */
    protected void setAttribute(final short tag, final Object value) {
        setAttribute(tag, value, false);
    }

    protected void setAttribute(final short tag, final Object value, boolean isUnsignedArray) {
        if (value != null &&
                !(value instanceof Byte || value instanceof Short || value instanceof Integer ||
                value instanceof String || value instanceof Character || value instanceof Float ||
                value instanceof byte[] || value instanceof short[] || value instanceof int[] ||
                        value instanceof float[])) {
            throw new SAMException("Attribute type " + value.getClass() + " not supported. Tag: " +
                    SAMTagUtil.getSingleton().makeStringTag(tag));
        }
        if (value == null) {
            if (this.mAttributes != null) this.mAttributes = this.mAttributes.remove(tag);
        }
        else {
            final SAMBinaryTagAndValue tmp;
            if(!isUnsignedArray) {
                tmp = new SAMBinaryTagAndValue(tag, value);
            }
            else {
                if (!value.getClass().isArray() || value instanceof float[]) {
                    throw new SAMException("Attribute type " + value.getClass() +
                            " cannot be encoded as an unsigned array. Tag: " +
                            SAMTagUtil.getSingleton().makeStringTag(tag));
                }
                tmp = new SAMBinaryTagAndUnsignedArrayValue(tag, value);
            }
            if (this.mAttributes == null) this.mAttributes = tmp;
            else this.mAttributes = this.mAttributes.insert(tmp);
        }
    }

    /**
     * Removes all attributes.
     */
    public void clearAttributes() {
        mAttributes = null;
    }

    /**
     * Replace any existing attributes with the given linked item.
     */
    protected void setAttributes(final SAMBinaryTagAndValue attributes) {
        mAttributes = attributes;
    }

    /**
     * @return Pointer to the first of the tags.  Returns null if there are no tags.
     */
    protected SAMBinaryTagAndValue getBinaryAttributes() {
        return mAttributes;
    }

    /**
     * Tag name and value of an attribute, for getAttributes() method.
     */
    public static class SAMTagAndValue {
        public final String tag;
        public final Object value;

        public SAMTagAndValue(final String tag, final Object value) {
            this.tag = tag;
            this.value = value;
        }
    }

    /**
     * @return list of {tag, value} tuples
     */
    public List<SAMTagAndValue> getAttributes() {
        SAMBinaryTagAndValue binaryAttributes = getBinaryAttributes();
        final List<SAMTagAndValue> ret = new ArrayList<SAMTagAndValue>();
        while (binaryAttributes != null) {
            ret.add(new SAMTagAndValue(SAMTagUtil.getSingleton().makeStringTag(binaryAttributes.tag),
                    binaryAttributes.value));
            binaryAttributes = binaryAttributes.getNext();
        }
        return ret;
    }

    Integer getIndexingBin() {
        return mIndexingBin;
    }

    /**
     * Used internally when writing BAMRecords.
     * @param mIndexingBin c.f. http://samtools.sourceforge.net/SAM1.pdf
     */
    void setIndexingBin(final Integer mIndexingBin) {
        this.mIndexingBin = mIndexingBin;
    }

    /**
     * Does not change state of this.
     * @return indexing bin based on alignment start & end.
     */
    int computeIndexingBin() {
        // reg2bin has zero-based, half-open API
        final int alignmentStart = getAlignmentStart()-1;
        int alignmentEnd = getAlignmentEnd();
        if (alignmentEnd <= 0) {
            // If alignment end cannot be determined (e.g. because this read is not really aligned),
            // then treat this as a one base alignment for indexing purposes.
            alignmentEnd = alignmentStart + 1;
        }
        return SAMUtils.reg2bin(alignmentStart, alignmentEnd);
    }

    public SAMFileHeader getHeader() {
        return mHeader;
    }

    /**
     * Setting header into SAMRecord facilitates conversion btw reference sequence names and indices
     * @param header contains sequence dictionary for this SAMRecord
     */
    public void setHeader(final SAMFileHeader header) {
        this.mHeader = header;
    }

    /**
     * If this record has a valid binary representation of the variable-length portion of a binary record stored,
     * return that byte array, otherwise return null.  This will never be true for SAMRecords.  It will be true
     * for BAMRecords that have not been eagerDecoded(), and for which none of the data in the variable-length
     * portion has been changed.
     */
    public byte[] getVariableBinaryRepresentation() {
        return null;
    }

    /**
     * Depending on the concrete implementation, the binary file size of attributes may be known without
     * computing them all.
     * @return binary file size of attribute, if known, else -1
     */
    public int getAttributesBinarySize() {
        return -1;
    }

    /**
     *
     * @return String representation of this.
     * @deprecated This method is not guaranteed to return a valid SAM text representation of the SAMRecord.
     * To get standard SAM text representation, use net.sf.samtools.SAMRecord#getSAMString().
     */
    public String format() {
        final StringBuilder buffer = new StringBuilder();
        addField(buffer, getReadName(), null, null);
        addField(buffer, getFlags(), null, null);
        addField(buffer, getReferenceName(), null, "*");
        addField(buffer, getAlignmentStart(), 0, "*");
        addField(buffer, getMappingQuality(), 0, "0");
        addField(buffer, getCigarString(), null, "*");
        addField(buffer, getMateReferenceName(), null, "*");
        addField(buffer, getMateAlignmentStart(), 0, "*");
        addField(buffer, getInferredInsertSize(), 0, "*");
        addField(buffer, getReadString(), null, "*");
        addField(buffer, getBaseQualityString(), null, "*");
        if (mAttributes != null) {
            SAMBinaryTagAndValue entry = getBinaryAttributes();
            while (entry != null) {
                addField(buffer, formatTagValue(entry.tag, entry.value));
                entry = entry.getNext();
            }
        }
        return buffer.toString();
    }

    private void addField(final StringBuilder buffer, final Object value, final Object defaultValue, final String defaultString) {
        if (safeEquals(value, defaultValue)) {
            addField(buffer, defaultString);
        } else if (value == null) {
            addField(buffer, "");
        } else {
            addField(buffer, value.toString());
        }
    }

    private void addField(final StringBuilder buffer, final String field) {
        if (buffer.length() > 0) {
            buffer.append('\t');
        }
        buffer.append(field);
    }

    private String formatTagValue(final short tag, final Object value) {
        final String tagString = SAMTagUtil.getSingleton().makeStringTag(tag);
        if (value == null || value instanceof String) {
            return tagString + ":Z:" + value;
        } else if (value instanceof Integer || value instanceof Long ||
                   value instanceof Short || value instanceof Byte) {
            return tagString + ":i:" + value;
        } else if (value instanceof Character) {
            return tagString + ":A:" + value;
        } else if (value instanceof Float) {
            return tagString + ":f:" + value;
        } else if (value instanceof byte[]) {
            return tagString + ":H:" + StringUtil.bytesToHexString((byte[]) value);
        } else {
            throw new RuntimeException("Unexpected value type for tag " + tagString +
                                       ": " + value + " of class " + value.getClass().getName());
        }
    }

    private boolean safeEquals(final Object o1, final Object o2) {
        if (o1 == o2) {
            return true;
        } else if (o1 == null || o2 == null) {
            return false;
        } else {
            return o1.equals(o2);
        }
    }

    /**
     * Force all lazily-initialized data members to be initialized.  If a subclass overrides this method,
     * typically it should also call  super method.
     */
    protected void eagerDecode() {
        getCigar();
        getCigarString();
    }

    /**
     * Returns blocks of the read sequence that have been aligned directly to the
     * reference sequence. Note that clipped portions of the read and inserted and
     * deleted bases (vs. the reference) are not represented in the alignment blocks.
     */
    public List<AlignmentBlock> getAlignmentBlocks() {
        if (this.mAlignmentBlocks != null) return this.mAlignmentBlocks;

        final Cigar cigar = getCigar();
        if (cigar == null) return Collections.emptyList();


        final List<AlignmentBlock> alignmentBlocks = new ArrayList<AlignmentBlock>();
        int readBase = 1;
        int refBase  = getAlignmentStart();

        for (final CigarElement e : cigar.getCigarElements()) {
            switch (e.getOperator()) {
                case H : break; // ignore hard clips
                case P : break; // ignore pads
                case S : readBase += e.getLength(); break; // soft clip read bases
                case N : refBase += e.getLength(); break;  // reference skip
                case D : refBase += e.getLength(); break;
                case I : readBase += e.getLength(); break;
                case M :
                case EQ :
                case X :
                    final int length = e.getLength();
                    alignmentBlocks.add(new AlignmentBlock(readBase, refBase, length));
                    readBase += length;
                    refBase  += length;
                    break;
                default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
            }
        }
        this.mAlignmentBlocks = Collections.unmodifiableList(alignmentBlocks);

        return this.mAlignmentBlocks;
    }

    /**
     * Run all validations of CIGAR.  These include validation that the CIGAR makes sense independent of
     * placement, plus validation that CIGAR + placement yields all bases with M operator within the range of the reference.
     * @param recordNumber For error reporting.  -1 if not known.
     * @return List of errors, or null if no errors.
     */
    public List<SAMValidationError> validateCigar(final long recordNumber) {
        List<SAMValidationError> ret = null;

        if (getValidationStringency() != SAMFileReader.ValidationStringency.SILENT && !this.getReadUnmappedFlag()) {
            // Don't know line number, and don't want to force read name to be decoded.
            ret = getCigar().isValid(getReadName(), recordNumber);
            if (getReferenceIndex() != NO_ALIGNMENT_REFERENCE_INDEX) {
                final SAMSequenceRecord sequence = getHeader().getSequence(getReferenceIndex());
                final int referenceSequenceLength = sequence.getSequenceLength();
                for (final AlignmentBlock alignmentBlock : getAlignmentBlocks()) {
                    if (alignmentBlock.getReferenceStart() + alignmentBlock.getLength() - 1 > referenceSequenceLength) {
                        if (ret == null) ret = new ArrayList<SAMValidationError>();
                        ret.add(new SAMValidationError(SAMValidationError.Type.CIGAR_MAPS_OFF_REFERENCE,
                                "CIGAR M operator maps off end of reference", getReadName(), recordNumber));
                        break;
                    }
                }
            }
        }
        return ret;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof SAMRecord)) return false;

        final SAMRecord samRecord = (SAMRecord) o;

        // First check all the elements that do not require decoding
        if (mAlignmentStart != samRecord.mAlignmentStart) return false;
        if (mFlags != samRecord.mFlags) return false;
        if (mInferredInsertSize != samRecord.mInferredInsertSize) return false;
        if (mMappingQuality != samRecord.mMappingQuality) return false;
        if (mMateAlignmentStart != samRecord.mMateAlignmentStart) return false;
        if (mIndexingBin != null ? !mIndexingBin.equals(samRecord.mIndexingBin) : samRecord.mIndexingBin != null)
            return false;
        if (mMateReferenceIndex != null ? !mMateReferenceIndex.equals(samRecord.mMateReferenceIndex) : samRecord.mMateReferenceIndex != null)
            return false;
        if (mReferenceIndex != null ? !mReferenceIndex.equals(samRecord.mReferenceIndex) : samRecord.mReferenceIndex != null)
            return false;

        eagerDecode();
        samRecord.eagerDecode();

        if (mReadName != null ? !mReadName.equals(samRecord.mReadName) : samRecord.mReadName != null) return false;
        if (mAttributes != null ? !mAttributes.equals(samRecord.mAttributes) : samRecord.mAttributes != null)
            return false;
        if (!Arrays.equals(mBaseQualities, samRecord.mBaseQualities)) return false;
        if (mCigar != null ? !mCigar.equals(samRecord.mCigar) : samRecord.mCigar != null)
            return false;
        if (mMateReferenceName != null ? !mMateReferenceName.equals(samRecord.mMateReferenceName) : samRecord.mMateReferenceName != null)
            return false;
        if (!Arrays.equals(mReadBases, samRecord.mReadBases)) return false;
        if (mReferenceName != null ? !mReferenceName.equals(samRecord.mReferenceName) : samRecord.mReferenceName != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        eagerDecode();
        int result = mReadName != null ? mReadName.hashCode() : 0;
        result = 31 * result + (mReadBases != null ? Arrays.hashCode(mReadBases) : 0);
        result = 31 * result + (mBaseQualities != null ? Arrays.hashCode(mBaseQualities) : 0);
        result = 31 * result + (mReferenceName != null ? mReferenceName.hashCode() : 0);
        result = 31 * result + mAlignmentStart;
        result = 31 * result + mMappingQuality;
        result = 31 * result + (mCigarString != null ? mCigarString.hashCode() : 0);
        result = 31 * result + mFlags;
        result = 31 * result + (mMateReferenceName != null ? mMateReferenceName.hashCode() : 0);
        result = 31 * result + mMateAlignmentStart;
        result = 31 * result + mInferredInsertSize;
        result = 31 * result + (mAttributes != null ? mAttributes.hashCode() : 0);
        result = 31 * result + (mReferenceIndex != null ? mReferenceIndex.hashCode() : 0);
        result = 31 * result + (mMateReferenceIndex != null ? mMateReferenceIndex.hashCode() : 0);
        result = 31 * result + (mIndexingBin != null ? mIndexingBin.hashCode() : 0);
        return result;
    }

    /**
     * Perform various validations of SAMRecord.
     * Note that this method deliberately returns null rather than Collections.emptyList() if there
     * are no validation errors, because callers tend to assume that if a non-null list is returned, it is modifiable.
     * @return null if valid.  If invalid, returns a list of error messages.
     */
    public List<SAMValidationError> isValid() {
        // ret is only instantiate if there are errors to report, in order to reduce GC in the typical case
        // in which everything is valid.  It's ugly, but more efficient.
        ArrayList<SAMValidationError> ret = null;
        if (!getReadPairedFlag()) {
            if (getProperPairFlagUnchecked()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_PROPER_PAIR, "Proper pair flag should not be set for unpaired read.", getReadName()));
            }
            if (getMateUnmappedFlagUnchecked()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_MATE_UNMAPPED, "Mate unmapped flag should not be set for unpaired read.", getReadName()));
            }
            if (getMateNegativeStrandFlagUnchecked()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_MATE_NEG_STRAND, "Mate negative strand flag should not be set for unpaired read.", getReadName()));
            }
            if (getFirstOfPairFlagUnchecked()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_FIRST_OF_PAIR, "First of pair flag should not be set for unpaired read.", getReadName()));
            }
            if (getSecondOfPairFlagUnchecked()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_SECOND_OF_PAIR, "Second of pair flag should not be set for unpaired read.", getReadName()));
            }
            if (getMateReferenceIndex() != NO_ALIGNMENT_REFERENCE_INDEX) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_MATE_REF_INDEX, "MRNM should not be set for unpaired read.", getReadName()));
            }
        } else {
            final List<SAMValidationError> errors = isValidReferenceIndexAndPosition(mMateReferenceIndex, mMateReferenceName,
                    getMateAlignmentStart(), true);
            if (errors != null) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.addAll(errors);
            }
            if (!hasMateReferenceName() && !getMateUnmappedFlag()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_MATE_UNMAPPED, "Mapped mate should have mate reference name", getReadName()));
            }
            if (!getFirstOfPairFlagUnchecked() && !getSecondOfPairFlagUnchecked()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND,
                        "Paired read should be marked as first of pair or second of pair.", getReadName()));
            }
/*
            TODO: PIC-97 This validation should be enabled, but probably at this point there are too many
            BAM files that have the proper pair flag set when read or mate is unmapped.
            if (getMateUnmappedFlag() && getProperPairFlagUnchecked()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_PROPER_PAIR, "Proper pair flag should not be set for unpaired read.", getReadName()));
            }
*/
        }
        if (getInferredInsertSize() > MAX_INSERT_SIZE || getInferredInsertSize() < -MAX_INSERT_SIZE) {
            if (ret == null) ret = new ArrayList<SAMValidationError>();
            ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_INSERT_SIZE, "Insert size out of range", getReadName()));
        }
        if (getReadUnmappedFlag()) {
            if (getNotPrimaryAlignmentFlag()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_NOT_PRIM_ALIGNMENT, "Not primary alignment flag should not be set for unmapped read.", getReadName()));
            }
            if (getMappingQuality() != 0) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_MAPPING_QUALITY, "MAPQ should be 0 for unmapped read.", getReadName()));
            }
            /* non-empty CIGAR on unmapped read is now allowed, because there are special reads when SAM is used to store assembly. */
/*
            TODO: PIC-97 This validation should be enabled, but probably at this point there are too many
            BAM files that have the proper pair flag set when read or mate is unmapped.
            if (getProperPairFlagUnchecked()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_PROPER_PAIR, "Proper pair flag should not be set for unmapped read.", getReadName()));
            }
*/
        } else {
            if (getMappingQuality() >= 256) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_MAPPING_QUALITY, "MAPQ should be < 256.", getReadName()));
            }
            if (getCigarLength() == 0) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_CIGAR, "CIGAR should have > zero elements for mapped read.", getReadName()));
            /* todo - will uncomment once unit tests are added
            } else if (getCigar().getReadLength() != getReadLength()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_CIGAR, "CIGAR read length " + getCigar().getReadLength() + " doesn't match read length " + getReadLength(), getReadName()));
            */
            }
            if (getHeader().getSequenceDictionary().size() == 0) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.MISSING_SEQUENCE_DICTIONARY, "Empty sequence dictionary.", getReadName()));
            }
            if (!hasReferenceName()) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_READ_UNMAPPED, "Mapped read should have valid reference name", getReadName()));
            }
/*
            Oops!  We know this is broken in older BAM files, so this having this validation will cause all sorts of
            problems!
            if (getIndexingBin() != null && getIndexingBin() != computeIndexingBin()) {
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_INDEXING_BIN,
                        "Indexing bin (" + getIndexingBin() + ") does not agree with computed value (" + computeIndexingBin() + ")",
                        getReadName()));

            }
*/
        }
        // Validate the RG ID is found in header
        final String rgId = (String)getAttribute(SAMTagUtil.getSingleton().RG);
        if (rgId != null && getHeader().getReadGroup(rgId) == null) {
            if (ret == null) ret = new ArrayList<SAMValidationError>();
            ret.add(new SAMValidationError(SAMValidationError.Type.READ_GROUP_NOT_FOUND,
                    "RG ID on SAMRecord not found in header: " + rgId, getReadName()));
        }
        final List<SAMValidationError> errors = isValidReferenceIndexAndPosition(mReferenceIndex, mReferenceName, getAlignmentStart(), false);
        if (errors != null) {
            if (ret == null) ret = new ArrayList<SAMValidationError>();
            ret.addAll(errors);
        }
        if (this.getReadLength() == 0 && !this.getNotPrimaryAlignmentFlag()) {
            Object fz = getAttribute(SAMTagUtil.getSingleton().FZ);
            if (fz == null) {
                String cq = (String)getAttribute(SAMTagUtil.getSingleton().CQ);
                String cs = (String)getAttribute(SAMTagUtil.getSingleton().CS);
                if (cq == null || cq.length() == 0 || cs == null || cs.length() == 0) {
                    if (ret == null) ret = new ArrayList<SAMValidationError>();
                    ret.add(new SAMValidationError(SAMValidationError.Type.EMPTY_READ,
                            "Zero-length read without FZ, CS or CQ tag", getReadName()));
                } else if (!getReadUnmappedFlag()) {
                    boolean hasIndel = false;
                    for (CigarElement cigarElement : getCigar().getCigarElements()) {
                        if (cigarElement.getOperator() == CigarOperator.DELETION ||
                                cigarElement.getOperator() == CigarOperator.INSERTION) {
                            hasIndel = true;
                            break;
                        }
                    }
                    if (!hasIndel) {
                        if (ret == null) ret = new ArrayList<SAMValidationError>();
                        ret.add(new SAMValidationError(SAMValidationError.Type.EMPTY_READ,
                                "Colorspace read with zero-length bases but no indel", getReadName()));
                    }
                }
            }
        }
        if (this.getReadLength() != getBaseQualities().length &&  !Arrays.equals(getBaseQualities(), NULL_QUALS)) {
            if (ret == null) ret = new ArrayList<SAMValidationError>();
            ret.add(new SAMValidationError(SAMValidationError.Type.MISMATCH_READ_LENGTH_AND_QUALS_LENGTH,
                    "Read length does not match quals length", getReadName()));
        }
        if (ret == null || ret.size() == 0) {
            return null;
        }
        return ret;
    }

    /**
     * Gets the source of this SAM record -- both the reader that retrieved the record and the position on disk from
     * whence it came. 
     * @return The file source.  Note that the reader will be null if not activated using SAMFileReader.enableFileSource().
     */
    public SAMFileSource getFileSource() {
        return mFileSource;
    }

    /**
     * Sets a marker providing the source reader for this file and the position in the file from which the read originated.
     * @param fileSource source of the given file.
     */
    protected void setFileSource(final SAMFileSource fileSource) {
        mFileSource = fileSource;
    }

    private List<SAMValidationError> isValidReferenceIndexAndPosition(final Integer referenceIndex, final String referenceName,
                                                          final int alignmentStart, final boolean isMate) {
        final boolean hasReference = hasReferenceName(referenceIndex, referenceName);

        // ret is only instantiate if there are errors to report, in order to reduce GC in the typical case
        // in which everything is valid.  It's ugly, but more efficient.
        ArrayList<SAMValidationError> ret = null;
        if (!hasReference) {
            if (alignmentStart != 0) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_ALIGNMENT_START, buildMessage("Alignment start should be 0 because reference name = *.", isMate), getReadName()));
            }
        } else {
            if (alignmentStart == 0) {
                if (ret == null) ret = new ArrayList<SAMValidationError>();
                ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_ALIGNMENT_START, buildMessage("Alignment start should != 0 because reference name != *.", isMate), getReadName()));
            }

            if (getHeader().getSequenceDictionary().size() > 0) {
                final SAMSequenceRecord sequence =
                        (referenceIndex != null? getHeader().getSequence(referenceIndex): getHeader().getSequence(referenceName));
                if (sequence == null) {
                    if (ret == null) ret = new ArrayList<SAMValidationError>();
                    ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_REFERENCE_INDEX, buildMessage("Reference sequence not found in sequence dictionary.", isMate), getReadName()));
                } else {
                    if (alignmentStart > sequence.getSequenceLength()) {
                        if (ret == null) ret = new ArrayList<SAMValidationError>();
                        ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_ALIGNMENT_START, buildMessage("Alignment start (" + alignmentStart + ") must be <= reference sequence length (" +
                                sequence.getSequenceLength() + ") on reference " + sequence.getSequenceName(), isMate), getReadName()));
                    }
                }
            }
        }
        return ret;
    }
    
    private String buildMessage(final String baseMessage, final boolean isMate) {
        return isMate ? "Mate " + baseMessage : baseMessage;
    }

    /**
     * Note that this does a shallow copy of everything, except for the attribute list, for which a copy of the list
     * is made, but the attributes themselves are copied by reference.  This should be safe because callers should
     * never modify a mutable value returned by any of the get() methods anyway.
     */
    @Override
    public Object clone() throws CloneNotSupportedException {
        final SAMRecord newRecord = (SAMRecord)super.clone();
        if (mAttributes != null) {
            newRecord.mAttributes = this.mAttributes.copy();
        }

        return newRecord;
    }

    /** Simple toString() that gives a little bit of useful info about the read. */
    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder(64);
        builder.append(getReadName());
        if (getReadPairedFlag()) {
            if (getFirstOfPairFlag()) {
                builder.append(" 1/2");
            }
            else {
                builder.append(" 2/2");
            }
        }

        builder.append(" ");
        builder.append(String.valueOf(getReadLength()));
        builder.append("b");

        if (getReadUnmappedFlag()) {
            builder.append(" unmapped read.");
        }
        else {
            builder.append(" aligned read.");
        }

        return builder.toString();
    }

    /**
	Returns the record in the SAM line-based text format.  Fields are
	separated by '\t' characters, and the String is terminated by '\n'.
    */
    public String getSAMString() {
	return SAMTextWriter.getSAMString(this);
    }
}

