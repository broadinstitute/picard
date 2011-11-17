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

import net.sf.samtools.util.StringUtil;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;


/**
 * Wrapper class for binary BAM records.
 * Delays unpacking all data binary until requested.
 */
public class BAMRecord extends SAMRecord {
    /**
     * Offset of the read name in the variable length section of the disk representation of BAMRecord
     */
    private static final int READ_NAME_OFFSET = 0;

    /**
     * Variable-length part of BAMRecord.  Lazily decoded.
     */
    private byte[] mRestOfBinaryData = null;

    // Various lengths are stored, because they are in the fixed-length part of the BAMRecord, and it is
    // more efficient to remember them than decode the element they store the length of.
    // The length becomes invalid if the element is changed with a set() method.
    private int mReadLength = 0;
    private boolean mReadLengthValid = true;
    private final short mReadNameLength;
    private boolean mReadNameLengthValid = true;
    private final int mCigarLength;
    private boolean mCigarLengthValid = true;

    // Whether or not the getter needs to decode the corresponding element.
    // For all the other variable length elements, null == not yet decoded.
    private boolean mAttributesDecoded = false;
    private boolean mCigarDecoded = false;

    /**
     * If any of the properties set from mRestOfBinaryData have been overridden by calls to setters,
     * this is set to true, indicating that mRestOfBinaryData cannot be used to write this record to disk.
     */
    private boolean mBinaryDataStale;

    protected BAMRecord(final SAMFileHeader header,
                        final int referenceID,
                        final int coordinate,
                        final short readNameLength,
                        final short mappingQuality,
                        final int indexingBin,
                        final int cigarLen,
                        final int flags,
                        final int readLen,
                        final int mateReferenceID,
                        final int mateCoordinate,
                        final int insertSize,
                        final byte[] restOfData) {
        super(header);
        setReferenceIndex(referenceID);
        setAlignmentStart(coordinate);
        mReadNameLength = readNameLength;
        setMappingQuality(mappingQuality);
        mCigarLength = cigarLen;
        setFlags(flags);
        mReadLength = readLen;
        setMateReferenceIndex(mateReferenceID);
        setMateAlignmentStart(mateCoordinate);
        setInferredInsertSize(insertSize);
        mRestOfBinaryData = restOfData;

        // Set these to null in order to mark them as being candidates for lazy initialization.
        // If this is not done, they will have non-null defaults.
        super.setReadName(null);
        super.setCigarString(null);
        super.setReadBases(null);
        super.setBaseQualities(null);

        // Do this after the above because setCigarString will clear it.
        setIndexingBin(indexingBin);

        // Mark the binary block as being valid for writing back out to disk
        mBinaryDataStale = false;
    }

    /**
     * Force all the lazily-initialized attributes to be decoded.
     */
    protected void eagerDecode() {
        getReadName();
        getCigar();
        getReadBases();
        getBaseQualities();
        getBinaryAttributes();
        super.eagerDecode();
        mRestOfBinaryData = null;
    }

    /**
     * If this record has a valid binary representation of the variable-length portion of a binary record stored,
     * return that byte array, otherwise return null.  This will never be true for SAMRecords.  It will be true
     * for BAMRecords that have not been eagerDecoded(), and for which none of the data in the variable-length
     * portion has been changed.
     */
    @Override
    public byte[] getVariableBinaryRepresentation() {
        if (mBinaryDataStale) {
            return null;
        }
        // This may have been set to null by eagerDecode()
        return mRestOfBinaryData;
    }

    /**
     * Depending on the concrete implementation, the binary file size of attributes may be known without
     * computing them all.
     *
     * @return binary file size of attribute, if known, else -1.
     */
    @Override
    public int getAttributesBinarySize() {
        if (mBinaryDataStale || mRestOfBinaryData == null) {
            return -1;
        }
        final int tagsOffset = readNameSize() + cigarSize() + basesSize() + qualsSize();
        return mRestOfBinaryData.length - tagsOffset;
    }

    @Override
    public void setReadName(final String value) {
        super.setReadName(value);
        mBinaryDataStale = true;
        mReadNameLengthValid = false;
    }

    @Override
    public void setCigar(final Cigar cigar) {
        super.setCigar(cigar);
        mBinaryDataStale = true;
        mCigarLengthValid = false;
        mCigarDecoded = true;
    }

    @Override
    public void setCigarString(final String value) {
        super.setCigarString(value);
        mBinaryDataStale = true;
        mCigarLengthValid = false;
        mCigarDecoded = true;
    }

    @Override
    public void setReadBases(final byte[] value) {
        super.setReadBases(value);
        mBinaryDataStale = true;
        mReadLengthValid = false;
    }

    @Override
    public void setBaseQualities(final byte[] value) {
        super.setBaseQualities(value);
        mBinaryDataStale = true;
    }

    @Override
    protected void setAttribute(final short tag, final Object value, final boolean isUnsignedArray) {
        // populate all the attributes from the binary block before overwriting one
        getBinaryAttributes();
        super.setAttribute(tag, value, isUnsignedArray);
        mBinaryDataStale = true;
    }

    /**
     * Removes all attributes.
     */
    @Override
    public void clearAttributes() {
        mAttributesDecoded = true;
        mBinaryDataStale = true;
        super.clearAttributes();
    }

    /**
     * Avoids decoding binary block to get read length.
     */
    @Override
    public int getReadLength() {
        if (mReadLengthValid) {
            return mReadLength;
        }
        return super.getReadLength();
    }

    @Override
    public String getReadName() {
        String result = super.getReadName();
        if (mRestOfBinaryData != null && result == null) {
            result = decodeReadName();
            super.setReadName(result);
        }
        return result;
    }

    /**
     * Avoids decoding read name to get read name length.  Do not include null terminator.
     */
    @Override
    public int getReadNameLength() {
        if (mReadNameLengthValid) {
            return mReadNameLength - 1;
        }
        return super.getReadNameLength();
    }

    @Override
    public Cigar getCigar() {
        if (mRestOfBinaryData != null && !mCigarDecoded) {
            final int cigarOffset = readNameSize();
            final ByteBuffer byteBuffer  = ByteBuffer.wrap(mRestOfBinaryData, cigarOffset, cigarSize());
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
            super.initializeCigar(BinaryCigarCodec.getSingleton().decode(byteBuffer));
            mCigarDecoded = true;
            if (getValidationStringency() != SAMFileReader.ValidationStringency.SILENT && !this.getReadUnmappedFlag()) {
                // Don't know line number, and don't want to force read name to be decoded.
                SAMUtils.processValidationErrors(validateCigar(-1L), -1, getValidationStringency());
            }
        }
        return super.getCigar();
    }

    /**
     * Avoids decoding CIGAR in order to get length.
     */
    @Override
    public int getCigarLength() {
        if (mCigarLengthValid) {
            return mCigarLength;
        } else {
            return super.getCigarLength();
        }
    }

    @Override
    public byte[] getReadBases() {
        byte[] result = super.getReadBases();
        if (mRestOfBinaryData != null && result == null) {
            result = decodeReadBases();
            super.setReadBases(result);
        }
        return result;
    }

    @Override
    public byte[] getBaseQualities() {
        byte[] ret = super.getBaseQualities();
        if (mRestOfBinaryData != null && ret == null) {
            ret = decodeBaseQualities();
            super.setBaseQualities(ret);
        }
        return ret;
    }

    @Override
    public Object getAttribute(final short tag) {
        if (!mAttributesDecoded) {
            decodeAttributes();
        }
        return super.getAttribute(tag);
    }

    @Override
    protected SAMBinaryTagAndValue getBinaryAttributes() {
        if (!mAttributesDecoded) {
            decodeAttributes();
        }
        return super.getBinaryAttributes();
    }

    private void decodeAttributes() {
        if (mAttributesDecoded) {
            return;
        }
        mAttributesDecoded = true;
        final int tagsOffset = readNameSize() + cigarSize() + basesSize() + qualsSize();
        final int tagsSize = mRestOfBinaryData.length - tagsOffset;
        final SAMBinaryTagAndValue attributes = BinaryTagCodec.readTags(mRestOfBinaryData, tagsOffset, tagsSize, getValidationStringency());
        setAttributes(attributes);
    }

    private byte[] decodeBaseQualities() {
        if (mReadLength == 0) {
            return SAMRecord.NULL_QUALS;
        }
        final int qualsOffset = readNameSize() + cigarSize() + basesSize();
        final byte[] ret = new byte[qualsSize()];
        System.arraycopy(mRestOfBinaryData, qualsOffset, ret, 0, qualsSize());
        if (ret.length > 0 && ret[0] == (byte) 0xFF) {
            // BAM files store missing qualities as an array of 0xFF bytes.
            // 0xFF is an illegal quality score value (it cannot be encoded in SAM)
            // and so the first byte is a suitable marker.
            // We hide this quirk of the BAM encoding so that the BAM interface looks the same as SAM.
            return NULL_QUALS;
        }
        return ret;
    }

    private String decodeReadName() {
        // Don't include terminating null
        return StringUtil.bytesToString(mRestOfBinaryData, READ_NAME_OFFSET, mReadNameLength-1);
    }

    private byte[] decodeReadBases() {
        if (mReadLength == 0) {
            return NULL_SEQUENCE;
        }
        final int basesOffset = readNameSize() + cigarSize();
        return SAMUtils.compressedBasesToBytes(mReadLength, mRestOfBinaryData, basesOffset);
    }

    /* methods for computing disk size of variably-sized elements, in order to locate
     * elements in mRestOfBinaryData */

    private int readNameSize() {
        return mReadNameLength;
    }

    private int cigarSize() {
        return mCigarLength * 4;
    }

    private int basesSize() {
        return (mReadLength + 1)/2;
    }

    private int qualsSize() {
        return mReadLength;
    }
}
