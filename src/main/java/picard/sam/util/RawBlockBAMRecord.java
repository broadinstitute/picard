package picard.sam.util;

import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMBinaryTagAndValue;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringUtil;

import java.util.Arrays;

/**
 * A BAMRecord that keeps the variable-length block it was decoded from until its read name, cigar, bases or
 * qualities are set, so that {@link ParallelBAMFileWriter} can copy those four fields from the block rather than
 * decode and re-encode them. Setting tags keeps the block: the writer encodes the tags from the record's tag list.
 * The arrays getReadBases and getBaseQualities return must not be changed in place.
 *
 * After {@link #deferTags()} the record's tags and cigar stay undecoded: getAttribute answers from the block, and
 * setAttribute with a null or String value is recorded as an edit, which the writer applies as it copies the tags.
 * Anything else that reads or changes the tags decodes them, and the cigar, as {@link BAMRecord} would have at SILENT
 * before the edits, then applies the edits in order.
 */
public class RawBlockBAMRecord extends BAMRecord {
    /** Makes a RawBlockBAMRecord for every BAM record and what DefaultSAMRecordFactory makes for every other. */
    public static final SAMRecordFactory FACTORY = new SAMRecordFactory() {
        @Override
        public SAMRecord createSAMRecord(final SAMFileHeader header) {
            return DefaultSAMRecordFactory.getInstance().createSAMRecord(header);
        }

        @Override
        public BAMRecord createBAMRecord(final SAMFileHeader header, final int referenceSequenceIndex,
                                         final int alignmentStart, final short readNameLength,
                                         final short mappingQuality, final int indexingBin, final int cigarLen,
                                         final int flags, final int readLen, final int mateReferenceSequenceIndex,
                                         final int mateAlignmentStart, final int insertSize,
                                         final byte[] variableLengthBlock) {
            return new RawBlockBAMRecord(header, referenceSequenceIndex, alignmentStart, readNameLength,
                    mappingQuality, indexingBin, cigarLen, flags, readLen, mateReferenceSequenceIndex,
                    mateAlignmentStart, insertSize, variableLengthBlock);
        }
    };

    private static final short CG = SAMTag.CG.getBinaryTag();
    // Largest BAM cigar operator code, X; BinaryCigarCodec throws on a larger one.
    private static final int MAX_CIGAR_OPERATOR = 8;

    private final byte[] bytes;
    private byte[] block;

    // While the tags are deferred: each tag's start and end offsets in bytes, in signed tag-code order, which is the
    // order BinaryTagCodec.readTags gives them. Null otherwise.
    private int[] tags = null;
    private int tagCount = 0;
    private int tagsStart = 0;
    // The edits, at most one per tag code, in signed tag-code order: the value of the last setAttribute for that code,
    // null for a removal.
    private short[] editCodes = null;
    private String[] editValues = null;
    private int editCount = 0;
    // Whether the cigar is still to be decoded at SILENT, as the decoding threads decode a validated record's cigar.
    private boolean silentCigar = false;

    protected RawBlockBAMRecord(final SAMFileHeader header, final int referenceSequenceIndex, final int alignmentStart,
                                final short readNameLength, final short mappingQuality, final int indexingBin,
                                final int cigarLen, final int flags, final int readLen,
                                final int mateReferenceSequenceIndex, final int mateAlignmentStart,
                                final int insertSize, final byte[] variableLengthBlock) {
        super(header, referenceSequenceIndex, alignmentStart, readNameLength, mappingQuality, indexingBin, cigarLen,
                flags, readLen, mateReferenceSequenceIndex, mateAlignmentStart, insertSize, variableLengthBlock);
        this.bytes = variableLengthBlock;
        this.block = variableLengthBlock;
    }

    /** The block this record was decoded from, or null once its read name, cigar, bases or qualities were set. */
    byte[] getOriginalBlock() {
        return block;
    }

    /** The record's tags, in the order BAMRecordCodec writes them. */
    SAMBinaryTagAndValue getTags() {
        return getBinaryAttributes();
    }

    /**
     * For a record whose tags and cigar would otherwise be decoded at SILENT now, before anything else reads or
     * changes it: leaves them undecoded and returns true when decoding them cannot fail: every cigar operator is
     * known, and the tags fill the rest of the block exactly, with no CG tag, no tag code twice, no H value and no B
     * subtype other than c, C, s, S, i, I and f. Otherwise returns false and changes nothing.
     */
    boolean deferTags() {
        if (tags != null) {
            return true;
        }
        if (block == null || bytes == null) {
            return false;
        }
        final int readLength = getReadLength();
        final int cigarLength = getCigarLength();
        final int cigarOffset = getReadNameLength() + 1;
        final long start = cigarOffset + 4L * cigarLength + (readLength + 1L) / 2 + readLength;
        if (readLength < 0 || cigarOffset < 1 || start > bytes.length) {
            return false;
        }
        for (int i = 0, position = cigarOffset; i < cigarLength; i++, position += 4) {
            if ((bytes[position] & 0xF) > MAX_CIGAR_OPERATOR) {
                return false;
            }
        }
        final int end = bytes.length;
        int[] found = new int[32];
        int count = 0;
        int position = (int) start;
        while (position < end) {
            if (end - position < 3) {
                return false;
            }
            final short code = codeAt(bytes, position);
            if (code == CG) {
                return false;
            }
            final int value = position + 3;
            final long next;
            switch (bytes[position + 2]) {
                case 'A':
                case 'c':
                case 'C':
                    next = value + 1L;
                    break;
                case 's':
                case 'S':
                    next = value + 2L;
                    break;
                case 'i':
                case 'I':
                case 'f':
                    next = value + 4L;
                    break;
                case 'Z': {
                    int nul = value;
                    while (nul < end && bytes[nul] != 0) {
                        nul++;
                    }
                    next = nul + 1L;
                    break;
                }
                case 'B': {
                    if (end - value < 5) {
                        return false;
                    }
                    final int size;
                    switch (bytes[value]) {
                        case 'c':
                        case 'C':
                            size = 1;
                            break;
                        case 's':
                        case 'S':
                            size = 2;
                            break;
                        case 'i':
                        case 'I':
                        case 'f':
                            size = 4;
                            break;
                        default:
                            return false;
                    }
                    final int elements = intAt(bytes, value + 1);
                    if (elements < 0) {
                        return false;
                    }
                    next = value + 5L + (long) elements * size;
                    break;
                }
                default:
                    return false;
            }
            if (next > end) {
                return false;
            }
            int slot = count;
            while (slot > 0 && codeAt(bytes, found[2 * slot - 2]) > code) {
                slot--;
            }
            if (slot > 0 && codeAt(bytes, found[2 * slot - 2]) == code) {
                return false;
            }
            if (2 * count + 2 > found.length) {
                found = Arrays.copyOf(found, 2 * found.length);
            }
            System.arraycopy(found, 2 * slot, found, 2 * slot + 2, 2 * (count - slot));
            found[2 * slot] = position;
            found[2 * slot + 1] = (int) next;
            count++;
            position = (int) next;
        }
        tags = found;
        tagCount = count;
        tagsStart = (int) start;
        silentCigar = true;
        return true;
    }

    /** Whether the tags are deferred, so that {@link #writeTags} writes them. */
    boolean hasDeferredTags() {
        return tags != null;
    }

    /** The most bytes {@link #writeTags} writes. */
    int maxTagsLength() {
        int length = bytes.length - tagsStart;
        for (int i = 0; i < editCount; i++) {
            if (editValues[i] != null) {
                length += 4 + editValues[i].length();
            }
        }
        return length;
    }

    /**
     * Writes the deferred tags with their edits at position in out, as BinaryTagCodec.writeTag writes the tag list
     * BAMRecordCodec would encode, and returns the position after them.
     */
    int writeTags(final byte[] out, int position) {
        int edit = 0;
        for (int i = 0; i < tagCount; i++) {
            final int start = tags[2 * i];
            final short code = codeAt(bytes, start);
            while (edit < editCount && editCodes[edit] < code) {
                position = writeEdit(out, position, edit++);
            }
            if (edit < editCount && editCodes[edit] == code) {
                position = writeEdit(out, position, edit++);
            } else {
                position = writeTag(out, position, start, tags[2 * i + 1]);
            }
        }
        while (edit < editCount) {
            position = writeEdit(out, position, edit++);
        }
        return position;
    }

    private int writeEdit(final byte[] out, final int position, final int edit) {
        final String value = editValues[edit];
        if (value == null) {
            return position;
        }
        out[position] = (byte) editCodes[edit];
        out[position + 1] = (byte) (editCodes[edit] >> 8);
        out[position + 2] = 'Z';
        final int length = value.length();
        // StringUtil.stringToBytes: each char's low byte.
        for (int i = 0; i < length; i++) {
            out[position + 3 + i] = (byte) value.charAt(i);
        }
        out[position + 3 + length] = 0;
        return position + 4 + length;
    }

    /** One tag from the block: integers in the type BinaryTagCodec.getIntegerType gives their value, others as they are. */
    private int writeTag(final byte[] out, final int position, final int start, final int end) {
        final int value = start + 3;
        final long integer;
        switch (bytes[start + 2]) {
            case 'c':
                integer = bytes[value];
                break;
            case 'C':
                integer = bytes[value] & 0xFF;
                break;
            case 's':
                integer = (short) ushortAt(bytes, value);
                break;
            case 'S':
                integer = ushortAt(bytes, value);
                break;
            case 'i':
                integer = intAt(bytes, value);
                break;
            case 'I':
                integer = intAt(bytes, value) & 0xFFFFFFFFL;
                break;
            default:
                System.arraycopy(bytes, start, out, position, end - start);
                return position + end - start;
        }
        out[position] = bytes[start];
        out[position + 1] = bytes[start + 1];
        final int width;
        if (integer > Integer.MAX_VALUE) {
            out[position + 2] = 'I';
            width = 4;
        } else if (integer > 65535) {
            out[position + 2] = 'i';
            width = 4;
        } else if (integer > 32767) {
            out[position + 2] = 'S';
            width = 2;
        } else if (integer > 255) {
            out[position + 2] = 's';
            width = 2;
        } else if (integer > 127) {
            out[position + 2] = 'C';
            width = 1;
        } else if (integer >= -128) {
            out[position + 2] = 'c';
            width = 1;
        } else if (integer >= -32768) {
            out[position + 2] = 's';
            width = 2;
        } else {
            out[position + 2] = 'i';
            width = 4;
        }
        for (int i = 0; i < width; i++) {
            out[position + 3 + i] = (byte) (integer >>> (8 * i));
        }
        return position + 3 + width;
    }

    /** Decodes the deferred tags and cigar at SILENT, as the decoding threads would have, then applies the edits. */
    private void decodeDeferred() {
        if (tags == null) {
            return;
        }
        tags = null;
        final ValidationStringency stringency = getValidationStringency();
        setValidationStringency(ValidationStringency.SILENT);
        try {
            super.getBinaryAttributes();
            if (silentCigar) {
                super.getCigar();
                silentCigar = false;
            }
        } finally {
            setValidationStringency(stringency);
        }
        for (int i = 0; i < editCount; i++) {
            super.setAttribute(editCodes[i], editValues[i], false);
        }
        editCodes = null;
        editValues = null;
        editCount = 0;
    }

    @Override
    public Object getAttribute(final short tag) {
        if (tags == null) {
            return super.getAttribute(tag);
        }
        for (int i = 0; i < editCount; i++) {
            if (editCodes[i] == tag) {
                return editValues[i];
            }
        }
        for (int i = 0; i < tagCount; i++) {
            final int start = tags[2 * i];
            if (codeAt(bytes, start) != tag) {
                continue;
            }
            // BinaryTagCodec.readSingleValue's values.
            final int value = start + 3;
            switch (bytes[start + 2]) {
                case 'A':
                    return (char) bytes[value];
                case 'c':
                    return (int) bytes[value];
                case 'C':
                    return bytes[value] & 0xFF;
                case 's':
                    return (int) (short) ushortAt(bytes, value);
                case 'S':
                    return ushortAt(bytes, value);
                case 'i':
                    return intAt(bytes, value);
                case 'I': {
                    final long unsigned = intAt(bytes, value) & 0xFFFFFFFFL;
                    if (unsigned <= Integer.MAX_VALUE) {
                        return (int) unsigned;
                    }
                    return unsigned;
                }
                case 'f':
                    return Float.intBitsToFloat(intAt(bytes, value));
                case 'Z':
                    return StringUtil.bytesToString(bytes, value, tags[2 * i + 1] - 1 - value);
                default:
                    decodeDeferred();
                    return super.getAttribute(tag);
            }
        }
        return null;
    }

    @Override
    protected void setAttribute(final short tag, final Object value, final boolean isUnsignedArray) {
        if (tags != null && !isUnsignedArray && tag != CG && (value == null || value instanceof String)) {
            int i = 0;
            while (i < editCount && editCodes[i] < tag) {
                i++;
            }
            if (i == editCount || editCodes[i] != tag) {
                if (editCodes == null) {
                    editCodes = new short[4];
                    editValues = new String[4];
                } else if (editCount == editCodes.length) {
                    editCodes = Arrays.copyOf(editCodes, 2 * editCount);
                    editValues = Arrays.copyOf(editValues, 2 * editCount);
                }
                System.arraycopy(editCodes, i, editCodes, i + 1, editCount - i);
                System.arraycopy(editValues, i, editValues, i + 1, editCount - i);
                editCodes[i] = tag;
                editCount++;
            }
            editValues[i] = (String) value;
            return;
        }
        decodeDeferred();
        super.setAttribute(tag, value, isUnsignedArray);
    }

    @Override
    protected SAMBinaryTagAndValue getBinaryAttributes() {
        decodeDeferred();
        return super.getBinaryAttributes();
    }

    @Override
    protected void setAttributes(final SAMBinaryTagAndValue attributes) {
        decodeDeferred();
        super.setAttributes(attributes);
    }

    @Override
    public void clearAttributes() {
        decodeDeferred();
        super.clearAttributes();
    }

    @Override
    public byte[] getVariableBinaryRepresentation() {
        return tags != null && editCount > 0 ? null : super.getVariableBinaryRepresentation();
    }

    @Override
    public int getAttributesBinarySize() {
        if (editCount > 0) {
            decodeDeferred();
        }
        return super.getAttributesBinarySize();
    }

    @Override
    public boolean isUnsignedArrayAttribute(final String tag) {
        decodeDeferred();
        return super.isUnsignedArrayAttribute(tag);
    }

    @Override
    public Cigar getCigar() {
        if (!silentCigar) {
            return super.getCigar();
        }
        final ValidationStringency stringency = getValidationStringency();
        setValidationStringency(ValidationStringency.SILENT);
        try {
            final Cigar cigar = super.getCigar();
            silentCigar = false;
            return cigar;
        } finally {
            setValidationStringency(stringency);
        }
    }

    // SAMRecord.clone copies the tag list without decoding it; equals and hashCode call eagerDecode first.
    @Override
    public boolean equals(final Object o) {
        decodeDeferred();
        if (o instanceof RawBlockBAMRecord) {
            ((RawBlockBAMRecord) o).decodeDeferred();
        }
        return super.equals(o);
    }

    @Override
    public int hashCode() {
        decodeDeferred();
        return super.hashCode();
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        decodeDeferred();
        return super.clone();
    }

    @Override
    public void setReadName(final String value) {
        block = null;
        super.setReadName(value);
    }

    @Override
    public void setCigar(final Cigar cigar) {
        block = null;
        silentCigar = false;
        super.setCigar(cigar);
    }

    @Override
    public void setCigarString(final String value) {
        block = null;
        silentCigar = false;
        super.setCigarString(value);
    }

    @Override
    public void setReadBases(final byte[] value) {
        block = null;
        super.setReadBases(value);
    }

    // SAMRecord.setReadString("*") assigns the bases without calling setReadBases.
    @Override
    public void setReadString(final String value) {
        block = null;
        super.setReadString(value);
    }

    @Override
    public void setBaseQualities(final byte[] value) {
        block = null;
        super.setBaseQualities(value);
    }

    private static short codeAt(final byte[] bytes, final int offset) {
        return (short) ((bytes[offset] & 0xFF) | bytes[offset + 1] << 8);
    }

    private static int ushortAt(final byte[] bytes, final int offset) {
        return (bytes[offset] & 0xFF) | (bytes[offset + 1] & 0xFF) << 8;
    }

    private static int intAt(final byte[] bytes, final int offset) {
        return (bytes[offset] & 0xFF) | (bytes[offset + 1] & 0xFF) << 8 | (bytes[offset + 2] & 0xFF) << 16 |
                (bytes[offset + 3] & 0xFF) << 24;
    }
}
