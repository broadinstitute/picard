package net.sf.samtools;

/**
 * Factory interface which allows plugging in of different classes for generating instances of
 * SAMRecord and BAMRecord when reading from SAM/BAM files.
 *
 * @author Tim Fennell
 */
public interface SAMRecordFactory {

    /** Create a new SAMRecord to be filled in */
    public SAMRecord createSAMRecord(SAMFileHeader header);

    /** Create a new BAM Record. */
    public BAMRecord createBAMRecord(final SAMFileHeader header,
                                     final int referenceSequenceIndex,
                                     final int alignmentStart,
                                     final short readNameLength,
                                     final short mappingQuality,
                                     final int indexingBin,
                                     final int cigarLen,
                                     final int flags,
                                     final int readLen,
                                     final int mateReferenceSequenceIndex,
                                     final int mateAlignmentStart,
                                     final int insertSize,
                                     final byte[] variableLengthBlock);
}
