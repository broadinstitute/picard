package net.sf.samtools;

/**
 * Default factory for creating SAM and BAM records used by the SAMFileReader classes.
 *
 * @author Tim Fennell
 */
public class DefaultSAMRecordFactory implements SAMRecordFactory {

    /** Create a new SAMRecord to be filled in */
    public SAMRecord createSAMRecord(final SAMFileHeader header) {
        return new SAMRecord(header);
    }

    /** Create a new BAM Record. */
    public BAMRecord createBAMRecord (final SAMFileHeader header,
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
                                      final byte[] variableLengthBlock) {

        return new BAMRecord(header,
                             referenceSequenceIndex,
                             alignmentStart,
                             readNameLength,
                             mappingQuality,
                             indexingBin,
                             cigarLen,
                             flags,
                             readLen,
                             mateReferenceSequenceIndex,
                             mateAlignmentStart,
                             insertSize,
                             variableLengthBlock);
    }
}
