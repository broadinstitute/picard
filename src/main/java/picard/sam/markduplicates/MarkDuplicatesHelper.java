package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;

public interface MarkDuplicatesHelper {

    void generateDuplicateIndexes(final boolean useBarcodes, final boolean indexOpticalDuplicates);
    ReadEndsForMarkDuplicates buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec, final boolean useBarcodes);
    short getReadDuplicateScore(final SAMRecord rec, final ReadEndsForMarkDuplicates pairedEnds);
}
