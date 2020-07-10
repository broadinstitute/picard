package picard.fastq;

import java.io.File;

/**
 * Deprecated since July 2020. Use SamToBfqWriter instead
 */
@Deprecated
public class BamToBfqWriter extends SamToBfqWriter {
    public BamToBfqWriter(final File bamFile, final File referenceSequence, final String outputPrefix, final Integer total, final Integer chunk, final boolean pairedReads, final String namePrefix, final boolean includeNonPfReads, final boolean clipAdapters, final Integer basesToWrite) {
        super(bamFile, referenceSequence, outputPrefix, total, chunk, pairedReads, namePrefix, includeNonPfReads, clipAdapters, basesToWrite);
    }

    public BamToBfqWriter(final File bamFile, final String outputPrefix, final boolean pairedReads, final String namePrefix, final boolean includeNonPfReads) {
        super(bamFile, outputPrefix, pairedReads, namePrefix, includeNonPfReads);
    }
}
