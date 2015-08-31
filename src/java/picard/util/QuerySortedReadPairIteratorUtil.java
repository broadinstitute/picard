package picard.util;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.PeekableIterator;

/**
 * A collection of helper utilities for iterating through reads that are in query-name sorted
 * read order as pairs
 */
public class QuerySortedReadPairIteratorUtil {

    public static class ReadPair {
        public SAMRecord read1 = null;
        public SAMRecord read2 = null;
    }

    /**
     * Get the next read pair (where both have the same read name).
     * If we encounter an unpaired read, the second read in the pair will be set to null.
     *
     * @param iterator iterator of reads
     * @return ReadPair object holding the reads, or null if there are no more reads in the iterator
     */
    public static ReadPair getNextReadPair(final PeekableIterator<SAMRecord> iterator) {

        final ReadPair readPair = new ReadPair();
        readPair.read1 = getNextUsableRead(iterator, false);
        if (readPair.read1 == null) {
            return null;
        }

        final SAMRecord peekedNextRead = getNextUsableRead(iterator, true);
        if (peekedNextRead != null && peekedNextRead.getReadName().equals(readPair.read1.getReadName())) {
            readPair.read2 = getNextUsableRead(iterator, false);
        }

        return readPair;
    }

    /**
     * Return the next usable read in the iterator
     *
     * @param iterator the iterator to pull from
     * @param justPeek if true, just peek the next usable read rather than pulling it (note: it may remove unusable reads from the iterator)
     * @return the next read or null if none are left
     */
    private static SAMRecord getNextUsableRead(final PeekableIterator<SAMRecord> iterator, final boolean justPeek) {

        while (iterator.hasNext()) {
            // trash the next read if it fails PF, is secondary, or is supplementary
            final SAMRecord nextRead = iterator.peek();
            if (nextRead.getReadFailsVendorQualityCheckFlag() || nextRead.isSecondaryOrSupplementary()) {
                iterator.next();
            }
            // otherwise, return it
            else {
                return justPeek ? nextRead : iterator.next();
            }
        }

        // no good reads left
        return null;
    }

}