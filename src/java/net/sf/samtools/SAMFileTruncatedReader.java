package net.sf.samtools;

import java.io.File;
import java.util.NoSuchElementException;

/**
 * A truncated form of a SAMFileReader that iterates over a limited number of records.
 *
 * @author mccowan@broadinstitute.org
 */
public class SAMFileTruncatedReader extends SAMFileReader {
    private class TruncatedIterator implements SAMRecordIterator {
        final SAMRecordIterator i;
        final long max;
        long currentRecord = 0;

        TruncatedIterator(final SAMRecordIterator i, final long max) {
            this.i = i;
            this.max = max;
        }

        public boolean hasNext() {
            return i.hasNext() && max != currentRecord;
        }

        public SAMRecord next() {
            if (this.hasNext()) {
                currentRecord += 1;
                return i.next();
            } else {
                throw new NoSuchElementException();
            }
        }

        public void remove() {
            i.remove();
        }

        public void close() {
            i.close();
        }

        public SAMRecordIterator assertSorted(final SAMFileHeader.SortOrder sortOrder) {
            return i.assertSorted(sortOrder);
        }
    }

    private final long maxRecordsToIterate;
    
    /**
     * @param input The SAM file
     * @param max   The maximum number of records to read from the file via iterator() methods
     */
    public SAMFileTruncatedReader(final File input, final long max) {
        super(input);
        this.maxRecordsToIterate = max;
    }
    
    @Override
    public SAMRecordIterator iterator() {
        return new TruncatedIterator(super.iterator(), maxRecordsToIterate);
    }

    @Override
    public SAMRecordIterator iterator(final SAMFileSpan chunks) {
        return new TruncatedIterator(super.iterator(chunks), maxRecordsToIterate);
    }
}
