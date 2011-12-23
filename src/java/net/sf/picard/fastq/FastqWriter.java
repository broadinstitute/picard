package net.sf.picard.fastq;

import java.io.Closeable;

/**
 * Simple interface for a class that can write out fastq records.
 *
 * @author Tim Fennell
 */
public interface FastqWriter {
    void write(final FastqRecord rec);
    void close();
}
