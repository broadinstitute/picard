package net.sf.picard.fastq;

import net.sf.samtools.util.AbstractAsyncWriter;

/**
 * Implementation of a FastqWriter that provides asynchronous output.
 * @author Tim Fennell
 */
public class AsyncFastqWriter extends AbstractAsyncWriter<FastqRecord> implements FastqWriter {
    private final FastqWriter writer;

    public AsyncFastqWriter(final FastqWriter out, final int queueSize) {
        super(queueSize);
        this.writer = out;
    }

    @Override protected String getThreadNamePrefix() { return "FastqWriterThread-"; }
    @Override protected void synchronouslyWrite(final FastqRecord item) { this.writer.write(item); }
    @Override protected void synchronouslyClose() { this.writer.close(); }
}
