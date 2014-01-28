package org.broadinstitute.variant.variantcontext.writer;

import net.sf.samtools.util.AbstractAsyncWriter;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFHeader;

/**
 * AsyncVariantContextWriter that can be wrapped around an underlying AsyncVariantContextWriter to provide asynchronous output. Records
 * added are placed into a queue, the queue is then drained into the underlying VariantContextWriter by a thread owned
 * by the instance.
 *
 * Exceptions experienced by the writer thread will be emitted back to the caller in subsequent calls to either
 * add() or close().
 *
 * @author George Grant
 */
public class AsyncVariantContextWriter extends AbstractAsyncWriter<VariantContext> implements VariantContextWriter {
    private final VariantContextWriter underlyingWriter;

    /**
     * Creates a new AsyncVariantContextWriter wrapping the provided VariantContextWriter.
     */
    public AsyncVariantContextWriter(final VariantContextWriter out) {
        this(out, DEFAULT_QUEUE_SIZE);
    }

    /**
     * Creates an AsyncVariantContextWriter wrapping the provided VariantContextWriter and using the specified
     * queue size for buffer VariantContexts.
     */
    public AsyncVariantContextWriter(final VariantContextWriter out, final int queueSize) {
        super(queueSize);
        this.underlyingWriter = out;
    }

    @Override protected void synchronouslyWrite(final VariantContext item) { this.underlyingWriter.add(item); }

    @Override protected void synchronouslyClose() { this.underlyingWriter.close();  }

    @Override protected final String getThreadNamePrefix() { return "VariantContextWriterThread-"; }

    public void add(final VariantContext vc) {
        write(vc);
    }

    public void writeHeader(final VCFHeader header) {
        this.underlyingWriter.writeHeader(header);
    }
}
