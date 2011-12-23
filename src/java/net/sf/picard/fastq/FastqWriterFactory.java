package net.sf.picard.fastq;

import net.sf.samtools.Defaults;

import java.io.File;

/**
 * Factory class for creating FastqWriter objects.
 *
 * @author Tim Fennell
 */
public class FastqWriterFactory {
    boolean useAsyncIo = Defaults.USE_ASYNC_IO;

    /** Sets whether or not to use async io (i.e. a dedicated thread per writer. */
    public void setUseAsyncIo(final boolean useAsyncIo) { this.useAsyncIo = useAsyncIo; }

    public FastqWriter newWriter(final File out) {
        final FastqWriter writer = new BasicFastqWriter(out);
        if (useAsyncIo) {
            return new AsyncFastqWriter(writer, AsyncFastqWriter.DEFAULT_QUEUE_SIZE);
        }
        else {
            return writer;
        }
    }
}
