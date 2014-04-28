package net.sf.samtools;


import net.sf.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;

/**
 * OutputStream implementation that writes output to an underlying output stream while also copying the
 * same bytes to a PipedOutputStream that routes the data back into an Indexer to generate a BAMIndex
 * by inflating and decoding the stream and feeding the SAMRecords to a BAMIndexer.
 */
class StreamInflatingIndexingOutputStream extends OutputStream {
    private final OutputStream s1;
    private final PipedOutputStream s2;
    private final Thread thread;

    public StreamInflatingIndexingOutputStream(final OutputStream s1, final File indexFile) {
        try {
            this.s1 = s1;
            this.s2 = new PipedOutputStream();
            final PipedInputStream pin = new PipedInputStream(this.s2, Defaults.NON_ZERO_BUFFER_SIZE);
            this.thread = new Thread(new Indexer(indexFile, pin), "BamIndexingThread");
            this.thread.start();
        }
        catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    @Override public final void write(final int b) throws IOException { this.s1.write(b); this.s2.write(b); }
    @Override public final void write(final byte[] b) throws IOException { this.s1.write(b); this.s2.write(b); }
    @Override public final void write(final byte[] b, final int off, final int len) throws IOException { this.s1.write(b, off, len); this.s2.write(b, off, len); }
    @Override public final void flush() throws IOException { this.s1.flush(); this.s2.flush(); }

    @Override public final void close() throws IOException {
        this.s1.close();
        this.s2.close();

        try { this.thread.join(); }
        catch (final InterruptedException ie) { throw new RuntimeException(ie); }
    }
}

/**
 * A little class that takes an InputStream from which it reads a BAM file, generates
 * a BAMIndex and then writes the index to the File provided.  All operations are designed
 * to be carried out in a separate thread.
 */
class Indexer implements Runnable {
    private final File index;
    private final InputStream stream;

    /** Constructs an indexer that reads from the stream provided and writes an index to the File provided. */
    Indexer(final File index, final InputStream stream) {
        this.index = index;
        this.stream = stream;
    }

    /** Runnable implementation that reads the entire stream and writes the index. */
    @Override
    public void run() {
        final SAMFileReader in = new SAMFileReader(this.stream);
        in.enableFileSource(true);
        in.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        in.enableCrcChecking(false);
        final BAMIndexer indexer = new BAMIndexer(this.index, in.getFileHeader());
        for (final SAMRecord rec : in) {
            indexer.processAlignment(rec);
        }

        indexer.finish();
        in.close();
    }
}
