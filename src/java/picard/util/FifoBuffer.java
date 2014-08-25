package picard.util;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * A program that is designed to act as a large memory buffer between processes that are
 * connected with unix pipes where one or more processes either produce or consume their
 * input or output in bursts.  By inserting a large memory buffer between such processes
 * each process can run at full speed and the bursts can be smoothed out.
 *
 * For example:
 *   java -jar SamToFastq.jar F=my.fastq INTERLEAVE=true | java -jar FifoBuffer | bwa mem -t 8 ...
 *
 * @author Tim Fennell
 */
public class FifoBuffer extends CommandLineProgram {
    @Usage
    public final String USAGE = "Provides a large, configurable, FIFO buffer that can be used to buffer input and output " +
            "streams between programs with a buffer size that is larger than that offered by native unix FIFOs (usually 64k).";

    @Option(doc="The size of the memory buffer in bytes.")
    public int BUFFER_SIZE = 512 * 1024 * 1024;

    @Option(doc="The size, in bytes, to read/write atomically to the input and output streams.")
    public int IO_SIZE = 64 * 1024; // 64k == most common unix pipe buffer size

    @Option(doc="How frequently, in seconds, to report debugging statistics. Set to zero for never.")
    public int DEBUG_FREQUENCY = 0;

    @Option(doc="Name to use for Fifo in debugging statements.", optional=true)
    public String NAME;

    // Standard log object
    private final Log log = Log.getInstance(FifoBuffer.class);

    /** The input stream that bytes will be read from and deposited into the byte buffer. */
    private final InputStream inputStream;

    /** The output [Print] stream which bytes read from the buffer will be emitted into. */
    private final PrintStream outputStream;

    /** Small custom exception handler for threads. */
    class LoggingExceptionHandler implements Thread.UncaughtExceptionHandler {
        public Throwable throwable;

        @Override public void uncaughtException(final Thread t, final Throwable e) {
            this.throwable = e;
            log.error(e, "Exception caught on thread ", t.getName());
        }
    }

    /**
     * Constructor that defaults to QUIET since Fifos don't do anything beyond buffering having their
     * start/end information logged is often undesirable.
     */
    public FifoBuffer(final InputStream in, final PrintStream out) {
        this.inputStream = in;
        this.outputStream = out;
        this.QUIET=true;
    }

    /**
     * Constructor that defaults to QUIET since Fifos don't do anything beyond buffering having their
     * start/end information logged is often undesirable.
     */
    public FifoBuffer() {
        this(System.in, System.out);
    }

    // Stock main method
    public static void main(final String[] args) {
        new FifoBuffer().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        final CircularByteBuffer fifo = new CircularByteBuffer(BUFFER_SIZE);

        // Input thread that reads from inputStream until it is closed and writes the contents
        // into the circular byte buffer.
        final Thread input = new Thread(new Runnable() {
            @Override public void run() {
                try {
                    final byte[] buffer = new byte[IO_SIZE];
                    int read = 0;
                    while ((read = inputStream.read(buffer)) > -1) {
                        int start = 0;
                        while (start < read) {
                            start += fifo.write(buffer, start, read-start);
                        }
                    }
                }
                catch (final IOException ioe) {
                    throw new RuntimeIOException(ioe);
                }
                finally {
                    fifo.close();
                }
            }
        });

        // Output thread that reads from the circular byte buffer until it is closed and writes
        // the results to the outputStream
        final Thread output = new Thread(new Runnable() {
            @Override public void run() {
                final byte[] buffer = new byte[IO_SIZE];
                int read = 0;
                while ((read = fifo.read(buffer, 0, buffer.length)) > 0 || !fifo.isClosed()) {
                    outputStream.write(buffer, 0, read);
                }
            }
        });

        try {
            // If debugging is turned on then start another thread that will report on the utilization of the
            // circular byte buffer every N seconds.
            if (DEBUG_FREQUENCY > 0) {
                final Thread debug = new Thread(new Runnable() {
                    @Override public void run() {
                        final NumberFormat pFmt = NumberFormat.getPercentInstance();
                        final NumberFormat iFmt = new DecimalFormat("#,##0");
                        while (true) {
                            final int capacity = fifo.getCapacity();
                            final int used     = fifo.getBytesAvailableToRead();
                            final double pct   = used / (double) capacity;
                            final String name = NAME == null ? "" : NAME + " ";
                            log.info("Fifo buffer ", name, "used ", iFmt.format(used), " / ", iFmt.format(capacity), " (", pFmt.format(pct), ").");
                            try { Thread.sleep(DEBUG_FREQUENCY * 1000); } catch (final InterruptedException ie) { /* do nothing */ }
                        }
                    }
                });

                debug.setName("BufferDebugThread");
                debug.setDaemon(true);
                debug.start();
            }

            // Start the input and output threads.
            final LoggingExceptionHandler inputExceptionHandler = new LoggingExceptionHandler();
            input.setUncaughtExceptionHandler(inputExceptionHandler);
            input.setName("Fifo Input Thread");
            input.start();

            final LoggingExceptionHandler outputExceptionHandler = new LoggingExceptionHandler();
            output.setUncaughtExceptionHandler(new LoggingExceptionHandler());
            output.setName("Fifo Output Thread");
            output.start();

            // Join on both the input and output threads to make sure that we've fully flushed the buffer
            input.join();
            output.join();

            // Double check that neither thread exited with an exception; If so propagate the exception
            if (inputExceptionHandler.throwable  != null) throw new PicardException("Exception on input thread.",  inputExceptionHandler.throwable);
            if (outputExceptionHandler.throwable != null) throw new PicardException("Exception on output thread.", outputExceptionHandler.throwable);
        }
        catch (final InterruptedException ie) {
            throw new PicardException("Interrupted!", ie);
        }

        return 0;
    }
}

/**
 * Implementation of a circular byte buffer that uses a large byte[] internally and supports basic
 * read/write operations from/to other byte[]s passed as arguments. Uses wait/nofity() to manage
 * cross-thread coordination when the buffer is either full or empty.
 */
class CircularByteBuffer {
    private final byte[] bytes;
    private final int capacity;

    private int nextWritePos = 0;
    private int bytesAvailableToWrite;
    private int nextReadPos = 0;
    private int bytesAvailableToRead = 0;

    private boolean closed = false;

    /** Constructs a buffer capable of holding the given number of bytes. */
    CircularByteBuffer(final int size) {
        this.bytes = new byte[size];
        this.capacity = this.bytes.length;
        this.bytesAvailableToWrite = this.capacity;
    }

    /**
     * Write bytes into the buffer from the supplied array.  Will attempt to read 'size' bytes beginning
     * at 'start' in the supplied array and copy them into the buffer.  If the buffer is near full or
     * cannot write 'size' bytes contiguously it may write fewer than 'size' bytes.
     *
     * @return the number of bytes read from the input array and copied into the buffer
     */
    synchronized int write(final byte[] bytes, final int start, final int size) {
        if (closed) throw new IllegalStateException("Cannot write to closed buffer.");

        try { if (this.bytesAvailableToWrite == 0) wait(); }
        catch (final InterruptedException ie) {throw new PicardException("Interrupted while waiting to write to fifo.", ie); }

        final int writePos      = this.nextWritePos;
        final int distanceToEnd = this.capacity - writePos;
        final int available     = distanceToEnd < this.bytesAvailableToWrite ? distanceToEnd : this.bytesAvailableToWrite;
        final int length        = available < size ? available : size;

        System.arraycopy(bytes, start, this.bytes, writePos, length);
        this.bytesAvailableToWrite -= length;
        this.bytesAvailableToRead  += length;
        this.nextWritePos = (writePos + length) % this.capacity;
        notify();
        return length;
    }

    /**
     * Read bytes from the buffer into the supplied array.  Will attempt to read 'size' bytes and write
     * them into the supplied array beginning at index 'start' in the supplied array.  If the buffer is
     * near empty or cannot read 'size' bytes contiguously it may write fewer than 'size' bytes.
     *
     * @return the number of bytes read from the buffer and copied into the input array
     */
    synchronized int read(final byte[] bytes, final int start, final int size) {
        try { if (this.bytesAvailableToRead == 0 && !closed) wait(); }
        catch (final InterruptedException ie) {throw new PicardException("Interrupted while waiting to read from fifo.", ie); }

        final int readPos       = this.nextReadPos;
        final int distanceToEnd = this.capacity - readPos;
        final int available     = distanceToEnd < this.bytesAvailableToRead ? distanceToEnd : this.bytesAvailableToRead;
        final int length        = available < size ? available : size;

        System.arraycopy(this.bytes, readPos, bytes, start, length);
        this.bytesAvailableToRead  -= length;
        this.bytesAvailableToWrite += length;
        this.nextReadPos = (readPos + length) % this.capacity;
        notify();
        return length;
    }

    /** Signals that the buffer is closed and no further writes will occur. */
    synchronized void close() {
        this.closed = true;
        notify();
    }

    /** Returns true if the buffer is closed, false otherwise. */
    synchronized boolean isClosed() {
        return this.closed;
    }

    /** Returns the total capacity of the buffer (empty+filled). */
    int getCapacity() { return this.capacity; }

    /** Returns the number of bytes that are in the buffer at the time of the method invocation. */
    synchronized int getBytesAvailableToRead() { return this.bytesAvailableToRead; }
}