package picard.util;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.programgroups.None;

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
@CommandLineProgramProperties(
        summary = "Provides a large, configurable, FIFO buffer that can be used to buffer input and output " +
                "streams between programs with a buffer size that is larger than that offered by native unix FIFOs (usually 64k).",
        oneLineSummary = "FIFO buffer used to buffer input and output streams with a customizable buffer size ",
        programGroup = None.class
)
public class FifoBuffer extends CommandLineProgram {
    @Argument(doc="The size of the memory buffer in bytes.")
    public int BUFFER_SIZE = 512 * 1024 * 1024;

    @Argument(doc="The size, in bytes, to read/write atomically to the input and output streams.")
    public int IO_SIZE = 64 * 1024; // 64k == most common unix pipe buffer size

    @Argument(doc="How frequently, in seconds, to report debugging statistics. Set to zero for never.")
    public int DEBUG_FREQUENCY = 0;

    @Argument(doc="Name to use for Fifo in debugging statements.", optional=true)
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
