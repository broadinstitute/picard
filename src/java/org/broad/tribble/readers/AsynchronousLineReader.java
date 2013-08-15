package org.broad.tribble.readers;

import net.sf.samtools.util.CloserUtil;
import org.broad.tribble.TribbleException;

import java.io.Reader;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

/**
 * A LineReader implementation that delegates the work of reading and fetching lines to another thread.  The thread terminates when it
 * encounters EOF in the underlying reader, or when this LineReader is closed.
 *
 * @author mccowan
 */
public class AsynchronousLineReader implements LineReader {
    public static final int DEFAULT_NUMBER_LINES_BUFFER = 100;
    
    private final LongLineBufferedReader bufferedReader;
    private final BlockingQueue<String> lineQueue;
    private final Thread worker;
    private volatile Exception workerException = null;
    private volatile boolean eofReached = false;

    public AsynchronousLineReader(final Reader reader, final int lineReadAheadSize) {
        bufferedReader = new LongLineBufferedReader(reader);
        lineQueue = new LinkedBlockingQueue<String>(lineReadAheadSize);
        worker = new Thread(new Worker());
        worker.start();
    }

    public AsynchronousLineReader(final Reader reader) {
        this(reader, DEFAULT_NUMBER_LINES_BUFFER);
    }

    @Override
    public String readLine() {
        try {
            // Continually poll until we get a result, unless the underlying reader is finished.
            for (; ; ) {
                checkAndThrowIfWorkerException();
                final String pollResult = this.lineQueue.poll(100, TimeUnit.MILLISECONDS); // Not ideal for small files.
                if (pollResult == null) {
                    if (eofReached) {
                        checkAndThrowIfWorkerException();
                        return lineQueue.poll(); // If there is nothing left, returns null as expected.  Otherwise, grabs next element.
                    }
                } else {
                    return pollResult;
                }
            }
        } catch (InterruptedException e) {
            throw new TribbleException("Line polling interrupted.", e);
        }
    }

    private void checkAndThrowIfWorkerException() {
        if (workerException != null) {
            throw new TribbleException("Exception encountered in worker thread.", workerException);
        }
    }

    @Override
    public void close() {
        this.worker.interrupt(); // Allow the worker to close gracefully.
    } 

    private class Worker implements Runnable {
        @Override
        public void run() {
            try {
                for (; ; ) {
                    final String line = bufferedReader.readLine();
                    if (line == null) {
                        eofReached = true;
                        break;
                    } else {
                        try {
                            lineQueue.put(line);
                        } catch (InterruptedException e) {
                            /**
                             * A thread interruption is not an exceptional state: it means a {@link AsynchronousLineReader#close();} has 
                             * been called, so shut down gracefully.
                             */
                            break;
                        }
                    }
                }
            } catch (Exception e) {
                AsynchronousLineReader.this.workerException = e;
            } finally {
                CloserUtil.close(AsynchronousLineReader.this.bufferedReader);
            }
        }
    }
}
