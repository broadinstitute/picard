package picard.util;

import picard.PicardException;

import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * This version of the thread pool executor will throw an exception if any of the internal jobs have throw exceptions
 * while executing
 */
public class ThreadPoolExecutorWithExceptions extends ThreadPoolExecutor {
    public Throwable exception = null;

    /**
     * Creates a fixed size thread pool executor that will rethrow exceptions from submitted jobs.
     *
     * @param threads The number of threads in the executor pool.
     */
    public ThreadPoolExecutorWithExceptions(final int threads) {
        super(threads, threads, 0, TimeUnit.SECONDS, new LinkedBlockingDeque<>());
    }

    @Override
    protected void afterExecute(final Runnable r, Throwable t) {
        if (t == null && r instanceof Future<?>) {
            try {
                final Future<?> future = (Future<?>) r;
                if (future.isDone()) {
                    future.get();
                }
            } catch (final CancellationException ce) {
                t = ce;
            } catch (final ExecutionException ee) {
                t = ee.getCause();
            } catch (final InterruptedException ie) {
                Thread.currentThread().interrupt(); // ignore/reset
            }
        }
        if (t != null) {
            exception = t;
            throw new PicardException(t.getMessage(), t);
        }
    }

    @Override
    protected void beforeExecute(Thread t, Runnable r) {
        super.beforeExecute(t, r);
        t.setUncaughtExceptionHandler((t1, e) -> {
            throw new PicardException("Uncaught exception in thread: " + t1.getName() + " : " + e.getMessage(), e);
        });
    }
}
