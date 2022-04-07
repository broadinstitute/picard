package picard.util;

import htsjdk.samtools.util.Log;

import java.time.Duration;
import java.util.concurrent.TimeUnit;

public class ThreadPoolExecutorUtil {
    private static final Log log = Log.getInstance(ThreadPoolExecutorUtil.class);

    public static void awaitThreadPoolTermination(final String executorName, final ThreadPoolExecutorWithExceptions executorService,
                                                  final Duration timeBetweenChecks) {
        try {
            while (!executorService.awaitTermination(timeBetweenChecks.getSeconds(), TimeUnit.SECONDS)) {
                if (executorService.hasError()) {
                    log.error(executorService.exception,
                            String.format("%s terminating because a worker thread had an exception.", executorName));
                    executorService.shutdownNow();
                    break;
                }
                log.info(String.format("%s waiting for job completion. Finished jobs - %d : Running jobs - %d : Queued jobs  - %d",
                        executorName, executorService.getCompletedTaskCount(), executorService.getActiveCount(),
                        executorService.getQueue().size()));
            }
        } catch (final InterruptedException e) {
            log.error("Interrupted exception caught: ", e);
        }
    }
}
