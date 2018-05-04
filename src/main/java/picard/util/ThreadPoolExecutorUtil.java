package picard.util;

import htsjdk.samtools.util.Log;

import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

public class ThreadPoolExecutorUtil {
    private static final Log log = Log.getInstance(ThreadPoolExecutorUtil.class);

    public static void awaitThreadPoolTermination(final String executorName, final ThreadPoolExecutor executorService) {
        try {
            while (!executorService.awaitTermination(300, TimeUnit.SECONDS)) {
                log.info(String.format("%s waiting for job completion. Finished jobs - %d : Running jobs - %d : Queued jobs  - %d",
                        executorName, executorService.getCompletedTaskCount(), executorService.getActiveCount(),
                        executorService.getQueue().size()));
            }
        } catch (final InterruptedException e) {
            log.error("Interrupted exception caught: ", e);
        }
    }
}
