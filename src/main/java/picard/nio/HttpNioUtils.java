package picard.nio;

import org.broadinstitute.http.nio.HttpFileSystemProvider;
import org.broadinstitute.http.nio.HttpFileSystemProviderSettings;
import org.broadinstitute.http.nio.RetryHandler;

import java.time.Duration;

/**
 * This class provides a way to easily configure the HttpNioProvider
 *
 * This class contains hard-coded setting that have been found to work for the access patterns that characterize genomics
 * work.
 *
 */
public final class HttpNioUtils {

    private HttpNioUtils() {}
    public static final Duration MAX_TIMEOUT = Duration.ofMillis(120_000);
    public static final int MAX_RETRIES = 20;

    public static void initialize() {
        final HttpFileSystemProviderSettings.RetrySettings retrySettings = new HttpFileSystemProviderSettings.RetrySettings(
                MAX_RETRIES,
                RetryHandler.DEFAULT_RETRYABLE_HTTP_CODES,
                RetryHandler.DEFAULT_RETRYABLE_EXCEPTIONS,
                RetryHandler.DEFALT_RETRYABLE_MESSAGES,
                e -> false);

        final HttpFileSystemProviderSettings settings = new HttpFileSystemProviderSettings(
                MAX_TIMEOUT,
                HttpFileSystemProviderSettings.DEFAULT_SETTINGS.redirect(),
                retrySettings);

        HttpFileSystemProvider.setSettings(settings);
    }
}
