/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.nio;

import com.google.cloud.http.HttpTransportOptions;
import com.google.cloud.storage.StorageOptions;
import com.google.cloud.storage.contrib.nio.CloudStorageConfiguration;
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider;
import shaded.cloud_nio.com.google.api.gax.retrying.RetrySettings;
import shaded.cloud_nio.org.threeten.bp.Duration;


/**
 * This class serves as a connection to google's implementation of nio support for GCS housed files.
 *
 * While the actual code required to connect isn't packaged with Picard (only compiled), the Readme.md file in the
 * github repository describes how it can be used. Additionally, Picard is bundled in GATK4, and its tools exposed via
 * the GATK engine, since the nio library _is_ included the GATK4 jar. NIO enabled tools in picard can connect to
 * GCS when called through GATK4.
 *
 * This class contains hard-coded setting that have been found to work for the access patterns that characterize genomics
 * work. In the future it would make sense to expose these parameters so that they can be controlled via the commandline.
 * However, as the list of Path-enabled tools in picard is small, there seems to be little impetus to do so right now.
 *
 *
 */
class GoogleStorageUtils {

    public static void initialize() {
        CloudStorageFileSystemProvider.setDefaultCloudStorageConfiguration(GoogleStorageUtils.getCloudStorageConfiguration(20));
        CloudStorageFileSystemProvider.setStorageOptions(GoogleStorageUtils.setGenerousTimeouts(StorageOptions.newBuilder()).build());
    }

    /** The config we want to use. **/
    private static CloudStorageConfiguration getCloudStorageConfiguration(int maxReopens) {
        return CloudStorageConfiguration.builder()
                // if the channel errors out, re-open up to this many times
                .maxChannelReopens(maxReopens)
                .build();
    }

    private static StorageOptions.Builder setGenerousTimeouts(StorageOptions.Builder builder) {
        return builder
                .setTransportOptions(HttpTransportOptions.newBuilder()
                        .setConnectTimeout(120_000)
                        .setReadTimeout(120_000)
                        .build())
                .setRetrySettings(RetrySettings.newBuilder()
                        .setMaxAttempts(15)
                        .setMaxRetryDelay(Duration.ofMillis(256_000L))
                        .setTotalTimeout(Duration.ofMillis(4000_000L))
                        .setInitialRetryDelay(Duration.ofMillis(1000L))
                        .setRetryDelayMultiplier(2.0)
                        .setInitialRpcTimeout(Duration.ofMillis(180_000L))
                        .setRpcTimeoutMultiplier(1.0)
                        .setMaxRpcTimeout(Duration.ofMillis(180_000L))
                        .build());
    }
}
