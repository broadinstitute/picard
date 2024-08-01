package picard.util;

import picard.nio.PicardHtsPath;

public final class GCloudTestUtils {
    /**
     *  This is a public requester pays bucket owned by the broad-gatk-test project.
     *  It must be owned by a different project than the service account doing the testing or the test may fail because it can access the
     *  file directly through alternative permissions.
     */
    public static final PicardHtsPath REQUESTER_PAYS_BUCKET_DEFAULT = new PicardHtsPath("gs://hellbender-requester-pays-test/");

    public static final String TEST_INPUTS_DEFAULT_GCLOUD_STR = "gs://hellbender/test/resources/";
    public static final String TEST_STAGING_DEFAULT_GCLOUD_STR = "gs://hellbender-test-logs/staging/";
    public static final String TEST_GCLOUD_PROJECT_DEFAULT_STR = "broad-dsde-dev";
    public static final PicardHtsPath TEST_INPUTS_DEFAULT_GCLOUD = new PicardHtsPath(TEST_INPUTS_DEFAULT_GCLOUD_STR);
    public static final PicardHtsPath TEST_STAGING_DEFAULT_GCLOUD = new PicardHtsPath(TEST_STAGING_DEFAULT_GCLOUD_STR);
    public static final PicardHtsPath TEST_OUTPUT_DEFAULT_GCLOUD = PicardHtsPath.resolve(TEST_STAGING_DEFAULT_GCLOUD, "picard/");

    /**
     * A publicly readable GCS bucket set as requester pays, this should not be owned by the same project that is set
     * as {@link #getTestProject()} or the tests for requester pays access may be invalid.
     *
     * @return PICARD_REQUESTER_PAYS_BUCKET env. var if defined, or {@value GCloudTestUtils#REQUESTER_PAYS_BUCKET_DEFAULT.getURIString()},
     * wrapped in PicardHtsPath.
     */
    public static PicardHtsPath getRequesterPaysBucket() {
        return new PicardHtsPath(getSystemProperty("PICARD_REQUESTER_PAYS_BUCKET", REQUESTER_PAYS_BUCKET_DEFAULT.getURIString()));
    }

    private static String getSystemProperty(final String variableName, final String defaultValue) {
        final String valueFromEnvironment = System.getProperty(variableName);
        return valueFromEnvironment == null || valueFromEnvironment.isEmpty()? defaultValue : valueFromEnvironment;
    }

    /**
     * name of the google cloud project that stores the data and will run the code
     *
     * @return PICARD_TEST_PROJECT env. var if defined or {@value TEST_GCLOUD_PROJECT_DEFAULT_STR}
     */
    public static String getTestProject() {
        return getSystemProperty("PICARD_TEST_PROJECT", TEST_GCLOUD_PROJECT_DEFAULT_STR);
    }

    /**
     * A writable GCS path where java files can be cached and temporary test files can be written,
     * of the form gs://bucket/, or gs://bucket/path/.
     *
     * @return PICARD_TEST_STAGING env. var if defined, or {@value #TEST_STAGING_DEFAULT.getURIString()}
     */
    public static String getTestStaging() {
        return getSystemProperty("PICARD_TEST_STAGING", TEST_STAGING_DEFAULT_GCLOUD.getURIString());
    }

    /**
     * A GCS path where the test inputs are stored.
     * The value of PICARD_TEST_INPUTS should end in a "/" (for example, "gs://hellbender/test/resources/")
     *
     * @return PICARD_TEST_INPUTS env. var if defined or {@value #TEST_INPUTS_DEFAULT.getURIString()}.
     */
    public static String getTestInputPath() {
        return getSystemProperty("PICARD_TEST_INPUTS", TEST_INPUTS_DEFAULT_GCLOUD.getURIString());
    }

}
