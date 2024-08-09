package picard.nio;

import htsjdk.io.IOPath;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.utils.ValidationUtils;
import picard.PicardException;

import java.util.UUID;

/**
 * Derived from BucketUtils.java in GATK
 */
public class PicardBucketUtils {
    public static final String GOOGLE_CLOUD_STORAGE_FILESYSTEM_SCHEME = "gs";
    public static final String HTTP_FILESYSTEM_PROVIDER_SCHEME = "http";
    public static final String HTTPS_FILESYSTEM_PROVIDER_SCHEME = "https";
    public static final String FILE_SCHEME = "file";

    // This Picard test staging bucket has a TTL of 180 days (DeleteAction with Age = 180)
    public static final String GCLOUD_PICARD_STAGING_DIRECTORY_STR = "gs://hellbender-test-logs/staging/picard/";
    public static final PicardHtsPath GCLOUD_PICARD_STAGING_DIRECTORY = new PicardHtsPath(GCLOUD_PICARD_STAGING_DIRECTORY_STR);


    // slashes omitted since hdfs paths seem to only have 1 slash which would be weirder to include than no slashes
    private PicardBucketUtils(){} //private so that no one will instantiate this class

    /**
     * Get a temporary file path based on the prefix and extension provided.
     * This file (and possible indexes associated with it) will be scheduled for deletion on shutdown.
     *
     * @param directory the directory where the temporary fill will be placed. May be null.
     *               For remote paths this should be a valid URI to root the temporary file in (e.g. gs://hellbender/staging/)
     *               If the directory is not a GCS or hadoop URL (or if it is null), the resulting temp file will be placed in the local tmp folder.
     * @param prefix the prefix to prepend before the randomly generated characters. Must have length >= 3 for local temp files.
     * @param extension an extension for the temporary file path. Must start with a period e.g. ".txt"
     * @return a new temporary path of the form [directory]/[prefix][random chars][.extension]
     *
     */
    public static IOPath getTempFilePath(final IOPath directory, String prefix, final String extension){
        ValidationUtils.validateArg(extension.startsWith("."), "The new extension must start with a period '.'");
        final String defaultPrefix = "tmp";

        if (directory == null){
            // If directory = null, we are creating a local temp file.
            // File#createTempFile requires that the prefix be at least 3 characters long
            prefix = prefix.length() >= 3 ? prefix : defaultPrefix;
            return new PicardHtsPath(PicardIOUtils.createTempFile(prefix, extension));
        } else if (PicardBucketUtils.isLocalPath(directory)) {
            // Assume the (non-null) directory points to a directory on a local filesystem
            prefix = prefix.length() >= 3 ? prefix : defaultPrefix;
            return new PicardHtsPath(PicardIOUtils.createTempFileInDirectory(prefix, extension, directory.toPath().toFile()));
        } else {
            if (isSupportedCloudFilesystem(directory)) {
                final IOPath path = randomRemotePath(directory, prefix, extension);
                PicardIOUtils.deleteOnExit(path.toPath());
                // Mark auxiliary files to be deleted
                PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.TRIBBLE_INDEX, true).toPath());
                PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.TABIX_INDEX, true).toPath());
                PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.BAI_INDEX, true).toPath()); // e.g. file.bam.bai
                PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.BAI_INDEX, false).toPath()); // e.g. file.bai
                PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, ".md5", true).toPath());
                return path;
            } else {
                throw new PicardException("Unsupported cloud filesystem: " + directory.getURIString());
            }
        }
    }

    /**
     * Creates a temporary file in a local directory.
     *
     * @see #getTempFilePath(IOPath, String, String)
     */
    public static IOPath getLocalTempFilePath(final String prefix, final String extension){
        return getTempFilePath(null, prefix, extension);
    }

    /**
     * Creates a PicardHtsPath object to a "directory" on a Google Cloud System filesystem with a randomly generated URI.
     *
     * Note that the notion of directories does not exist in GCS. Thus, by "directory,"
     * we mean a path object with a randomly generated URI ending in "/", which
     * the caller can use as a root URI/path for other files to be created e.g. via PicardHtsPath::resolve.
     *
     * Note that this method does *not* create an actual directory/file on GCS that one can write to, delete, or otherwise manipulate.
     *
     * See: https://stackoverflow.com/questions/51892343/google-gsutil-create-folder
     *
     * @param relativePath The relative location for the new "directory" under the harcoded staging bucket with a TTL set e.g. "test/RevertSam/".
     * @return A PicardHtsPath object to a randomly generated "directory" e.g. "gs://hellbender-test-logs/staging/picard/test/RevertSam/{randomly-generated-string}/"
     */
    public static IOPath getRandomGCSDirectory(final String relativePath){
        ValidationUtils.validateArg(relativePath.endsWith("/"), "relativePath must end in backslash '/': " + relativePath);

        return PicardBucketUtils.randomRemotePath(PicardHtsPath.resolve(GCLOUD_PICARD_STAGING_DIRECTORY, relativePath), "", "/");
    }

    /**
     * Picks a random name, by putting some random letters between "prefix" and "suffix".
     *
     * @param stagingLocation The folder where you want the file to be. Must start with "gs://" or "hdfs://"
     * @param prefix The beginning of the file name
     * @param suffix The end of the file name, e.g. ".tmp"
     */
    public static IOPath randomRemotePath(final IOPath stagingLocation, final String prefix, final String suffix) {
        if (isGcsUrl(stagingLocation)) {
            return PicardHtsPath.resolve(stagingLocation, prefix + UUID.randomUUID() + suffix);
        } else {
            throw new IllegalArgumentException("Staging location is not remote: " + stagingLocation);
        }
    }

    /**
     *
     * Return true if this {@code PicardHTSPath} represents a gcs URI.
     * @param pathSpec specifier to inspect
     * @return true if this {@code PicardHTSPath} represents a gcs URI.
     */
    public static boolean isGcsUrl(final IOPath pathSpec) {
        GATKUtils.nonNull(pathSpec);
        return pathSpec.getScheme().equals(GOOGLE_CLOUD_STORAGE_FILESYSTEM_SCHEME);
    }

    /**
     * @param path specifier to inspect
     * @return true if this {@code IOPath} represents a remote storage system which may benefit from prefetching (gcs or http(s))
     */
    public static boolean isEligibleForPrefetching(final IOPath path) {
        GATKUtils.nonNull(path);
        final String scheme = path.getScheme();
        return scheme != null
                && (scheme.equals(GOOGLE_CLOUD_STORAGE_FILESYSTEM_SCHEME)
                || scheme.equals(HTTP_FILESYSTEM_PROVIDER_SCHEME)
                || scheme.equals(HTTPS_FILESYSTEM_PROVIDER_SCHEME));
    }

    public static boolean isLocalPath(final IOPath path){
        return path.getScheme().equals(FILE_SCHEME);
    }

    /**
     * As of August 2024, we only support Google Cloud.
     * Will add other filesystems (e.g. Azure, AWS) when ready.
     *
     * @return whether the cloud filesystem is currently supported by Picard.
     */
    public static boolean isSupportedCloudFilesystem(final IOPath path){
        ValidationUtils.validateArg(! isLocalPath(path), "isSupportedCloudFilesystem should be called on a cloud path but was given: " +
                path.getURIString());
        return isGcsUrl(path);
    }
}
