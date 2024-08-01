package picard.nio;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import com.google.cloud.storage.contrib.nio.CloudStoragePath;
import htsjdk.io.IOPath;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.utils.ValidationUtils;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.UUID;


/**
 * Derived from BucketUtils.java in GATK
 */
public class PicardBucketUtils {
    public static final String GOOGLE_CLOUD_STORAGE_FILESYSTEM_SCHEME = "gs";
    public static final String HTTP_FILESYSTEM_PROVIDER_SCHEME = "http";
    public static final String HTTPS_FILESYSTEM_PROVIDER_SCHEME = "https";
    public static final String HDFS_SCHEME = "hdfs";
    public static final String FILE_SCHEME = "file";

    // This Picard test staging bucket has a TTL of 180 days (DeleteAction with Age = 180)
    public static final String GCLOUD_PICARD_STAGING_DIRECTORY = "gs://hellbender-test-logs/staging/picard/";

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
    public static PicardHtsPath getTempFilePath(final String directory, String prefix, final String extension){
        ValidationUtils.validateArg(extension.startsWith("."), "The new extension must start with a period '.'");
        final String defaultPrefix = "tmp";

        if (directory == null){
            // File#createTempFile requires that the prefix be at least 3 characters long
            prefix = prefix.length() >= 3 ? prefix : defaultPrefix;
            return new PicardHtsPath(PicardIOUtils.createTempFile(prefix, extension));
        }

        if (isGcsUrl(directory) || isHadoopUrl(directory)){
            final PicardHtsPath path = PicardHtsPath.fromPath(randomRemotePath(directory, prefix, extension));
            PicardIOUtils.deleteOnExit(path.toPath());
            // Mark auxiliary files to be deleted
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.TRIBBLE_INDEX, true).toPath());
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.TABIX_INDEX, true).toPath());
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.BAI_INDEX, true).toPath()); // e.g. file.bam.bai
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.BAI_INDEX, false).toPath()); // e.g. file.bai
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, ".md5", true).toPath());
            return path;
        } else {
            // Assume the (non-null) directory points to a directory on a local filesystem
            prefix = prefix.length() >= 3 ? prefix : defaultPrefix;
            return new PicardHtsPath(PicardIOUtils.createTempFileInDirectory(prefix, extension, new File(directory)));
        }
    }

    /**
     * This overload of getTempFilePath takes the directory of type PicardHtsPath instead of String.
     *
     * @see #getTempFilePath(String, String, String)
     *
     */
    public static PicardHtsPath getTempFilePath(final IOPath directory, String prefix, final String extension){
        return getTempFilePath(directory.getURIString(), prefix, extension);
    }

    /**
     * Calls getTempFilePath with the empty string as the prefix.
     *
     * @see #getTempFilePath(String, String, String)
     */
    public static PicardHtsPath getTempFilePath(String directory, String extension){
        return getTempFilePath(directory, "", extension);
    }

    /**
     * Creates a temporary file in a local directory.
     *
     * @see #getTempFilePath(String, String, String)
     */
    public static PicardHtsPath getLocalTempFilePath(final String prefix, final String extension){
        return getTempFilePath((String) null, prefix, extension);
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
    public static PicardHtsPath getRandomGCSDirectory(final String relativePath){
        ValidationUtils.validateArg(relativePath.endsWith("/"), "relativePath must end in backslash '/': " + relativePath);

        return PicardHtsPath.fromPath(PicardBucketUtils.randomRemotePath(GCLOUD_PICARD_STAGING_DIRECTORY + relativePath, "", "/"));
    }

    /**
     * Picks a random name, by putting some random letters between "prefix" and "suffix".
     *
     * @param stagingLocation The folder where you want the file to be. Must start with "gs://" or "hdfs://"
     * @param prefix The beginning of the file name
     * @param suffix The end of the file name, e.g. ".tmp"
     */
    public static Path randomRemotePath(String stagingLocation, String prefix, String suffix) {
        if (isGcsUrl(stagingLocation)) {
            return getPathOnGcs(stagingLocation).resolve(prefix + UUID.randomUUID() + suffix);
        } else if (isHadoopUrl(stagingLocation)) {
            return Paths.get(stagingLocation, prefix + UUID.randomUUID() + suffix);
        } else {
            throw new IllegalArgumentException("Staging location is not remote: " + stagingLocation);
        }
    }

    /**
     * String -> Path. This *should* not be necessary (use Paths.get(URI.create(...)) instead) , but it currently is
     * on Spark because using the fat, shaded jar breaks the registration of the GCS FilesystemProvider.
     * To transform other types of string URLs into Paths, use IOUtils.getPath instead.
     */
    private static CloudStoragePath getPathOnGcs(String gcsUrl) {
        // use a split limit of -1 to preserve empty split tokens, especially trailing slashes on directory names
        final String[] split = gcsUrl.split("/", -1);
        final String BUCKET = split[2];
        final String pathWithoutBucket = String.join("/", Arrays.copyOfRange(split, 3, split.length));
        return CloudStorageFileSystem.forBucket(BUCKET).getPath(pathWithoutBucket);
    }

    /**
     *
     * @param path path to inspect
     * @return true if this path represents a gcs location
     */
    private static boolean isGcsUrl(final String path) {
        GATKUtils.nonNull(path);
        return path.startsWith(GOOGLE_CLOUD_STORAGE_FILESYSTEM_SCHEME + "://");
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
     * @param pathSpec specifier to inspect
     * @return true if this {@code GATKPath} represents a remote storage system which may benefit from prefetching (gcs or http(s))
     */
    public static boolean isEligibleForPrefetching(final IOPath pathSpec) {
        GATKUtils.nonNull(pathSpec);
        return isEligibleForPrefetching(pathSpec.getScheme());
    }

    /**
     * @param path path to inspect
     * @return true if this {@code Path} represents a remote storage system which may benefit from prefetching (gcs or http(s))
     */
    public static boolean isEligibleForPrefetching(final Path path) {
        GATKUtils.nonNull(path);
        return isEligibleForPrefetching(path.toUri().getScheme());
    }

    private static boolean isEligibleForPrefetching(final String scheme){
        return scheme != null
                && (scheme.equals(GOOGLE_CLOUD_STORAGE_FILESYSTEM_SCHEME)
                || scheme.equals(HTTP_FILESYSTEM_PROVIDER_SCHEME)
                || scheme.equals(HTTPS_FILESYSTEM_PROVIDER_SCHEME));
    }

    /**
     * Returns true if the given path is a HDFS (Hadoop filesystem) URL.
     */
    private static boolean isHadoopUrl(String path) {
        return path.startsWith(HDFS_SCHEME + "://");
    }
}
