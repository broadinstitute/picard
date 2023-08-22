package picard.nio;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import com.google.cloud.storage.contrib.nio.CloudStoragePath;
import htsjdk.samtools.util.FileExtensions;
import org.apache.commons.io.FilenameUtils;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.UUID;


/**
 * Copied from BucketUtils.java in GATK
 * To be replaced once the original GATK BucketUtils.java is ported to htsjdk
 */
public class PicardBucketUtils {
    // In GATK these are accessed as e.g. GoogleCloudStorageFileSystem.SCHEME
    public static final String GOOGLE_CLOUD_STORAGE_FILESYSTEM_SCHEME = "gs";
    public static final String HTTP_FILESYSTEM_PROVIDER_SCHEME = "http";
    public static final String HTTPS_FILESYSTEM_PROVIDER_SCHEME = "https";
    public static final String HDFS_SCHEME = "hdfs";


    // chris: should be able to get rid of these '://' Instead, compare getScheme()
    private static final String GCS_PREFIX = "gs://";
    private static final String HTTP_PREFIX = "http://";
    private static final String HTTPS_PREFIX = "https://";
    private static final String HDFS_PREFIX = HDFS_SCHEME + "://";

    // slashes omitted since hdfs paths seem to only have 1 slash which would be weirder to include than no slashes
    public static final String FILE_PREFIX = "file:";

    private PicardBucketUtils(){} //private so that no one will instantiate this class

    /**
     * Get a temporary file path based on the prefix and extension provided.
     * This file (and possible indexes associated with it) will be scheduled for deletion on shutdown
     *
     * @param prefix a prefix for the file name
     *               for remote paths this should be a valid URI to root the temporary file in (e.g. gs://hellbender/staging/)
     *               there is no guarantee that this will be used as the root of the tmp file name, a local prefix may be placed in the tmp folder for example
     * @param extension and extension for the temporary file path, the resulting path will end in this. It should include the period
     *                  e.g. ".txt" extension of "txt" results in a filename "prefixtxt" without a period.
     * @return a path to use as a temporary file, on remote file systems which don't support an atomic tmp file reservation a path is chosen with a long randomized name
     *
     */
    public static PicardHtsPath getTempFilePath(String prefix, String extension){
        if (isGcsUrl(prefix) || (isHadoopUrl(prefix))){
            // tsato: should we convert to PicardHtsPath or stay in Path...
            final PicardHtsPath path = PicardHtsPath.fromPath(randomRemotePath(prefix, "", extension));
            PicardIOUtils.deleteOnExit(path.toPath());
            // Mark auxiliary files to be deleted
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.TRIBBLE_INDEX, true).toPath()); // tsato: check append or not
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.TABIX_INDEX, true).toPath());
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.BAI_INDEX, true).toPath()); // e.g. file.bam.bai
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, ".md5", true).toPath());
            PicardIOUtils.deleteOnExit(PicardHtsPath.replaceExtension(path, FileExtensions.BAI_INDEX, false).toPath()); // e.g. file.bai
            return path;
        } else {
            return new PicardHtsPath(PicardIOUtils.createTempFile(prefix, extension));
        }
    }



    /**
     * Picks a random name, by putting some random letters between "prefix" and "suffix".
     *
     * @param stagingLocation The folder where you want the file to be. Must start with "gs://" or "hdfs://"
     * @param prefix The beginning of the file name
     * @param suffix The end of the file name, e.g. ".tmp"
     */
    private static Path randomRemotePath(String stagingLocation, String prefix, String suffix) {
        if (isGcsUrl(stagingLocation)) {
            // Go through URI because Path.toString isn't guaranteed to include the "gs://" prefix.
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
    public static CloudStoragePath getPathOnGcs(String gcsUrl) { // tsato: return cloud storage path?
        // use a split limit of -1 to preserve empty split tokens, especially trailing slashes on directory names
        final String[] split = gcsUrl.split("/", -1);
        final String BUCKET = split[2];
        final String pathWithoutBucket = String.join("/", Arrays.copyOfRange(split, 3, split.length));
        return CloudStorageFileSystem.forBucket(BUCKET).getPath(pathWithoutBucket);
    }

    /**
     * @param path path to inspect
     * @return true if this path represents a gcs location
     */
    // TODO: to be deleted. Or, since using it as part of makeTempFile is defensible; just make private for now
    private static boolean isGcsUrl(final String path) {
        GATKUtils.nonNull(path);
        return path.startsWith(GCS_PREFIX); // tsato: red flag
    } // tsato: should be a method in PicardHtsPath

    /**
     *
     * Return true if this {@code PicardHTSPath} represents a gcs URI.
     * @param pathSpec specifier to inspect
     * @return true if this {@code PicardHTSPath} represents a gcs URI.
     */
    public static boolean isGcsUrl(final PicardHtsPath pathSpec) {
        GATKUtils.nonNull(pathSpec);
        return pathSpec.getScheme().equals(GOOGLE_CLOUD_STORAGE_FILESYSTEM_SCHEME);
    }

    /**
     * @param pathSpec specifier to inspect
     * @return true if this {@code GATKPath} represents a remote storage system which may benefit from prefetching (gcs or http(s))
     */
    public static boolean isEligibleForPrefetching(final PicardHtsPath pathSpec) {
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
        return path.startsWith(HDFS_PREFIX);
    }
}
