package picard.nio;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.nio.file.*;
import java.util.HashMap;

public class GATKIOUtils {
    /**
     * Schedule a file or directory to be deleted on JVM shutdown.
     *
     * This calls {@link GATKIOUtils#deleteRecursively(Path)} on {@code fileToDelete }as a shutdown hook.
     * @param fileToDelete file or directory to be deleted recursively at JVM shutdown.
     */
    public static void deleteOnExit(final Path fileToDelete){
        DeleteRecursivelyOnExitPathHook.add(fileToDelete);
    }

    /**
     * Converts the given URI to a {@link Path} object. If the filesystem cannot be found in the usual way, then attempt
     * to load the filesystem provider using the thread context classloader. This is needed when the filesystem
     * provider is loaded using a URL classloader (e.g. in spark-submit).
     *
     * Also makes an attempt to interpret the argument as a file name if it's not a URI.
     *
     * @param uriString the URI to convert.
     * @return the resulting {@code Path}
     * @throws UserException if an I/O error occurs when creating the file system
     */
    public static Path getPath(String uriString) {
        GATKUtils.nonNull(uriString);
        URI uri;
        try {
            uri = URI.create(uriString);
        } catch (IllegalArgumentException x) {
            // not a valid URI. Caller probably just gave us a file name.
            return Paths.get(uriString);
        }
        try {
            // special case GCS, in case the filesystem provider wasn't installed properly but is available.
            if (CloudStorageFileSystem.URI_SCHEME.equals(uri.getScheme())) {
                return GATKBucketUtils.getPathOnGcs(uriString);
            }
            // Paths.get(String) assumes the default file system
            // Paths.get(URI) uses the scheme
            return uri.getScheme() == null ? Paths.get(uriString) : Paths.get(uri);
        } catch (FileSystemNotFoundException e) {
            try {
                ClassLoader cl = Thread.currentThread().getContextClassLoader();
                if ( cl == null ) {
                    throw e;
                }
                return FileSystems.newFileSystem(uri, new HashMap<>(), cl).provider().getPath(uri);
            }
            catch (ProviderNotFoundException x) {
                // TODO: this creates bogus Path on the current file system for schemes such as gendb, nonexistent, gcs
                // TODO: we depend on this code path to allow IntervalUtils to all getPath on a string that may be either
                // a literal interval or a feature file containing intervals
                // not a valid URI. Caller probably just gave us a file name or "chr1:1-2".
                return Paths.get(uriString);
            }
            catch ( IOException io ) {
                // UserException in GATK but Picard does not differentiate between e.g. PicardException vs UserException
                // Might be useful to add the UserException class in the long term
                throw new PicardException(uriString + " is not a supported path", io);
            }
        }
    }

    /**
     * Creates a temp file that will be deleted on exit
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     * @param name Prefix of the file; {@link File#createTempFile(String, String, File)} requires that this be >= 3 characters
     * @param extension Extension to concat to the end of the file.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFile(String name, String extension) {
        return createTempFileInDirectory(name, extension, null);
    }

    /**
     * Creates a temp file in a target directory that will be deleted on exit
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     * @param name Prefix of the file; {@link File#createTempFile(String, String, File)} requires that this be >= 3 characters
     * @param extension Extension to concat to the end of the file name.
     * @param targetDir Directory in which to create the temp file. If null, the default temp directory is used.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFileInDirectory(final String name, String extension, final File targetDir) {
        try {

            if (!extension.startsWith(".")) {
                extension = "." + extension;
            }

            final File file = File.createTempFile(name, extension, targetDir);
            file.deleteOnExit();

            // Mark corresponding indices for deletion on exit as well just in case an index is created for the temp file:
            new File(file.getAbsolutePath() + FileExtensions.TRIBBLE_INDEX).deleteOnExit();
            new File(file.getAbsolutePath() + FileExtensions.TABIX_INDEX).deleteOnExit();
            new File(file.getAbsolutePath() + ".bai").deleteOnExit();
            new File(file.getAbsolutePath() + ".md5").deleteOnExit();
            new File(file.getAbsolutePath().replaceAll(extension + "$", ".bai")).deleteOnExit();

            return file;
        } catch (IOException ex) {
            throw new PicardException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

    /**
     * Delete rootPath recursively
     * @param rootPath is the file/directory to be deleted
     */
    public static void deleteRecursively(final Path rootPath) {
        IOUtil.recursiveDelete(rootPath);
    }

}
