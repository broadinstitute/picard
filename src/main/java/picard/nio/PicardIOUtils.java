package picard.nio;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

/**
 * Ported from GATKIOUtils.java
 */
public class PicardIOUtils {
    /**
     * Schedule a file or directory to be deleted on JVM shutdown.
     *
     * This calls {@link PicardIOUtils#deleteRecursively(Path)} on {@code fileToDelete }as a shutdown hook.
     * @param fileToDelete file or directory to be deleted recursively at JVM shutdown.
     */
    public static void deleteOnExit(final Path fileToDelete){
        DeleteRecursivelyOnExitPathHook.add(fileToDelete);
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
     * @param extension Extension to concat to the end of the file name e.g. ".txt"
     * @param targetDir Directory in which to create the temp file. If null, the default temp directory is used.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     *
     * TODO: consolidate this and BucketUtils::getTempFilePath
     */
    public static File createTempFileInDirectory(final String name, String extension, final File targetDir) {
        try {
            final File file = File.createTempFile(name, extension, targetDir);
            file.deleteOnExit();

            // Mark corresponding indices for deletion on exit as well just in case an index is created for the temp file:
            new File(file.getAbsolutePath() + FileExtensions.TRIBBLE_INDEX).deleteOnExit();
            new File(file.getAbsolutePath() + FileExtensions.TABIX_INDEX).deleteOnExit();
            new File(file.getAbsolutePath() + FileExtensions.BAI_INDEX).deleteOnExit();
            new File(file.getAbsolutePath().replaceAll(extension + "$", FileExtensions.BAI_INDEX)).deleteOnExit();
            new File(file.getAbsolutePath() + ".md5").deleteOnExit();
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
