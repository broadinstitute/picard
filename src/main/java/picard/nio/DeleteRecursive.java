package picard.nio;

import htsjdk.samtools.util.Log;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;

/**
 *
 * Adapted from GATK.
 *
 * Class to hold a set of {@link Path} to be delete on the JVM exit through a shutdown hook.
 *
 * <p>This class is a modification of {@link htsjdk.samtools.util.nio.DeleteOnExitPathHook}
 *
 * This class should be considered an implementation detail of {@link GATKIOUtils#deleteOnExit(Path)} and not used directly.
 */
class DeleteRecursivelyOnExitPathHook {
    private static final Log LOG = Log.getInstance(DeleteRecursivelyOnExitPathHook.class);
    private static LinkedHashSet<Path> paths = new LinkedHashSet<>();
    static {
        Runtime.getRuntime().addShutdownHook(new Thread(DeleteRecursivelyOnExitPathHook::runHooks));
    }

    private DeleteRecursivelyOnExitPathHook() {}

    /**
     * Adds a {@link Path} for deletion on JVM exit.
     *
     * @param path path to be deleted.  This path may be a non-empty directory and the entire directory structure will
     *            be deleted.
     *
     * @throws IllegalStateException if the shutdown hook is in progress.
     */
    public static synchronized void add(Path path) {
        if(paths == null) {
            // DeleteOnExitHook is running. Too late to add a file
            throw new IllegalStateException("Shutdown in progress");
        }

        paths.add(path);
    }

    static void runHooks() {
        LinkedHashSet<Path> thePaths;

        synchronized (DeleteRecursivelyOnExitPathHook.class) {
            thePaths = paths;
            paths = null;
        }

        ArrayList<Path> toBeDeleted = new ArrayList<>(thePaths);

        // reverse the list to maintain previous jdk deletion order.
        // Last in first deleted.
        Collections.reverse(toBeDeleted);
        for (Path path : toBeDeleted) {
            try {
                PicardIOUtils.deleteRecursively(path);
            } catch (final Exception e) {
                // do nothing if itcannot be deleted, because it is a shutdown hook
                LOG.debug(e, "Could not recursively delete ", path.toString(), " during JVM shutdown because we encountered the following exception:");
            }
        }
    }
}
