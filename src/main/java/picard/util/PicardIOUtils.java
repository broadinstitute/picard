package picard.util;

import htsjdk.samtools.SamReader;

import java.nio.file.Path;

/**
 * These utility methods should probably go to htsjdk, but implement here temporarily in order to avoid waiting for
 * picard-htsjdk synchronization.
 */
public class PicardIOUtils {

    public static boolean isBamFile(final Path path) {
        // tsato: is path.toString correct? Or getName()?
        return ((path != null) && SamReader.Type.BAM_TYPE.hasValidFileExtension(path.toString()));
    }
}
