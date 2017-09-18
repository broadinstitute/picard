package picard.cmdline.argumentcollections;

import java.io.File;

/**
 * Base interface for a reference argument collection.
 */
public interface ReferenceArgumentCollection {
    /**
     * @return The reference provided by the user, if any, or the default, if any.
     */
    File getReferenceFile();
}
