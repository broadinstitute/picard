package picard.cmdline.argumentcollections;

import java.io.File;

/**
 * Base interface for an interval argument collection.
 */
public interface IntervalArgumentCollection {
    /**
     * @return The interval file provided by the user, if any, or the default, if any.
     */
    File getIntervalFile();
}
