package net.sf.samtools;

import java.util.ArrayList;
import java.util.List;

/**
 * The linear index associated with a given reference in a BAM index.
 *
 * @author mhanna
 * @version 0.1
 */
public class LinearIndexEntry {

    public final int window;
    public final long offset;

    public LinearIndexEntry(int window, long offset) {
        this.window = window;
        this.offset = offset;
    }
}
