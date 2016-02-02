/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates.util;

import htsjdk.samtools.util.Log;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.ReadNameParser;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Contains methods for finding optical/co-localized/sequencing duplicates.
 *
 * @author Tim Fennell
 * @author Nils Homer
 */
public class OpticalDuplicateFinder extends ReadNameParser {

    public int opticalDuplicatePixelDistance;

    public static final int DEFAULT_OPTICAL_DUPLICATE_DISTANCE = 100;

    /**
     * Uses the default duplicate distance {@value DEFAULT_OPTICAL_DUPLICATE_DISTANCE} and the default read name regex
     * {@link ReadNameParser#DEFAULT_READ_NAME_REGEX}.
     */
    public OpticalDuplicateFinder() {
        super();
        this.opticalDuplicatePixelDistance = DEFAULT_OPTICAL_DUPLICATE_DISTANCE;
    }

    /**
     *
     * @param readNameRegex see {@link ReadNameParser#DEFAULT_READ_NAME_REGEX}.
     * @param opticalDuplicatePixelDistance the optical duplicate pixel distance
     * @param log the log to which to write messages.
     */
    public OpticalDuplicateFinder(final String readNameRegex, final int opticalDuplicatePixelDistance, final Log log) {
        super(readNameRegex, log);
        this.opticalDuplicatePixelDistance = opticalDuplicatePixelDistance;
    }

    /**
     * Finds which reads within the list of duplicates that are likely to be optical/co-localized duplicates of
     * one another. Within each cluster of optical duplicates that is found, one read remains un-flagged for
     * optical duplication and the rest are flagged as optical duplicates.  The set of reads that are considered
     * optical duplicates are indicated by returning "true" at the same index in the resulting boolean[] as the
     * read appeared in the input list of physical locations.
     *
     * @param list a list of reads that are determined to be duplicates of one another
     * @param keeper a single PhysicalLocation that is the one being kept as non-duplicate, and thus should never be
     *               annotated as an optical duplicate. May in some cases be null, or a PhysicalLocation not
     *               contained within the list!
     * @return a boolean[] of the same length as the incoming list marking which reads are optical duplicates
     */
    public boolean[] findOpticalDuplicates(final List<? extends PhysicalLocation> list, final PhysicalLocation keeper) {
        // If there is only one or zero reads passed in, then just return an array of all false
        if (list.size() < 2) return new boolean[list.size()];

        final int length = list.size();
        final boolean[] opticalDuplicateFlags = new boolean[length];
        final int distance = this.opticalDuplicatePixelDistance;

        final PhysicalLocation actualKeeper = keeperOrNull(list, keeper);

        // First go through and compare all the reads to the keeper
        if (actualKeeper != null) {
            for (int i=0; i<length; ++i) {
                final PhysicalLocation other = list.get(i);
                opticalDuplicateFlags[i] = closeEnough(actualKeeper, other, distance);
            }
        }

        // Now go through and do each pairwise comparison not involving the actualKeeper
        for (int i=0; i<length; ++i) {
            final PhysicalLocation lhs = list.get(i);
            if (lhs == actualKeeper) continue; // no comparisons to actualKeeper since those are all handled above

            for (int j =i+1; j<length; ++j) {
                final PhysicalLocation rhs = list.get(j);
                if (rhs == actualKeeper) continue; // no comparisons to actualKeeper since those are all handled above
                if (opticalDuplicateFlags[i] && opticalDuplicateFlags[j]) continue; // both already marked, no need to check

                if (closeEnough(lhs, rhs, distance)) {
                    // At this point we want to mark either lhs or rhs as duplicate. Either could have been marked
                    // as a duplicate of the keeper (but not both - that's checked above), so be careful about which
                    // one to now mark as a duplicate.
                    final int index = opticalDuplicateFlags[j] ? i : j;
                    opticalDuplicateFlags[index] = true;
                }
            }
        }
        
        return opticalDuplicateFlags;
    }

    /** Returns the keeper if it is contained within the list and has location information, otherwise null. */
    private PhysicalLocation keeperOrNull(final List<? extends PhysicalLocation> list, final PhysicalLocation keeper) {
        if (keeper != null && keeper.hasLocation()) {
            for (final PhysicalLocation loc : list) {
                if (loc == keeper) return keeper;
            }
        }
        return null;
    }

    /** Simple method to test whether two physical locations are close enough to each other to be deemed optical dupes. */
    private boolean closeEnough(final PhysicalLocation lhs, final PhysicalLocation rhs, final int distance) {
        return lhs != rhs &&                                    // no comparing an object to itself (checked using object identity)!
               lhs.hasLocation() && rhs.hasLocation() &&        // no comparing objects without locations
               lhs.getReadGroup() == rhs.getReadGroup() &&      // must be in the same RG to be optical duplicates
               lhs.getTile()      == rhs.getTile()      &&      // and the same tile
               Math.abs(lhs.getX() - rhs.getX()) <= distance &&
               Math.abs(lhs.getY() - rhs.getY()) <= distance;
    }
}
