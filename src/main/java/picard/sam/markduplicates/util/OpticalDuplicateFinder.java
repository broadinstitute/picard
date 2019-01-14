/*
 * The MIT License
 *
 * Copyright (c) 2014-2016 The Broad Institute
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
import htsjdk.samtools.util.ProgressLogger;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.ReadNameParser;
import picard.util.GraphUtils;

import java.io.Serializable;
import java.util.*;

/**
 * Contains methods for finding optical/co-localized/sequencing duplicates.
 *
 * @author Tim Fennell
 * @author Nils Homer
 */
public class OpticalDuplicateFinder extends ReadNameParser implements Serializable {

    public int opticalDuplicatePixelDistance;

    public static final int DEFAULT_OPTICAL_DUPLICATE_DISTANCE = 100;
    public static final int DEFAULT_BIG_DUPLICATE_SET_SIZE = 1000;
    public static final int DEFAULT_MAX_DUPLICATE_SET_SIZE = 300000; // larger than this number will generate over 100 billion comparisons in the n^2 algorithm below

    private int bigDuplicateSetSize = DEFAULT_BIG_DUPLICATE_SET_SIZE;
    private long maxDuplicateSetSize = DEFAULT_MAX_DUPLICATE_SET_SIZE;

    /**
     * Sets the size of a set that is big enough to log progress about.
     * Defaults to {@value picard.sam.markduplicates.util.OpticalDuplicateFinder#DEFAULT_BIG_DUPLICATE_SET_SIZE}
     *
     * @param bigDuplicateSetSize the size of a set that is big enough to log progress about
     */
    public void setBigDuplicateSetSize(final int bigDuplicateSetSize) {
        this.bigDuplicateSetSize = bigDuplicateSetSize;
    }

    /**
     * Sets the size of a set that is too big to process.
     * Defaults to {@value picard.sam.markduplicates.util.OpticalDuplicateFinder#DEFAULT_MAX_DUPLICATE_SET_SIZE}
     *
     * @param maxDuplicateSetSize the size of a set that is too big enough to process
     */
    public void setMaxDuplicateSetSize(final long maxDuplicateSetSize) {
        if (maxDuplicateSetSize < 1) {
            this.maxDuplicateSetSize = Long.MAX_VALUE;
        }
        this.maxDuplicateSetSize = maxDuplicateSetSize;
    }

    /**
     * Uses the default duplicate distance {@link OpticalDuplicateFinder#DEFAULT_OPTICAL_DUPLICATE_DISTANCE}
     * ({@value picard.sam.markduplicates.util.OpticalDuplicateFinder#DEFAULT_OPTICAL_DUPLICATE_DISTANCE}) and the default read name regex
     * {@link ReadNameParser#DEFAULT_READ_NAME_REGEX}.
     */
    public OpticalDuplicateFinder() {
        super();
        this.opticalDuplicatePixelDistance = DEFAULT_OPTICAL_DUPLICATE_DISTANCE;
    }

    /**
     * @param readNameRegex                 see {@link ReadNameParser#DEFAULT_READ_NAME_REGEX}.
     * @param opticalDuplicatePixelDistance the optical duplicate pixel distance
     * @param log                           the log to which to write messages.
     */
    public OpticalDuplicateFinder(final String readNameRegex, final int opticalDuplicatePixelDistance, final Log log) {
        super(readNameRegex, log);
        this.opticalDuplicatePixelDistance = opticalDuplicatePixelDistance;
    }

    /**
     * @param readNameRegex                 see {@link ReadNameParser#DEFAULT_READ_NAME_REGEX}.
     * @param opticalDuplicatePixelDistance the optical duplicate pixel distance
     * @param maxDuplicateSetSize           the size of a set that is too big enough to process
     * @param log                           the log to which to write messages.
     */
    public OpticalDuplicateFinder(final String readNameRegex, final int opticalDuplicatePixelDistance, final long maxDuplicateSetSize, final Log log) {
        super(readNameRegex, log);
        this.opticalDuplicatePixelDistance = opticalDuplicatePixelDistance;
        this.maxDuplicateSetSize = maxDuplicateSetSize;
    }

    /**
     * Finds which reads within the list of duplicates that are likely to be optical/co-localized duplicates of
     * one another. Within each cluster of optical duplicates that is found, one read remains un-flagged for
     * optical duplication and the rest are flagged as optical duplicates.  The set of reads that are considered
     * optical duplicates are indicated by returning "true" at the same index in the resulting boolean[] as the
     * read appeared in the input list of physical locations.
     *
     * @param list   a list of reads that are determined to be duplicates of one another
     * @param keeper a single PhysicalLocation that is the one being kept as non-duplicate, and thus should never be
     *               annotated as an optical duplicate. May in some cases be null, or a PhysicalLocation not
     *               contained within the list!
     * @return a boolean[] of the same length as the incoming list marking which reads are optical duplicates
     */
    public boolean[] findOpticalDuplicates(final List<? extends PhysicalLocation> list, final PhysicalLocation keeper) {
        final int length = list.size();
        final boolean[] opticalDuplicateFlags = new boolean[length];

        // If there is only one or zero reads passed in (so there are obviously no optical duplicates),
        // or if there are too many reads (so we don't want to try to run this expensive n^2 algorithm),
        // then just return an array of all false
        if (this.readNameRegex == null || length < 2 || length > maxDuplicateSetSize) {
            return opticalDuplicateFlags;
        }

        final PhysicalLocation actualKeeper = keeperOrNull(list, keeper);

        final Log log;
        final ProgressLogger progressLoggerForKeeper, progressLoggerForRest;
        final boolean logProgress = length > bigDuplicateSetSize;

        if (logProgress) {
            log = Log.getInstance(OpticalDuplicateFinder.class);
            progressLoggerForKeeper = new ProgressLogger(log, 10000, "compared", "ReadEnds to keeper");
            progressLoggerForRest = new ProgressLogger(log, 1000, "compared", "ReadEnds to others");

            log.info("Large duplicate set. size = " + length);
            log.debug("About to compare to keeper:" + actualKeeper);
        } else {
            log = null;
            progressLoggerForKeeper = null;
            progressLoggerForRest = null;
        }
        if (length >= (keeper == null ? 3 : 4)) {
            return getOpticalDuplicatesFlagWithGraph(list, actualKeeper, opticalDuplicateFlags, log, progressLoggerForKeeper, progressLoggerForRest, logProgress);
        } else {
            return getOpticalDuplicatesFlagFast(list, actualKeeper, opticalDuplicateFlags, log, progressLoggerForKeeper, progressLoggerForRest, logProgress);
        }
    }

    /**
     * Compute optical duplicates quickly in the standard case where we know that there won't be any transitive distances to worry about.
     *
     * Note, this is guaranteed to be correct when there are at most 2x reads from a readgroup or 3x with the keeper present
     */
    private boolean[] getOpticalDuplicatesFlagFast(List<? extends PhysicalLocation> list, PhysicalLocation actualKeeper, boolean[] opticalDuplicateFlags, Log log, ProgressLogger progressLoggerForKeeper, ProgressLogger progressLoggerForRest, boolean logProgress) {
        final int length = list.size();

        // First go through and compare all the reads to the keeper
        if (actualKeeper != null) {
            for (int i = 0; i < length; ++i) {
                final PhysicalLocation other = list.get(i);
                opticalDuplicateFlags[i] = closeEnough(actualKeeper, other, this.opticalDuplicatePixelDistance);
                // The main point of adding this log and if statement (also below) is a workaround a bug in the JVM
                // which causes a deep exception (https://github.com/broadinstitute/picard/issues/472).
                // It seems that this is related to https://bugs.openjdk.java.net/browse/JDK-8033717 which
                // was closed due to non-reproducibility. We came across a bam file that evoked this error
                // every time we tried to duplicate-mark it. The problem seemed to be a duplicate-set of size 500,000,
                // and this loop seemed to kill the JVM for some reason. This logging statement (and the one in the
                // loop below) solved the problem.
            }
        }
        if (logProgress) log.debug("Done with comparing to keeper, now the rest.");
        // Now go through and do each pairwise comparison not involving the actualKeeper
        for (int i = 0; i < length; ++i) {
            final PhysicalLocation lhs = list.get(i);
            if (lhs == actualKeeper) continue; // no comparisons to actualKeeper since those are all handled above
            // logging here for same reason as above
            if (logProgress) progressLoggerForRest.record(String.format("%d", lhs.getReadGroup()), lhs.getX());
            for (int j = i + 1; j < length; ++j) {
                final PhysicalLocation rhs = list.get(j);
                if (rhs == actualKeeper) continue; // no comparisons to actualKeeper since those are all handled above
                if (opticalDuplicateFlags[i] && opticalDuplicateFlags[j])
                    continue; // both already marked, no need to check
                if (closeEnough(lhs, rhs, this.opticalDuplicatePixelDistance)) {
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

    /**
     * Compute the optical duplicates correctly in the case where the duplicate group could end up with transitive optical duplicates
     */
    private boolean[] getOpticalDuplicatesFlagWithGraph(List<? extends PhysicalLocation> list, PhysicalLocation keeper, boolean[] opticalDuplicateFlags, Log log, ProgressLogger progressLoggerForKeeper, ProgressLogger progressLoggerForRest, boolean logProgress) {
        // Make a graph where the edges are reads that lie within the optical duplicate pixel distance from each other,
        // we will then use the union-find algorithm to cluster the graph and find optical duplicate groups
        final GraphUtils.Graph<Integer> opticalDistanceRelationGraph = new GraphUtils.Graph<>();
        if (logProgress) {
            log.debug("Building adjacency graph for duplicate group");
        }

        final Map<Integer, List<Integer>> tileRGmap = new HashMap<>();

        int keeperIndex = -1;
        for (int i = 0; i < list.size(); i++) {
            PhysicalLocation currentLoc = list.get(i);
            if (currentLoc == keeper) {
                keeperIndex = i;
            }
            if (currentLoc.hasLocation()) {
                final int key = ((int) currentLoc.getReadGroup() << 16) + currentLoc.getTile();

                if (tileRGmap.containsKey(key)) {
                    tileRGmap.get(key).add(i);
                } else {
                    final List<Integer> pLocation = new ArrayList<>();
                    pLocation.add(i);
                    tileRGmap.put(key, pLocation);
                }
            }
            opticalDistanceRelationGraph.addNode(i);
        }

        // Since because finding adjacent optical duplicates is an O(n^2) operation, we can subdivide the input into its
        // readgroups in order to reduce the amount of redundant checks across readgroups between reads.
        for (List<Integer> tileGroup : tileRGmap.values()) {
            if (tileGroup.size() > 1) {
                fillGraphFromAGroup(list, tileGroup, logProgress, progressLoggerForKeeper,  this.opticalDuplicatePixelDistance, opticalDistanceRelationGraph);
            }
        }

        if (logProgress) {
            log.debug("Finished building adjacency graph for duplicate group, moving onto clustering");
        }

        // Keep a map of the reads and their cluster assignments
        final Map<Integer, Integer> opticalDuplicateClusterMap = opticalDistanceRelationGraph.cluster();
        final Map<Integer, Integer> clusterToRepresentativeRead = new HashMap<>();
        Integer keeperCluster = null;

        // Specially mark the keeper as specifically not a duplicate if it exists
        if (keeperIndex >= 0) {
            clusterToRepresentativeRead.put(opticalDuplicateClusterMap.get(keeperIndex), keeperIndex);
            keeperCluster = opticalDuplicateClusterMap.get(keeperIndex);
        }

        for (final Map.Entry<Integer, Integer> entry : opticalDuplicateClusterMap.entrySet()) {
            // logging here for same reason as above
            final int recordIndex = entry.getKey();
            final int recordAssignedCluster = entry.getValue();
            if (logProgress) {
                progressLoggerForRest.record(String.format("%d", list.get(recordIndex).getReadGroup()), list.get(recordIndex).getX());
            }

            // If its not the first read we've seen for this cluster, mark it as an optical duplicate
            if (clusterToRepresentativeRead.containsKey(recordAssignedCluster) && recordIndex != keeperIndex) {
                final PhysicalLocation representativeLoc = list.get(clusterToRepresentativeRead.get(recordAssignedCluster));
                final PhysicalLocation currentRecordLoc = list.get(recordIndex);

                // If not in the keeper cluster, then keep the minX -> minY valued duplicate (note the tile must be equal for reads to cluster together)
                if (!(keeperIndex >= 0 && recordAssignedCluster == keeperCluster) && // checking we don't accidentally set the keeper as an optical duplicate
                        (currentRecordLoc.getX() < representativeLoc.getX() || (currentRecordLoc.getX() == representativeLoc.getX() && currentRecordLoc.getY() < representativeLoc.getY()))) {
                    // Mark the old min as an optical duplicate, and save the new min
                    opticalDuplicateFlags[clusterToRepresentativeRead.get(recordAssignedCluster)] = true;
                    clusterToRepresentativeRead.put(recordAssignedCluster, recordIndex);
                } else {
                    // If a smaller read has already been visited, mark the test read as an optical duplicate
                    opticalDuplicateFlags[recordIndex] = true;
                }
            } else {
                clusterToRepresentativeRead.put(recordAssignedCluster, recordIndex);
            }
        }

        return opticalDuplicateFlags;
    }

    private void fillGraphFromAGroup(final List<? extends PhysicalLocation> wholeList, final List<Integer> groupList, final boolean logProgress, final ProgressLogger progressLoggerForKeeper, final int distance, final GraphUtils.Graph<Integer> opticalDistanceRelationGraph) {

        for (int i = 0; i < groupList.size(); i++) {
            final int iIndex = groupList.get(i);
            final PhysicalLocation currentLoc = wholeList.get(iIndex);
            // The main point of adding this log and if statement (also below) is a workaround a bug in the JVM
            // which causes a deep exception (https://github.com/broadinstitute/picard/issues/472).
            // It seems that this is related to https://bugs.openjdk.java.net/browse/JDK-8033717 which
            // was closed due to non-reproducibility. We came across a bam file that evoked this error
            // every time we tried to duplicate-mark it. The problem seemed to be a duplicate-set of size 500,000,
            // and this loop seemed to kill the JVM for some reason. This logging statement (and the one in the
            // loop below) solved the problem.
            if (logProgress) {
                progressLoggerForKeeper.record(String.format("%d", currentLoc.getReadGroup()), currentLoc.getX());
            }

            for (int j = i + 1; j < groupList.size(); j++) {
                final int jIndex = groupList.get(j);
                final PhysicalLocation other = wholeList.get(jIndex);

                if (closeEnoughShort(currentLoc, other, distance)) {
                    opticalDistanceRelationGraph.addEdge(iIndex, jIndex);
                }
            }
        }
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

    private boolean closeEnoughShort(final PhysicalLocation lhs, final PhysicalLocation rhs, final int distance) {
        return lhs != rhs &&
                Math.abs(lhs.getX() - rhs.getX()) <= distance &&
                Math.abs(lhs.getY() - rhs.getY()) <= distance;
    }
}
