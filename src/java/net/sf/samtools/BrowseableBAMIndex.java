package net.sf.samtools;

/**
 * An index interface with additional functionality for querying and inspecting the structure of a BAM index.
 *
 * @author mhanna
 * @version 0.1
 */
public interface BrowseableBAMIndex extends BAMIndex {

    /**
     * Gets the size (number of bins in) a given level of a BAM index.
     * @param levelNumber Level for which to inspect the size.
     * @return Size of the given level.
     */
    public int getLevelSize(final int levelNumber);

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin  for which to determine the level.
     * @return the level associated with the given bin number.
     */
    public int getLevelForBin(final Bin bin);
    
    /**
     * Gets the first locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    int getFirstLocusInBin(final Bin bin);

    /**
     * Gets the last locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    int getLastLocusInBin(final Bin bin);

    /**
     * Get a list of bins in the BAM file that may contain SAMRecords for the given range.
     * @param referenceIndex sequence of desired SAMRecords
     * @param startPos 1-based start of the desired interval, inclusive
     * @param endPos 1-based end of the desired interval, inclusive
     * @return a list of bins that contain relevant data.
     */
    BinList getBinsOverlapping(final int referenceIndex, final int startPos, final int endPos);

    /**
     * Perform an overlapping query of all bins bounding the given location.
     * @param bin The bin over which to perform an overlapping query.
     * @return The file pointers
     */
    BAMFileSpan getSpanOverlapping(final Bin bin);    
}
