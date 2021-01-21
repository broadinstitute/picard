package picard.fastq;

import picard.illumina.parser.ClusterData;

/**
 * @author mccowan
 */
public interface ReadNameEncoder {
    /**
     * Generates a read name string for the provided cluster. 
     *
     * @param cluster The cluster whose reads are having its name generated
     * @param pairNumber 1 if this is the first of the pair, 2 if it is the second, or null if this not a paired read.
     * @return The read name
     */
    String generateReadName(final ClusterData cluster, final Integer pairNumber);

    /**
     * Generates a short read name that includes a minimal amount of information, this is used primarily for read
     * sorting.
     *
     * @param cluster The cluster to generate the short read name from
     * @return The short read name
     */
    String generateShortName(ClusterData cluster);
}
