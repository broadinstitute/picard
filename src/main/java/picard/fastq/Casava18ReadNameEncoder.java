package picard.fastq;

import picard.illumina.parser.ClusterData;

/**
 * A read name encoder conforming to the standard described by Illumina Casava 1.8.
 * 
 * @see <a href="http://biowulf.nih.gov/apps/CASAVA1_8_Changes.pdf">Casava 1.8 update</a>
 * @author mccowan
 */
public class Casava18ReadNameEncoder implements ReadNameEncoder {
    private static final char CONTROL_FIELD_VALUE = '0';
    private static final char SEPARATOR = ':';

    /* The following is to make generation of string representations of integers fast for a small subset of ints. */
    private static final int INT_CACHE_LIMIT = 5000;
    private static final String[] INT_STRINGS = new String[INT_CACHE_LIMIT];
    static {
        for (int i = 0; i < INT_STRINGS.length; ++i) {
            INT_STRINGS[i] = Integer.toString(i);
        }
    }

    private final String nameBase;
    private int bufferSize;

    public Casava18ReadNameEncoder(final String instrumentName, final String runId, final String flowcellId) {
        this.nameBase           = instrumentName + SEPARATOR + runId + SEPARATOR + flowcellId + SEPARATOR;
        this.bufferSize     = this.nameBase.length();
    }

    /** Converts an int to a string, with cached results for some ints. */
    private static String encodeInt(final int i) {
        if (i >= 0 && i < INT_CACHE_LIMIT) return INT_STRINGS[i];
        else return Integer.toString(i);
    }

    @Override
    public String generateReadName(final ClusterData cluster, final Integer pairNumber) {
        final StringBuilder builder = new StringBuilder(bufferSize);
        builder.append(this.nameBase);
        builder.append(encodeInt(cluster.getLane()));
        builder.append(SEPARATOR);
        builder.append(encodeInt(cluster.getTile()));
        builder.append(SEPARATOR);
        builder.append(encodeInt(cluster.getX()));
        builder.append(SEPARATOR);
        builder.append(encodeInt(cluster.getY()));

        builder.append(' ');

        if (pairNumber != null) builder.append(encodeInt(pairNumber));
        builder.append(SEPARATOR);
        builder.append(cluster.isPf() ? 'N' : 'Y');  // encoded in read name as Y == fails filter
        builder.append(SEPARATOR);
        builder.append(CONTROL_FIELD_VALUE);
        builder.append(SEPARATOR);
        if (cluster.getMatchedBarcode() != null) builder.append(cluster.getMatchedBarcode());

        // Update the buffer size so that next time through we don't need to grow the builder
        if (builder.length() > this.bufferSize) {
            this.bufferSize = builder.length();
        }

        return builder.toString();
    }

    @Override
    public String generateShortName(ClusterData cluster) {
        final StringBuilder builder = new StringBuilder(bufferSize);
        builder.append(this.nameBase);
        builder.append(encodeInt(cluster.getLane()));
        builder.append(SEPARATOR);
        builder.append(encodeInt(cluster.getTile()));
        builder.append(SEPARATOR);
        builder.append(encodeInt(cluster.getX()));
        builder.append(SEPARATOR);
        builder.append(encodeInt(cluster.getY()));

        return builder.toString();
    }
}
