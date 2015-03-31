package picard.illumina.quality;

import htsjdk.samtools.metrics.MetricBase;
import java.lang.Override;
import java.lang.String;

public class LocalDuplicationSummaryMetrics extends MetricBase {
    /** The Tile that is described by this metric. Can be a string (like "All") to mean some marginal over tiles. * */
    public String TILE = null;

    /** The total number of PF reads on this tile. */
    public long READS = 0;

    /** Local duplicates in this tile.  In a bunch of N clusters, N - 1 are duplicates. */
    public long LOCAL_DUPLICATES = 0;

    /** The rate (not the percentage!) of local duplication. */
    public double PCT_LOCAL_DUPLICATES = 0.0;

    public LocalDuplicationSummaryMetrics(final String tile) {
        TILE = tile;
    }

    /** This constructor is necessary for reading metrics from file. */
    public LocalDuplicationSummaryMetrics() { }

    public void merge(final LocalDuplicationSummaryMetrics metric) {
        READS += metric.READS;
        LOCAL_DUPLICATES += metric.LOCAL_DUPLICATES;
    }

    public void calculateDerivedFields() {
        if (READS != 0) PCT_LOCAL_DUPLICATES = ((double) LOCAL_DUPLICATES) / READS;
    }
}
