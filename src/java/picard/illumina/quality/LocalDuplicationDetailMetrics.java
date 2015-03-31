package picard.illumina.quality;

import htsjdk.samtools.metrics.MetricBase;
import java.lang.Override;
import java.lang.String;
import java.util.List;

/** a metric class for describing a single bunch of local duplicates **/
public class LocalDuplicationDetailMetrics extends MetricBase {
    /** The Tile that is described by this metric. */
    public Integer TILE;

    /** The sequence of bases common to duplicates in this bunch. */
    public String BASES;

    /** The number of reads (clusters) in this bunch. */
    public int SIZE;

    /**All the points in this bunch in a space-free format x1,y1;x2,y2; etc. */
    public String POINTS_STRING;

    public LocalDuplicationDetailMetrics(final Integer tile, final String bases, final List<Point> points) {
        TILE = tile;
        BASES = bases;
        SIZE = points.size();

        StringBuilder builder = new StringBuilder();
        for (final Point p : points) {
            builder.append(p.getX());
            builder.append(',');
            builder.append(p.getY());
            builder.append(';');
        }
        POINTS_STRING = builder.toString();
    }

    /** This constructor is necessary for reading metrics from file */
    public LocalDuplicationDetailMetrics() { }
}