package picard.illumina;

import htsjdk.samtools.metrics.MetricBase;
import java.lang.Double;
import java.lang.Long;
import java.lang.String;

/**
 * Embodies characteristics that describe a lane.  
 * @author mccowan
 */
public class IlluminaLaneMetrics extends MetricBase {
    /** The number of clusters per unit area on the this lane expressed in units of [cluster / mm^2]. */
    public Double CLUSTER_DENSITY;
    
    /** This lane's number. */
    public Long LANE;
    
    /** This property is not exposed in a field to avoid complications with MetricBase's dependency on reflection. */
    public static String getExtension() {
        return "illumina_lane_metrics";
    }
}
