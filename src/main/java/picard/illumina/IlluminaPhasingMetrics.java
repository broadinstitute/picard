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

package picard.illumina;

import htsjdk.samtools.metrics.MetricBase;
import picard.illumina.parser.Tile;
import picard.illumina.parser.TileTemplateRead;

import java.util.ArrayList;
import java.util.Collection;


/**
 * Metrics for Illumina Basecalling that stores median phasing and prephasing percentages on a per-template-read, per-lane basis.
 * Phasing refers to the fraction of molecules that fall behind or jump ahead (prephasing) during a read cycle.
 * For each lane/template read # (i.e. FIRST, SECOND) combination we will store the median values of both the phasing and prephasing
 * values for every tile in that lane/template read pair.
 *
 * @author jgentry
 */
public class IlluminaPhasingMetrics extends MetricBase {
    /** Illumina flowcell lane number */
    public long LANE;
    /** Defines an Illumina template read number (first or second)  */
    public String TYPE_NAME;
    /** Median phasing value across all tiles in a lane, applied to the first and second template reads */
    public double PHASING_APPLIED;
    /** Median pre-phasing value across all tiles in a lane, applied to the first and second template reads */
    public double PREPHASING_APPLIED;
    /** Calculate the median phasing & prephasing values for a lane's tiles and create the appropriate IlluminaPhasingMetrics for them */
    public static Collection<IlluminaPhasingMetrics> getPhasingMetricsForTiles(final long lane, final Collection<Tile> tilesForLane, final boolean usePercentage) {
        final LanePhasingMetricsCollector lanePhasingMetricsCollector = new LanePhasingMetricsCollector(tilesForLane, usePercentage);
        final Collection<IlluminaPhasingMetrics> phasingMetrics = new ArrayList<IlluminaPhasingMetrics>();
        for (final TileTemplateRead tileTemplateRead : lanePhasingMetricsCollector.getMedianPhasingMap().keySet()) {
            final IlluminaPhasingMetrics phasingMetric = new IlluminaPhasingMetrics();
            phasingMetric.LANE = lane;
            phasingMetric.TYPE_NAME = tileTemplateRead.toString();
            phasingMetric.PHASING_APPLIED = lanePhasingMetricsCollector.getMedianPhasingMap().get(tileTemplateRead);
            phasingMetric.PREPHASING_APPLIED = lanePhasingMetricsCollector.getMedianPrePhasingMap().get(tileTemplateRead);
            phasingMetrics.add(phasingMetric);
        }

        return phasingMetrics;
    }

    /** This property is not exposed in a field to avoid complications with MetricBase's dependency on reflection. */
    public static String getExtension() {
        return "illumina_phasing_metrics";
    }
}
