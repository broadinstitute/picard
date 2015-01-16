/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
import java.lang.Override;
import java.lang.String;

/***
 * Metric for Illumina Basecalling that stores means and standard deviations on a per-barcode per-lane basis.  Averages
 * and means are taken over all tiles.
 */
public class IlluminaBasecallingMetrics extends MetricBase
{
    /** The lane for which the metrics were calculated. */
    public String LANE;
    /** The barcode sequence for which the metrics were calculated. */
    public String MOLECULAR_BARCODE_SEQUENCE_1;
    /** The barcode name for which the metrics were calculated. */
    public String MOLECULAR_BARCODE_NAME;
    /** The total number of bases assigned to the index. */
    public long TOTAL_BASES;
    /** The total number of passing-filter bases assigned to the index. */
    public long PF_BASES;
    /** The total number of reads assigned to the index. */
    public long TOTAL_READS;
    /** The total number of passing-filter reads assigned to the index. */
    public long PF_READS;
    /** The total number of clusters assigned to the index. */
    public long TOTAL_CLUSTERS;
    /** The total number of PF clusters assigned to the index. */
    public long PF_CLUSTERS;
    /** The mean number of clusters per tile. */
    public double MEAN_CLUSTERS_PER_TILE = 0d;
    /** The standard deviation of clusters per tile. */
    public double SD_CLUSTERS_PER_TILE = 0d;
    /** The mean percentage of pf clusters per tile. */
    public double MEAN_PCT_PF_CLUSTERS_PER_TILE = 0d;
    /** The standard deviation in percentage of pf clusters per tile. */
    public double SD_PCT_PF_CLUSTERS_PER_TILE = 0d;
    /** The mean number of pf clusters per tile. */
    public double MEAN_PF_CLUSTERS_PER_TILE = 0d;
    /** The standard deviation in number of pf clusters per tile. */
    public double SD_PF_CLUSTERS_PER_TILE = 0d;

    @Override
    public String toString() {
        return String.format("IlluminaBasecallingMetric(Lane:%s,Barcode:%s,Name:%s,MEAN_CLUSTERS_PER_TILE:%s,SD_CLUSTERS_PER_TILE:%s," +
                              "MEAN_PCT_PF_CLUSTERS_PER_TILE:%s,SD_PCT_PF_CLUSTERS_STD_PER_TILE:%s," +
                              "MEAN_PF_CLUSTERS_PER_TILE:%s,SD_PF_CLUSTERS_PER_TILE:%s",
                LANE
                ,MOLECULAR_BARCODE_SEQUENCE_1
                ,MOLECULAR_BARCODE_NAME
                ,MEAN_CLUSTERS_PER_TILE
                ,SD_CLUSTERS_PER_TILE
                ,MEAN_PCT_PF_CLUSTERS_PER_TILE
                ,SD_PCT_PF_CLUSTERS_PER_TILE
                ,MEAN_PF_CLUSTERS_PER_TILE
                ,SD_PF_CLUSTERS_PER_TILE);
    }
}
