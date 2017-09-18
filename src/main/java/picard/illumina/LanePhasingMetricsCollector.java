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

import htsjdk.samtools.util.CollectionUtil;
import picard.illumina.parser.Tile;
import picard.illumina.parser.TileTemplateRead;
import picard.util.MathUtil;

import java.util.Collection;
import java.util.Collections;
import java.util.Map;
import java.util.TreeMap;

/** Helper class used to transform tile data for a lane into a collection of IlluminaPhasingMetrics */
public class LanePhasingMetricsCollector {
    private final Map<TileTemplateRead, Float> medianPhasingMap;
    private final Map<TileTemplateRead, Float> medianPrePhasingMap;

    /** Constructor takes a lane's collection of Tiles and calculates the median phasing/prephasing for the
     * first and second (if available) reads
     */
    public LanePhasingMetricsCollector(final Collection<Tile> laneTiles, final boolean usePercentage) {
        final Map<TileTemplateRead, Float> medianPhasingMap = new TreeMap<TileTemplateRead, Float>();
        final Map<TileTemplateRead, Float> medianPrePhasingMap = new TreeMap<TileTemplateRead, Float>();

        final CollectionUtil.MultiMap<TileTemplateRead, Float> phasingValues = new CollectionUtil.MultiMap<TileTemplateRead, Float>();
        final CollectionUtil.MultiMap<TileTemplateRead, Float> prePhasingValues = new CollectionUtil.MultiMap<TileTemplateRead, Float>();

        // Collect the phasing/prephasing values from all of the tiles, sorted by template read #
        for (final Tile tile : laneTiles) {
            for (final TileTemplateRead tileTemplateRead : tile.getPhasingMap().keySet()) {
                phasingValues.append(tileTemplateRead, tile.getPhasingMap().get(tileTemplateRead));
                prePhasingValues.append(tileTemplateRead, tile.getPrePhasingMap().get(tileTemplateRead));
            }
        }

        // Calculate the medians for the collected data
        for (final TileTemplateRead tileTemplateRead : phasingValues.keySet()) {
            medianPhasingMap.put(tileTemplateRead, median(phasingValues.get(tileTemplateRead), usePercentage));
            medianPrePhasingMap.put(tileTemplateRead, median(prePhasingValues.get(tileTemplateRead), usePercentage));
        }

        this.medianPhasingMap = Collections.unmodifiableMap(medianPhasingMap);
        this.medianPrePhasingMap = Collections.unmodifiableMap(medianPrePhasingMap);
    }

    public Map<TileTemplateRead, Float> getMedianPhasingMap() {
        return medianPhasingMap;
    }

    public Map<TileTemplateRead, Float> getMedianPrePhasingMap() {
        return medianPrePhasingMap;
    }


    private static float median(final Collection<Float> phaseValues, boolean usePercentage) {
        final double[] values = new double[phaseValues.size()];
        int i = 0;
        for (Float phaseValue : phaseValues) {
            values[i] = (double)phaseValue;
            i++;
        }
        float median = (float) MathUtil.median(values);
        return usePercentage ? median * 100 : median;
    }

}
