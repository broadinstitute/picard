/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.fingerprint;

import picard.util.ReflectionUtil;

/**
 * A metric class to hold the result of {@link ClusterCrosscheckMetrics} fingerprints.
 *
 * @author Yossi Farjoun
 */

public class ClusteredCrosscheckMetric extends CrosscheckMetric {
    /**
     The cluster identifier to which the two groups within this metric ({@link #LEFT_GROUP_VALUE} and {@link #RIGHT_GROUP_VALUE}) belong.
      */
    public Integer CLUSTER;

    /**
     * The number of different groups that are assigned to this cluster. Should be the same value for all rows within the same
     * {@link #CLUSTER}.
     */
    public Integer CLUSTER_SIZE;

    public ClusteredCrosscheckMetric() {
        super();
    }

    public ClusteredCrosscheckMetric(CrosscheckMetric metric) {
        super();
        ReflectionUtil.copyFromBaseClass(metric, this);
    }
}
