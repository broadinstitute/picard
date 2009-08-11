/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

package net.sf.picard.metrics;

/**
 * Collector that aggregates a number of other collectors.
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class AggregateMetricCollector<T extends MetricBase, R> implements MetricCollector<T, R> {
    private final MetricCollector<T, R>[] collectors;

    public AggregateMetricCollector(final MetricCollector<T, R>... collectors) {
        if (collectors.length == 0) {
            throw new IllegalArgumentException("Must supply at least one collector.");
        }
        this.collectors = collectors;
    }

    public void addRecord(final R record) {
        for (final MetricCollector<T, R> collector : this.collectors) {
            collector.addRecord(record);
        }
    }

    public void onComplete() {
        for (final MetricCollector<T, R> collector : this.collectors) {
            collector.onComplete();
        }
    }

    public void setMetrics(final T metrics) {
        for (final MetricCollector<T, R> collector : this.collectors) {
            collector.setMetrics(metrics);
        }
    }
    
    public T getMetrics() {
        return this.collectors[0].getMetrics();
    }
}