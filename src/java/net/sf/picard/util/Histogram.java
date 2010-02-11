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

package net.sf.picard.util;

import net.sf.picard.util.Histogram.Bin;

import java.util.*;

/**
 * Class for computing and accessing histogram type data.  Stored internally in
 * a sorted Map so that keys can be iterated in order.
 *
 * @author Tim Fennell
 */
public class Histogram<K extends Comparable> extends TreeMap<K, Bin> {
    private String binLabel   = "BIN";
    private String valueLabel = "VALUE";
    private Double mean;

    /** Constructs a new Histogram with default bin and value labels. */
    public Histogram() { }

    /** Constructs a new Histogram with supplied bin and value labels. */
    public Histogram(final String binLabel, final String valueLabel) {
        this.binLabel = binLabel;
        this.valueLabel = valueLabel;
    }

    /** Constructs a new Histogram that'll use the supplied comparator to sort keys. */
    public Histogram(final Comparator<K> comparator) {
        super(comparator);
    }

    /** Constructor that takes labels for the bin and values and a comparator to sort the bins. */
    public Histogram(final String binLabel, final String valueLabel, final Comparator<K> comparator) {
        this(comparator);
        this.binLabel = binLabel;
        this.valueLabel = valueLabel;
    }

    /** Copy constructor for a histogram. */
    public Histogram(final Histogram<K> in) {
        this.binLabel = in.binLabel;
        this.valueLabel = in.valueLabel;
        this.mean = in.mean;

        for (final Bin bin : in.values()) {
            increment(bin.id, bin.value);
        }
    }

    /** Represents a bin in the Histogram. */
    public class Bin {
        private final K id;
        private double value = 0;

        /** Constructs a new bin with the given ID. */
        private Bin(final K id) { this.id = id; }

        /** Gets the ID of this bin. */
        public K getId() { return id; }

        /** Gets the value in the bin. */
        public double getValue() { return value; }

        /** Returns the String format for the value in the bin. */
        public String toString() { return String.valueOf(this.value); }

        /** Checks the equality of the bin by ID and value. */
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final Bin bin = (Bin) o;

            if (Double.compare(bin.value, value) != 0) return false;
            if (!id.equals(bin.id)) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result;
            final long temp;
            result = id.hashCode();
            temp = value != +0.0d ? Double.doubleToLongBits(value) : 0L;
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }

        public double getIdValue() {
            if (id instanceof Number) {
                return ((Number) id).doubleValue();
            } else {
                throw new UnsupportedOperationException("getIdValue only supported for Histogram<? extends Number>");
            }
        }
    }

    /** Prefill the histogram with the supplied set of bins. */
    public void prefillBins(final K... ids) {
        for (final K id : ids) {
            put(id, new Bin(id));
        }
    }

    /** Increments the value in the designated bin by 1. */
    public void increment(final K id) {
        increment(id, 1d);
    }

    /** Increments the value in the designated bin by the supplied increment. */
    public void increment(final K id, final double increment) {
        Bin bin = get(id);
        if (bin == null) {
            bin = new Bin(id);
            put(id, bin);
        }

        bin.value += increment;
        mean = null;
    }

    public String getBinLabel() { return binLabel; }
    public void setBinLabel(final String binLabel) { this.binLabel = binLabel; }

    public String getValueLabel() { return valueLabel; }
    public void setValueLabel(final String valueLabel) { this.valueLabel = valueLabel; }

    /** Checks that the labels and values in the two histograms are identical. */
    public boolean equals(final Object o) {
        return o != null &&
                (o instanceof Histogram) &&
                ((Histogram) o).binLabel.equals(this.binLabel) &&
                ((Histogram) o).valueLabel.equals(this.valueLabel) &&
                super.equals(o);
    }

    public double getMean() {
        if (mean == null) {
            mean = getSum() / getCount();
        }

        return mean;
    }

    /**
     * Returns the sum of the products of the histgram bin ids and the number of entries in each bin.
     */
    public double getSum() {
        double total = 0;
        for (final Bin bin : values()) {
            total += bin.getValue() * bin.getIdValue();
        }

        return total;
    }

    /**
     * Returns the sum of the number of entries in each bin.
     */
    public double getSumOfValues() {
        double total = 0;
        for (final Bin bin : values()) {
            total += bin.getValue();
        }

        return total;
    }

    public double getStandardDeviation() {
        double total = 0;
        for (final Bin bin : values()) {
            total += bin.getValue() * bin.getIdValue() * bin.getIdValue();
        }

        return Math.sqrt((total / getCount()) - (getMean() * getMean()));
    }

    /**
     * Gets the bin in which the given percentile falls.
     *
     * @param percentile a value between 0 and 1
     * @return the bin value in which the percentile falls
     */
    public double getPercentile(double percentile) {
        if (percentile <= 0) throw new IllegalArgumentException("Cannot query percentiles of 0 or below");
        if (percentile >= 1) throw new IllegalArgumentException("Cannot query percentiles of 1 or above");

        double total = getCount();
        double sofar = 0;
        for (Bin bin : values()) {
            sofar += bin.getValue();
            if (sofar / total >= percentile) return bin.getIdValue();
        }

        throw new IllegalStateException("Could not find percentile: " + percentile);
    }

    /**
     * Returns the cumulative probability of observing a value <= v when sampling the
     * distribution represented by this histogram.
     */
    public double getCumulativeProbability(final double v) {
        double count = 0;
        double total = 0;

        for (final Bin bin : values()) {
            final double binValue = bin.getIdValue();
            if (binValue <= v) count += bin.getValue();
            total += bin.getValue();
        }

        return count / total;
    }

    public double getMedian() {
        double total = 0;
        final double halfCount = getCount() / 2;
        for (final Bin bin : values()) {
            total += bin.getValue();
            if (total >= halfCount) {
                return bin.getIdValue();
            }
        }
        return 0;
    }

    /** Returns id of the Bin that's the mode of the distribution (i.e. the largest bin). */
    public double getMode() {

        return getModeBin().getIdValue();
    }

    /** Returns the Bin that's the mode of the distribution (i.e. the largest bin). */
    private Bin getModeBin() {
        Bin modeBin = null;

        for (final Bin bin : values()) {
            if (modeBin == null || modeBin.value < bin.value) {
                modeBin = bin;
            }
        }

        return modeBin;
    }


    public double getMin() {
        return firstEntry().getValue().getIdValue();
    }

    public double getMax() {
        return lastEntry().getValue().getIdValue();
    }

    public double getCount() {
        double count = 0;
        for (final Bin bin : values()) {
            count += bin.value;
        }

        return count;
    }

    /**
     * Trims the histogram when the bins in the tail of the distribution contain fewer than mode/tailLimit items
     */
    public void trim(final int tailLimit) {
        if (isEmpty()) {
            return;
        }

        final Bin modeBin = getModeBin();
        final double mode = modeBin.getIdValue();
        final double sizeOfModeBin = modeBin.getValue();
        final double minimumBinSize = sizeOfModeBin/tailLimit;
        Histogram<K>.Bin lastBin = null;

        final List<K> binsToKeep = new ArrayList<K>();
        for (Histogram<K>.Bin bin : values()) {
            double binId = ((Number)bin.getId()).doubleValue();

            if (binId <= mode) {
                binsToKeep.add(bin.getId());
            }
            else if ((lastBin != null && ((Number)lastBin.getId()).doubleValue() != binId - 1) || bin.getValue() < minimumBinSize) {
                break;
            }
            else {
                binsToKeep.add(bin.getId());
            }
            lastBin = bin;
        }

        final Object keys[] = keySet().toArray();
        for (Object binId : keys) {
            if (!binsToKeep.contains((K)binId)) {
                remove(binId);
            }
        }
    }


}
