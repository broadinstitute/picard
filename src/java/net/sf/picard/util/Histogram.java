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

import java.util.TreeMap;

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
    public Histogram(String binLabel, String valueLabel) {
        this.binLabel = binLabel;
        this.valueLabel = valueLabel;
    }

    /** Copy constructor for a histogram. */
    public Histogram(Histogram<K> in) {
        this.binLabel = in.binLabel;
        this.valueLabel = in.valueLabel;
        this.mean = in.mean;

        for (Bin bin : in.values()) {
            increment(bin.id, bin.value);
        }
    }

    /** Represents a bin in the Histogram. */
    public class Bin {
        private final K id;
        private double value = 0;

        /** Constructs a new bin with the given ID. */
        private Bin(K id) { this.id = id; }

        /** Gets the ID of this bin. */
        public K getId() { return id; }

        /** Gets the value in the bin. */
        public double getValue() { return value; }

        /** Returns the String format for the value in the bin. */
        public String toString() { return String.valueOf(this.value); }

        /** Checks the equality of the bin by ID and value. */
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Bin bin = (Bin) o;

            if (Double.compare(bin.value, value) != 0) return false;
            if (!id.equals(bin.id)) return false;

            return true;
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
    public void prefillBins(K... ids) {
        for (K id : ids) {
            put(id, new Bin(id));
        }
    }

    /** Increments the value in the designated bin by 1. */
    public void increment(K id) {
        increment(id, 1d);
    }

    /** Increments the value in the designated bin by the supplied increment. */
    public void increment(K id, double increment) {
        Bin bin = get(id);
        if (bin == null) {
            bin = new Bin(id);
            put(id, bin);
        }

        bin.value += increment;
        mean = null;
    }

    public String getBinLabel() { return binLabel; }
    public void setBinLabel(String binLabel) { this.binLabel = binLabel; }

    public String getValueLabel() { return valueLabel; }
    public void setValueLabel(String valueLabel) { this.valueLabel = valueLabel; }

    /** Checks that the labels and values in the two histograms are identical. */
    public boolean equals(Object o) {
        return o != null &&
                (o instanceof Histogram) &&
                ((Histogram) o).binLabel.equals(this.binLabel) &&
                ((Histogram) o).valueLabel.equals(this.valueLabel) &&
                super.equals(o);
    }

    public double getMean() {
        if (mean == null) {
            double total = 0;
            for (Bin bin : values()) {
                total += bin.getValue() * bin.getIdValue();
            }
    
            mean = total / getCount();
        }
        
        return mean;
    }
    
    public double getStandardDeviation() {
        double total = 0;
        for (Bin bin : values()) {
            total += bin.getValue() * bin.getIdValue() * bin.getIdValue();
        }

        return Math.sqrt((total / getCount()) - (getMean() * getMean()));
    }
    
    public double getMedian() {
        double total = 0;
        double halfCount = getCount() / 2;
        for (Bin bin : values()) {
            total += bin.getValue();
            if (total >= halfCount) {
                return bin.getIdValue();
            }
        }
        return 0;
    }

    /** Gets the mode of the distribution (i.e. the largest bin). */
    public double getMode() {
        Bin mode = null;

        for (Bin bin : values()) {
            if (mode == null || mode.value < bin.value) {
                mode = bin;
            }
        }

        return mode.getIdValue();
    }

    public double getMin() {
        return firstEntry().getValue().getIdValue();
    }
    
    public double getMax() {
        return lastEntry().getValue().getIdValue();
    }
    
    public double getCount() {
        double count = 0;
        for (Bin bin : values()) {
            count += bin.value;
        }

        return count;
    }
}
