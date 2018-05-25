/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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
package picard.sam.BamErrorMetric;

/**
 * Simple Pair class.
 *
 * @author Yossi Farjoun
 */

public class Pair<X extends Comparable<X>, Y extends Comparable<Y>> implements Comparable<Pair<X, Y>> {
    private X first;
    private Y second;

    public Pair(final X x, final Y y) {
        first = x;
        second = y;
    }

    public X getFirst() {
        return first;
    }

    public Y getSecond() {
        return second;
    }

    /**
     * Calculate whether this pair object is equal to another object.
     *
     * @param o The other object (hopefully a pair).
     * @return True if the two are equal; false otherwise.
     */
    @Override
    public boolean equals(final Object o) {
        // fast return in trivial case
        if (this == o) return true;

        if (o == null) return false;

        if (!(o instanceof Pair)) return false;

        final Pair other = (Pair) o;

        // Check to see whether one is null but not the other.
        if (this.first == null && other.first != null) return false;
        if (this.second == null && other.second != null) return false;

        // Check to see whether the values are equal.
        //  If the param of equals is null, it should by contract return false.
        if (this.first != null && !this.first.equals(other.first)) return false;
        if (this.second != null && !this.second.equals(other.second)) return false;

        return true;
    }

    /**
     * Basic hashcode function.  Assume hashcodes of first and second are
     * randomly distributed and return the XOR of the two.
     *
     * @return hashcode of the pair.
     */
    @Override
    public int hashCode() {
        if (second == null && first == null)
            return 0;
        if (second == null)
            return first.hashCode();
        if (first == null)
            return second.hashCode();
        return first.hashCode() ^ second.hashCode();
    }

    public String toString() {
        return first.toString() + "," + second.toString();
    }

    @Override
    public int compareTo(final Pair<X, Y> o) {
        if (this.first == null || o.first == null)
            throw new NullPointerException("Found null objects in comparing " + toString() + " and " + o.toString());

        final int firstCompare = this.first.compareTo(o.first);
        if (firstCompare != 0) return firstCompare;

        if (this.second == null || o.second == null)
            throw new NullPointerException("Found null objects in comparing " + toString() + " and " + o.toString());

        return this.second.compareTo(o.second);
    }
}
