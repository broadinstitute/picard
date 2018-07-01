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
package picard.sam.util;

import java.util.Objects;

/**
 * Simple Pair class.
 *
 * reimplemented since the commons one is final, and I wanted a different toString function...
 */

public class Pair<X extends Comparable<X>, Y extends Comparable<Y>> implements Comparable<Pair<X, Y>> {
    private X left;
    private Y right;

    public Pair(final X left, final Y right) {
        this.left = left;
        this.right = right;
    }

    public X getLeft() {
        return left;
    }

    public Y getRight() {
        return right;
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

        if (!(o instanceof Pair)) return false;

        final Pair other = (Pair) o;

        return Objects.equals(left,other.left) &&
               Objects.equals(right,other.right);
    }

    /**
     * Basic hashcode function.  Assume hashcodes of left and right are
     * randomly distributed and return the XOR of the two.
     *
     * @return hashcode of the pair.
     */
    @Override
    public int hashCode() {
        return Objects.hash(left, right);
    }

    public String toString() {
        return left + "," + right;
    }

    @Override
    public int compareTo(final Pair<X, Y> o) {
        final int leftCompare;
        if (this.left == null && o.left == null) {
            leftCompare = 0;
        } else {
            leftCompare = this.left.compareTo(o.left);
        }

        if (leftCompare != 0) return leftCompare;

        if (this.right == null && o.right == null) {
            return 0;
        }

        return this.right.compareTo(o.right);
    }
}
