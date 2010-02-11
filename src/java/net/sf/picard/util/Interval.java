/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

import net.sf.picard.PicardException;
import net.sf.samtools.util.CoordMath;

import java.util.Collection;

/**
 * Represents a simple interval on a sequence.  Coordinates are 1-based closed ended.
 *
 * @author Tim Fennell
 */
public class Interval implements Comparable<Interval>, Cloneable {
    private final String sequence;
    private final int start;
    private final int end;
    private final boolean negativeStrand;
    private final String name;

    /**
     * Constructs an interval with the supplied sequence and start and end. If the end
     * position is less than the start position an exception is thrown.
     *
     * @param sequence the name of the sequence
     * @param start the start position of the interval on the sequence
     * @param end the end position of the interval on the sequence
     */
    public Interval(final String sequence, final int start, final int end) {
        this(sequence, start, end, false, null);

    }

    /**
     * Constructs an interval with the supplied sequence and start, end, strand and name.
     * If the end position is less than the start position an exception is thrown.
     *
     * @param sequence the name of the sequence
     * @param start the start position of the interval on the sequence
     * @param end the end position of the interval on the sequence
     * @param negative true to indicate negative strand, false otherwise
     * @param name the name (possibly null) of the interval
     *
     */
    public Interval(final String sequence, final int start, final int end, final boolean negative, final String name) {
        this.sequence = sequence;
        this.start = start;
        this.end = end;
        this.negativeStrand = negative;
        this.name = name;
        if (this.end < this.start) {
            throw new IllegalArgumentException("start must be less than or equal to end!");
        }
    }

    /** Gets the name of the sequence on which the interval resides. */
    public String getSequence() { return sequence; }

    /** Gets the 1-based start position of the interval on the sequence. */
    public int getStart() { return start; }

    /** Gets the 1-based closed-ended end position of the interval on the sequence. */
    public int getEnd() { return end; }

    /** Returns true if the interval is on the negative strand, otherwise false. */
    public boolean isNegativeStrand() { return this.negativeStrand; }

    /** Returns true if the interval is on the positive strand, otherwise false. */
    public boolean isPositiveStrand() { return !this.negativeStrand; }

    /** Returns the name of the interval, possibly null. */
    public String getName() { return this.name; }

    /** Returns true if this interval overlaps the other interval, otherwise false. */
    public boolean intersects(final Interval other) {
        return  (this.getSequence().equals(other.getSequence()) &&
                 CoordMath.overlaps(this.start, this.end, other.start, other.end));
    }

    public int getIntersectionLength(final Interval other) {
        if (this.intersects(other)) {
            return (int)CoordMath.getOverlap(this.getStart(), this.getEnd(), other.getStart(), other.getEnd());
        }
        return 0;
    }


    /** Returns a new Interval that represents the intersection between the two intervals. */
    public Interval intersect(final Interval that) {
        if (!intersects(that)) throw new IllegalArgumentException(that + " does not intersect " + this);
        return new Interval(this.sequence,
                            Math.max(this.start, that.start),
                            Math.min(this.end, that.end),
                            this.negativeStrand,
                            this.name + " intersection " + that.name);
    }


    /** Returns true if this interval overlaps the other interval, otherwise false. */
    public boolean abuts(final Interval other) {
        return this.getSequence().equals(other.getSequence()) &&
               (this.start == other.end + 1 || other.start == this.end + 1);
    }

    /** Gets the length of this interval. */
    public int length() { return this.end - this.start + 1; }

    /** Counts the total number of bases a collection of intervals. */
    public static long countBases(final Collection<Interval> intervals) {
        long total = 0;
        for (final Interval i : intervals) {
            total += i.length();
        }

        return total;
    }


    /**
     * Sort based on sequence.compareTo, then start pos, then end pos
     * with null objects coming lexically last
     */
    public int compareTo(final Interval that) {
        if (that == null) return -1; // nulls last

        int result = this.sequence.compareTo(that.sequence);
        if (result == 0) {
            if (this.start == that.start) {
                result = this.end - that.end;
            }
            else {
                result = this.start - that.start;
            }
        }

        return result;
    }

    /** Equals method that agrees with {@link #compareTo(Interval)}. */
    public boolean equals(final Object other) {
        if (!(other instanceof Interval)) return false;
        else if (this == other) return true;
        else {
            Interval that = (Interval)other;
            return (this.compareTo(that) == 0);
        }
    }

    @Override
    public int hashCode() {
        int result = sequence.hashCode();
        result = 31 * result + start;
        result = 31 * result + end;
        return result;
    }

    public String toString() {
        return getSequence() + ":" + start + "-" + end;
    }

    @Override
    public Interval clone() {
        try { return (Interval) super.clone(); }
        catch (CloneNotSupportedException cnse) { throw new PicardException("That's unpossible", cnse); }
    }
}
