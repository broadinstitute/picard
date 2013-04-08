/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

package org.broad.tribble.index.interval;

import org.broad.tribble.index.Block;


/**
 *  Quick and dirty interval class
 *  Describes a genomic interval and where in a file information for that
 *  interval can be obtained
 */
public class Interval implements Comparable {
    /**
     * Start of the interval in genomic coordinates -- this is exposed on purpose, getters have a significant
     * performance penalty for this field.
     */
     final int start;

    /**
     * End of the interval in genomic coordinates -- this is exposed on purpose, getters have a significant
     * performance penalty for this field.
     */
     final int end;

    /**
     * File block  (position, size) containing the data for this interval
     */
    private Block block;

    public Interval(int start, int end) {
        assert start <= end;
        this.start = start;
        this.end = end;
    }


    public Interval(int start, int end, Block block) {
        assert start <= end;
        this.start = start;
        this.end = end;
        this.block = block;
    }


    public boolean equals(Object other) {
        if (this == other)
            return true;
        if (this.getClass().equals(other.getClass())) {
            Interval otherInterval = (Interval) other;
            return (this.start == otherInterval.start &&
                    this.end == otherInterval.end);
        }
        return false;
    }


    public int hashCode() {
        return start;
    }


    public int compareTo(Object o) {
        Interval other = (Interval) o;
        if (this.start < other.start)
            return -1;
        if (this.start > other.start)
            return 1;

        if (this.end < other.end)
            return -1;
        if (this.end > other.end)
            return 1;

        return 0;
    }

    public String toString() {
        return "Interval[" + this.start + ", " + this.end + "]";
    }


    /**
     * @return whether this interval overlaps the other.
     */
    public boolean overlaps(Interval other) {
        return (this.start <= other.end &&
                other.start <= this.end);
    }


    /**
     * @return The file block for this interval
     */
    public Block getBlock() {
        return block;
    }
}

