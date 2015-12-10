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

package picard.sam.markduplicates.util;

import htsjdk.samtools.util.Log;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.ReadNameParser;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Contains methods for finding optical duplicates.
 *
 * @author Tim Fennell
 * @author Nils Homer
 */
public class OpticalDuplicateFinder extends ReadNameParser {

    public OpticalDuplicateFinder() {
        this(DEFAULT_READ_NAME_REGEX, DEFAULT_OPTICAL_DUPLICATE_DISTANCE);
    }

    public OpticalDuplicateFinder(final int opticalDuplicatePixelDistance) {
        this(DEFAULT_READ_NAME_REGEX, opticalDuplicatePixelDistance);
    }

    public OpticalDuplicateFinder(final String readNameRegex) {
        this(readNameRegex, DEFAULT_OPTICAL_DUPLICATE_DISTANCE);
    }

    public OpticalDuplicateFinder(final String readNameRegex, final int opticalDuplicatePixelDistance) {
        this(readNameRegex, opticalDuplicatePixelDistance, null);
    }

    public OpticalDuplicateFinder(final String readNameRegex, final int opticalDuplicatePixelDistance, final Log log) {
        super(log);
        this.readNameRegex = readNameRegex;
        this.opticalDuplicatePixelDistance = opticalDuplicatePixelDistance;
    }

    /**
     * Finds which reads within the list of duplicates are likely to be optical duplicates of
     * one another.
     * <p/>
     * Note: this method will perform a sort() of the list; if it is imperative that the list be
     * unmodified a copy of the list should be passed to this method.
     *
     * @param list a list of reads that are determined to be duplicates of one another
     * @return a boolean[] of the same length as the incoming list marking which reads are optical duplicates
     */
    public boolean[] findOpticalDuplicates(final List<? extends PhysicalLocation> list) {
        final int length = list.size();
        final boolean[] opticalDuplicateFlags = new boolean[length];

        Collections.sort(list, new Comparator<PhysicalLocation>() {
            public int compare(final PhysicalLocation lhs, final PhysicalLocation rhs) {
                int retval = lhs.getReadGroup() - rhs.getReadGroup();
                if (retval == 0) retval = lhs.getTile() - rhs.getTile();
                if (retval == 0) retval = lhs.getX() - rhs.getX();
                if (retval == 0) retval = lhs.getY() - rhs.getY();
                return retval;
            }
        });

        outer:
        for (int i = 0; i < length; ++i) {
            final PhysicalLocation lhs = list.get(i);
            if (lhs.getTile() < 0) continue;

            for (int j = i + 1; j < length; ++j) {
                final PhysicalLocation rhs = list.get(j);

                if (opticalDuplicateFlags[j]) continue;
                if (lhs.getReadGroup() != rhs.getReadGroup()) continue outer;
                if (lhs.getTile() != rhs.getTile()) continue outer;
                if (rhs.getX() > lhs.getX() + this.opticalDuplicatePixelDistance) continue outer;

                if (Math.abs(lhs.getY() - rhs.getY()) <= this.opticalDuplicatePixelDistance) {
                    opticalDuplicateFlags[j] = true;
                }
            }
        }
        return opticalDuplicateFlags;
    }
}
