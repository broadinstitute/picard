/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package picard.illumina.parser;

import picard.PicardException;

import java.util.List;
import java.util.Set;
import java.util.TreeMap;

/**
 * For per cycle files.  Maps a Cycle -> Tile -> List<File>
 *
 * @author jburke@broadinstitute.org
 */
class CycleIlluminaFileMap extends TreeMap<Integer, IlluminaFileMap> {
    /**
     * Return a CycleIlluminaFileMap with only the tiles listed and all of the cycles provided.
     * Important NOTE: this DOES NOT eliminate cycles from the cycles parameter passed in that are missing in the cyclesFileIterator of any given lane in the CycleIlluminaFileMap
     */
    public CycleIlluminaFileMap keep(final List<Integer> tilesToKeep, final Set<Integer> cycles) {
        final CycleIlluminaFileMap ciMap = new CycleIlluminaFileMap();

        if (cycles != null) {
            for (final int cycle : cycles) {
                final IlluminaFileMap template = this.get(cycle);
                if (template != null) {
                    ciMap.put(cycle, template.keep(tilesToKeep));
                }
            }
        }

        return ciMap;
    }

    /**
     * Assert that this map has an iterator for all of the expectedTiles and each iterator has expectedCycles number
     * of files.  Also, assert that each cycle file for a given tile is the same size
     *
     * @param expectedTiles  A list of tiles that should be in this map
     * @param expectedCycles The total number of files(cycles) that should be in each CycledFilesIterator
     */
    public void assertValid(final List<Integer> expectedTiles, final int[] expectedCycles) {
        if (size() != expectedCycles.length) {
            throw new PicardException("Expected CycledIlluminaFileMap to contain " + expectedCycles.length + " cycles but only " + size() + " were found!");
        }
        if (this.firstEntry().getValue().size() != expectedTiles.size()) {
            throw new PicardException("Expected CycledIlluminaFileMap to contain " + expectedTiles.size()
                    + " tiles but only " + this.firstEntry().getValue().size() + " were found!");
        }
    }

}
