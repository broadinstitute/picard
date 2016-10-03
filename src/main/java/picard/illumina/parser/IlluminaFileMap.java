/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

/**
 * For "non-cycle" files (files that have multiple cycles per file).  Maps a Tile -> File
 * @author jburke@broadinstitute.org
 */
class IlluminaFileMap extends TreeMap<Integer, File> {

    /** Return a file map that includes only the tiles listed */
    public IlluminaFileMap keep(final List<Integer> tilesToKeep) {
        final IlluminaFileMap fileMap = new IlluminaFileMap();
        for(final Integer tile : tilesToKeep) {
            final File file = this.get(tile);
            if(file != null) {
                fileMap.put(tile, file);
            }
        }
        return fileMap;
    }

    /**
     * Return the List of Files in order starting at the given tile and containing all files with tile numbers greater than startingTile that
     * are within this map
     * @param startingTile The first File in the returned list will correspond to this tile
     * @return A List of files for all tiles >= startingTile that are contained in this FileMap
     */
    public List<File> getFilesStartingAt(final int startingTile) {
        return new ArrayList<File>(this.tailMap(startingTile).values());
    }
}
