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
package net.sf.picard.illumina.parser.readers;

import net.sf.picard.PicardException;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *  The position files of Illumina are nearly the same form:  Pos files consist of text based tabbed
 *  x-y coordinate float pairs, locs files are binary x-y float pairs, clocs are compressed binary
 *  x-y float pairs.  Each of these file types we read sequentially and are really concerned with
 *  iterating over the coordinates and returning them as as they would appear in a QSeq file.
 *  Therefore, this abstract base class provides the basic functionality for iterating over
 *  the values found in these files and converting them into qseq style coordinates.
 *
 *  Currently these readers also return lane/tile but this will be unnecessary in future releases.
 */

public abstract class AbstractIlluminaPositionFileReader implements CloseableIterator<AbstractIlluminaPositionFileReader.PositionInfo> {
    public static final float MAX_POS = 9999999.99f;

    public class PositionInfo {
        /** The x-position as it occurs in the file being read */
        public final float xPos;

        /** The y-position as it occurs in the file being read*/
        public final float yPos;

        /** The lane, which is determined from the file name*/
        public final int lane;

        /** The tile, which is determined from the file name*/
        public final int tile;

        /** The QSeq style x-coordinat, an integer = Math.round(xPos*10 + 1000)*/
        public final int xQseqCoord;

        /** The QSeq style y-coordinates, an integer = Math.round(yPos*10 + 1000)*/
        public final int yQseqCoord;

        public PositionInfo(final float x, final float y, final int lane, final int tile) {
            if(x < 0 || y < 0) {
                throw new IllegalArgumentException("X and Y position values must be positive (x,y)=(" + x + ", y) lane(" + lane + ") tile(" + tile + ")");
            }

            if(x > MAX_POS || y > MAX_POS) {
                throw new IllegalArgumentException("X and Y position values must be less than " + MAX_POS + " (x,y)=(" + x + ", y) lane(" + lane + ") tile(" + tile + ")");
            }

            this.xPos = x;
            this.yPos = y;
            this.xQseqCoord = posToQSeqCoord(x);
            this.yQseqCoord = posToQSeqCoord(y);
            this.lane = lane;
            this.tile = tile;
        }

        /** Convert a value in float form as it occurs in pos,locs,and clocs files into integer as it is found in QSeqs */
        private int posToQSeqCoord(final float pos) {
            return Math.round(pos * 10 + 1000);
        }

        public boolean equals(final Object other) {
            if(other == null || other.getClass() != AbstractIlluminaPositionFileReader.PositionInfo.class) {
                return false;
            }
            if(other == this) return true;
            final PositionInfo otherPi = (PositionInfo) other;
            return this.xPos == otherPi.xPos  && this.yPos == otherPi.yPos &&
                   this.lane == otherPi.lane  && this.tile == otherPi.tile &&
                   this.xQseqCoord == otherPi.xQseqCoord && this.yQseqCoord == otherPi.yQseqCoord;
        }
    }

    //Note: Perhaps use the IlluminaFileUtil to do this part
    private static Pattern FileNamePattern = Pattern.compile("^s_(\\d+)_(\\d+)(_pos\\.txt|\\.locs|\\.clocs|_pos\\.txt.gz|_pos\\.txt.bz2)$");

    private final File file;
    private final int lane;
    private final int tile;

    public AbstractIlluminaPositionFileReader(final File file) {
        this.file = file;

        int [] laneAndTile = fileNameToLaneAndTile(file.getName());
        lane = laneAndTile[0];
        tile = laneAndTile[1];
    }

    public int getTile() {
        return tile;
    }

    public int getLane() {
        return lane;
    }

    public File getFile() {
        return file;
    }

    /** Extract the lane/tile from the given filename **/
    private int [] fileNameToLaneAndTile(final String fileName) {
        final String [] tokens = fileName.split(File.pathSeparator);
        final Matcher matcher = FileNamePattern.matcher(tokens[tokens.length - 1]);
        if(!matcher.matches()) {
            throw new PicardException("File name not of the right structure: <filePath>/s_<lane>_<tile>(_pos.txt|_pos.txt.gz|_pos.txt.bz2.locs|.clocs).  File name (" + fileName + ")");
        }

        return new int[]{Integer.parseInt(matcher.group(1)), Integer.parseInt(matcher.group(2))};
    }

    /** Return the next set of coordinates in a given file. **/
    public final PositionInfo next() {
        if(! hasNext()) {
            throw new NoSuchElementException("No such cluster, cluster count(" + makeExceptionMsg() +")");
        }
        return unsafeNextInfo();
    }

    /** Returns the next position info.  Implementations of this method do not need to call hasNext since
     * it is called in next() */
    protected abstract PositionInfo unsafeNextInfo();

    /** Create a string that will be included in any NoSuchElementException thrown by the next() method */
    protected abstract String makeExceptionMsg();

    /** Return true if the file has more elements to return, false otherwise */
    public abstract boolean hasNext();

    public void remove() {
        throw new UnsupportedOperationException();
    }
}
