/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
 * FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.util;

import picard.PicardException;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Contains class for figuring out the location of reads.
 *
 * @author Tim Fennell
 * @author Nils Homer
 * @author Yossi Farjoun
 */

/**
 * Small interface that provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  Tile should only allow
 * non-zero positive integers, x and y coordinates must be non-negative.
 * This is different from OpticalDuplicateFinder.PhysicalLocation in that the x and y positions are ints, not shorts
 * thus, they do not overflow within a HiSeqX tile.
 */
public class PhysicalLocation {
                                                        //FLOWCELL----:LANE-:TILE----:X_COORD-:Y_COORD-UNK
    public static final String DEFAULT_READ_NAME_REGEX = "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*";

    private final String readNameRegex;

    public PhysicalLocation() {this(DEFAULT_READ_NAME_REGEX);}

    public PhysicalLocation(final String readNameRegExp) {this.readNameRegex = readNameRegExp;}

    private Pattern readNamePattern;

    private short tile = -1;
    private int x = -1, y = -1;


    public short getTile() { return tile; }

    public void setTile(final short tile) { this.tile = tile; }

    public int getX() { return x; }

    public void setX(final int x) { this.x = x; }

    public int getY() { return y; }

    public void setY(final int y) { this.y = y; }


    private final int[] tmpLocationFields = new int[10]; // for optimization of addLocationInformation

    /**
     * Method used to extract tile/x/y from the read name and add it to the PhysicalLocation so that it
     * can be used later to determine optical duplication
     *
     * @param readName the name of the read/cluster
     * @param loc      the object to add tile/x/y to
     * @return true if the read name contained the information in parsable form, false otherwise
     */
    public boolean addLocationInformation(final String readName, final PhysicalLocation loc) {
        // Optimized version if using the default read name regex (== used on purpose):
        if (readNameRegex == DEFAULT_READ_NAME_REGEX) {
            final int fields = ReadNameParsingUtils.getRapidDefaultReadNameRegexSplit(readName, ':', tmpLocationFields);
            if (!(fields == 5 || fields == 7)) {
                throw new PicardException(String.format(" READ_NAME_REGEX '%s' did not match read name '%s'.  " ,
                        this.readNameRegex, readName));
            }

            final int offset = fields == 7 ? 2 : 0;
            loc.setTile((short) tmpLocationFields[offset + 2]);
            loc.setX(tmpLocationFields[offset + 3]);
            loc.setY(tmpLocationFields[offset + 4]);
            return true;
        } else if (readNameRegex == null) {
            return false;
        } else {
            // Standard version that will use the regex
            if (readNamePattern == null) readNamePattern = Pattern.compile(readNameRegex);

            final Matcher m = readNamePattern.matcher(readName);
            if (m.matches()) {
                loc.setTile((short) Integer.parseInt(m.group(1)));
                loc.setX(Integer.parseInt(m.group(2)));
                loc.setY(Integer.parseInt(m.group(3)));
                return true;
            } else {
                throw new PicardException(String.format("READ_NAME_REGEX '%s' did not match read name '%s'. ", readNameRegex, readName));
            }
        }
    }
}
