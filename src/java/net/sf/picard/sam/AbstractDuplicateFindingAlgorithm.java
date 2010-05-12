package net.sf.picard.sam;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Abstract class that holds parameters and methods common to classes that perform duplicate
 * detection and/or marking within SAM/BAM files.
 *
 * @author Tim Fennell
 */
public abstract class AbstractDuplicateFindingAlgorithm extends CommandLineProgram {
    @Option(doc="Regular expression that can be used to parse read names in the incoming SAM file. Read names are " +
            "parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used " +
            "to estimate the rate of optical duplication in order to give a more accurate estimated library size. " +
            "The regular expression should contain three capture groups for the three variables, in order.")
    public String READ_NAME_REGEX = "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*";
    
    @Option(doc="The maximum offset between two duplicte clusters in order to consider them optical duplicates. This " +
            "should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the " +
            "Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.")
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = 100;

    private Pattern READ_NAME_PATTERN = null;

    /**
     * Small interface that provides access to the physical location information about a cluster.
     * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile should only allow
     * non-zero positive integers, x and y coordinates may be negative.
     */
    public static interface PhysicalLocation {
        short getReadGroup();
        void  setReadGroup(short rg);
        byte  getTile();
        void  setTile(byte tile);
        short getX();
        void  setX(short x);
        short getY();
        void  setY(short y);
    }
    
    /**
     * Method used to extract tile/x/y from the read name and add it to the PhysicalLocation so that it
     * can be used later to determine optical duplication
     *
     * @param readName the name of the read/cluster
     * @param loc the object to add tile/x/y to
     * @return true if the read name contained the information in parsable form, false otherwise
     */
    boolean addLocationInformation(final String readName, final PhysicalLocation loc) {
        if (READ_NAME_PATTERN == null) {
            if (READ_NAME_REGEX == null) return false;
            READ_NAME_PATTERN = Pattern.compile(READ_NAME_REGEX);
        }
        final Matcher m = READ_NAME_PATTERN.matcher(readName);
        if (m.matches()) {
            loc.setTile((byte) Integer.parseInt(m.group(1)));
            loc.setX((short) Integer.parseInt(m.group(2)));
            loc.setY((short) Integer.parseInt(m.group(3)));
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Finds which reads within the list of duplicates are likely to be optical duplicates of
     * one another.
     *
     * Note: this method will perform a sort() of the list; if it is imperative that the list be
     * unmodified a copy of the list should be passed to this method.
     *
     * @param list a list of reads that are determined to be duplicates of one another
     * @param maxDistance maximum distance in x and y directions for reads to be considered optical duplicates
     * @return a boolean[] of the same length as the incoming list marking which reads are optical duplicates
     */
    boolean[] findOpticalDuplicates(final List<? extends PhysicalLocation> list, final int maxDistance) {
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

        outer: for (int i=0; i<length; ++i) {
            PhysicalLocation lhs = list.get(i);
            if (lhs.getTile() < 0) continue;

            for (int j=i+1; j<length; ++j) {
                PhysicalLocation rhs = list.get(j);

                if (lhs.getReadGroup() != rhs.getReadGroup()) continue outer;
                if (lhs.getTile() != rhs.getTile()) continue outer;
                if (rhs.getX() > lhs.getX() + maxDistance) continue outer;

                if (Math.abs(lhs.getY()  - rhs.getY()) <= maxDistance) {
                    opticalDuplicateFlags[j] = true;
                }
            }
        }
        return opticalDuplicateFlags;
    }
}
