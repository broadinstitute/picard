package picard.sam;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;

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
    private static Log LOG = Log.getInstance(AbstractDuplicateFindingAlgorithm.class);

    private static final String DEFAULT_READ_NAME_REGEX = "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*".intern();

    @Option(doc="Regular expression that can be used to parse read names in the incoming SAM file. Read names are " +
            "parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used " +
            "to estimate the rate of optical duplication in order to give a more accurate estimated library size. " +
            "The regular expression should contain three capture groups for the three variables, in order. " +
            "It must match the entire read name. " +
            "Note that if the default regex is specified, a regex match is not actually done, but instead the read name " +
            " is split on colon character. " +
            "For 5 element names, the 2nd, 3rd and 4th elements are assumed to be tile, x and y values. " +
            "For 7 element names (CASAVA 1.8), the 4th, 5th, and 6th elements are assumed to be tile, x and y values.")
    public String READ_NAME_REGEX = DEFAULT_READ_NAME_REGEX;
    
    @Option(doc="The maximum offset between two duplicte clusters in order to consider them optical duplicates. This " +
            "should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the " +
            "Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.")
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = 100;

    private Pattern READ_NAME_PATTERN = null;

    private boolean warnedAboutRegexNotMatching = false;

    /**
     * Small interface that provides access to the physical location information about a cluster.
     * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile should only allow
     * non-zero positive integers, x and y coordinates may be negative.
     */
    public static interface PhysicalLocation {
        short getReadGroup();
        void  setReadGroup(short rg);
        short  getTile();
        void  setTile(short tile);
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
    private final String[] tmpLocationFields = new String[10];
    boolean addLocationInformation(final String readName, final PhysicalLocation loc) {
        // Optimized version if using the default read name regex (== used on purpose):
        if (READ_NAME_REGEX == DEFAULT_READ_NAME_REGEX) {
            final int fields = StringUtil.split(readName, tmpLocationFields, ':');

            if (!(fields == 5 || fields == 7)) {
                if (!warnedAboutRegexNotMatching) {
                    LOG.warn(String.format("Default READ_NAME_REGEX '%s' did not match read name '%s'.  " +
                            "You may need to specify a READ_NAME_REGEX in order to correctly identify optical duplicates.  " +
                            "Note that this message will not be emitted again even if other read names do not match the regex.",
                            READ_NAME_REGEX, readName));
                    warnedAboutRegexNotMatching = true;
                }
                return false;
            }

            final int offset = fields == 7 ? 2 : 0;

            loc.setTile((short) rapidParseInt(tmpLocationFields[offset + 2]));
            loc.setX((short) rapidParseInt(tmpLocationFields[offset + 3]));
            loc.setY((short) rapidParseInt(tmpLocationFields[offset + 4]));
            return true;
        }
        else if (READ_NAME_REGEX == null) {
            return false;
        }
        else {
            // Standard version that will use the regex
            if (READ_NAME_PATTERN == null) READ_NAME_PATTERN = Pattern.compile(READ_NAME_REGEX);

            final Matcher m = READ_NAME_PATTERN.matcher(readName);
            if (m.matches()) {
                loc.setTile((short) Integer.parseInt(m.group(1)));
                loc.setX((short) Integer.parseInt(m.group(2)));
                loc.setY((short) Integer.parseInt(m.group(3)));
                return true;
            }
            else {
                if (!warnedAboutRegexNotMatching) {
                    LOG.warn(String.format("READ_NAME_REGEX '%s' did not match read name '%s'.  Your regex may not be correct.  " +
                            "Note that this message will not be emitted again even if other read names do not match the regex.",
                            READ_NAME_REGEX, readName));
                    warnedAboutRegexNotMatching = true;
                }
                return false;
            }
        }
    }

    /**
     * Very specialized method to rapidly parse a sequence of digits from a String up until the first
     * non-digit character. Does not handle negative numbers.
     */
    private final int rapidParseInt(final String input) {
        final int len = input.length();
        int val = 0;

        for (int i=0; i<len; ++i) {
            final char ch = input.charAt(i);
            if (Character.isDigit(ch)) {
                val = (val*10) + (ch-48);
            }
        }

        return val;
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
