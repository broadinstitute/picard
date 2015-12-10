package picard.sam.util;

import htsjdk.samtools.util.Log;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Small class that provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile should only allow
 * non-zero positive integers, x and y coordinates may be negative.
 */
public class ReadNameParser {

    public static final String DEFAULT_READ_NAME_REGEX = "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*".intern();

    public static final int DEFAULT_OPTICAL_DUPLICATE_DISTANCE = 100;

    private final int[] tmpLocationFields = new int[3]; // for optimization of addLocationInformation

    public String readNameRegex;

    public int opticalDuplicatePixelDistance;

    private Pattern readNamePattern;

    private boolean warnedAboutRegexNotMatching = false;

    private final Log log;

    public ReadNameParser(final Log log) {
        this.log = log;
    }

    /**
     * Method used to extract tile/x/y from the read name and add it to the PhysicalLocationShort so that it
     * can be used later to determine optical duplication
     *
     * @param readName the name of the read/cluster
     * @param loc the object to add tile/x/y to
     * @return true if the read name contained the information in parsable form, false otherwise
     */
    public boolean addLocationInformation(final String readName, final PhysicalLocation loc) {
        // Optimized version if using the default read name regex (== used on purpose):
        if (this.readNameRegex == this.DEFAULT_READ_NAME_REGEX) {
            final int fields = getLastThreeFields(readName, ':', tmpLocationFields);
            if (!(fields == 5 || fields == 7)) {
                if (null != log && !this.warnedAboutRegexNotMatching) {
                    this.log.warn(String.format("Default READ_NAME_REGEX '%s' did not match read name '%s'.  " +
                                    "You may need to specify a READ_NAME_REGEX in order to correctly identify optical duplicates.  " +
                                    "Note that this message will not be emitted again even if other read names do not match the regex.",
                            this.readNameRegex, readName));
                    this.warnedAboutRegexNotMatching = true;
                }
                return false;
            }
            loc.setTile((short) tmpLocationFields[0]);
            loc.setX(tmpLocationFields[1]);
            loc.setY(tmpLocationFields[2]);
            return true;
        } else if (this.readNameRegex == null) {
            return false;
        } else {
            // Standard version that will use the regex
            if (this.readNamePattern == null) this.readNamePattern = Pattern.compile(this.readNameRegex);

            final Matcher m = this.readNamePattern.matcher(readName);
            if (m.matches()) {
                loc.setTile((short) Integer.parseInt(m.group(1)));
                loc.setX(Integer.parseInt(m.group(2)));
                loc.setY(Integer.parseInt(m.group(3)));
                return true;
            } else {
                if (null != log && !this.warnedAboutRegexNotMatching) {
                    this.log.warn(String.format("READ_NAME_REGEX '%s' did not match read name '%s'.  Your regex may not be correct.  " +
                                    "Note that this message will not be emitted again even if other read names do not match the regex.",
                            this.readNameRegex, readName));
                    warnedAboutRegexNotMatching = true;
                }
                return false;
            }
        }
    }

    /**
     * Single pass method to parse the read name for the default regex.  Examines the last three fields as split by the delimiter.
     */
    public static int getLastThreeFields(final String readName, final char delim, final int[] tokens) {
        int tokensIdx = 2; // start at the last token
        int numFields = 0;
        int i, endIdx;
        endIdx = readName.length();
        // find the last three tokens only
        for (i = readName.length() - 1; 0 <= i && 0 <= tokensIdx; i--) {
            if (readName.charAt(i) == delim || 0 == i) {
                numFields++;
                tokens[tokensIdx] = rapidParseInt(readName.substring(i+1, endIdx));
                tokensIdx--;
                endIdx = i;
            }
        }
        // continue to find the # of fields
        while (0 <= i) {
            if (readName.charAt(i) == delim || 0 == i) numFields++;
            i--;
        }
        if (numFields < 3) {
            tokens[0] = tokens[1] = tokens[2] = -1;
            return -1;
        }
        else {
            return numFields;
        }
    }

    /**
     * Very specialized method to rapidly parse a sequence of digits from a String up until the first
     * non-digit character.
     */
    public static int rapidParseInt(final String input) {
        final int len = input.length();
        int val = 0;
        int i = 0;
        boolean isNegative = false;

        if (0 < len && '-' == input.charAt(0)) {
            i = 1;
            isNegative = true;
        }

        for (; i < len; ++i) {
            final char ch = input.charAt(i);
            if (Character.isDigit(ch)) {
                val = (val * 10) + (ch - 48);
            } else {
                break;
            }
        }

        if (isNegative) val = -val;

        return val;
    }
}
