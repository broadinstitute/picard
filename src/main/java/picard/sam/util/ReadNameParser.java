package picard.sam.util;

import htsjdk.samtools.util.Log;

import java.io.Serializable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile should only allow
 * non-zero positive integers, x and y coordinates may be negative.
 */
public class ReadNameParser implements Serializable {

    /**
     * The read name regular expression (regex) is used to extract three pieces of information from the read name: tile, x location,
     * and y location.  Any read name regex should parse the read name to produce these and only these values.  An example regex is:
     *  (?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$
     * which assumes that fields in the read name are delimited by ':' and the last three fields correspond to the tile, x and y locations,
     * ignoring any trailing non-digit characters.
     *
     * The default regex is optimized for fast parsing (see {@link #getLastThreeFields(String, char, int[])}) by searching for the last
     * three fields, ignoring any trailing non-digit characters, assuming the delimiter ':'.  This should consider correctly read names
     * where we have 5 or 7 field with the last three fields being tile/x/y, as is the case for the majority of read names produced by
     * Illumina technology.
     */
    public static final String DEFAULT_READ_NAME_REGEX = "<optimized capture of last three ':' separated fields as numeric values>".intern();

    /**
     * Opt-in variant of {@link #DEFAULT_READ_NAME_REGEX} for read names that carry a trailing UMI field, e.g.
     * {@code FC:1:1101:12345:67890:ACGT-TTGC}. The three numeric fields before the optional UMI field are the tile, x and
     * y values. The UMI field may contain only A, C, G, T, N, + and -. Parsing is done without a regex and gives the same
     * results as the regex {@code (?:.*:)?([0-9]+):([0-9]+):([0-9]+)(:[ACGTN+-]+)?$}.
     */
    public static final String ILLUMINA_UMI_READ_NAME_REGEX =
            "<optimized capture of last three ':' separated numeric fields, with an optional trailing UMI field>".intern();

    private String readNameStored = null;

    private PhysicalLocation physicalLocationStored = new PhysicalLocationInt();

    private final int[] tmpLocationFields = new int[3]; // for optimization of addLocationInformation

    private final boolean useOptimizedDefaultParsing; // was the regex default?

    private final boolean useOptimizedUmiParsing; // was the regex ILLUMINA_UMI_READ_NAME_REGEX?

    protected final String readNameRegex;

    private Pattern readNamePattern;

    private boolean warnedAboutRegexNotMatching = false;

    private transient final Log log;

    /**
     * Creates are read name parser using the default read name regex and optical duplicate distance.   See {@link #DEFAULT_READ_NAME_REGEX}
     * for an explanation on how the read name is parsed.
     */
    public ReadNameParser() {
        this(DEFAULT_READ_NAME_REGEX);
    }

    /**
     * Creates are read name parser using the given read name regex.  See {@link #DEFAULT_READ_NAME_REGEX} for an explanation on how to
     * format the regular expression (regex) string.
     * @param readNameRegex the read name regular expression string to parse read names, null to never parse location information.
     */
    public ReadNameParser(final String readNameRegex) {
        this(readNameRegex, null);
    }

    /**
     * Creates are read name parser using the given read name regex.  See {@link #DEFAULT_READ_NAME_REGEX} for an explanation on how to
     * format the regular expression (regex) string.
     * @param readNameRegex the read name regular expression string to parse read names, null to never parse location information..
     * @param log the log to which to write messages.
     */
    public ReadNameParser(final String readNameRegex, final Log log) {
        this.useOptimizedDefaultParsing = DEFAULT_READ_NAME_REGEX.equals(readNameRegex);
        this.useOptimizedUmiParsing = ILLUMINA_UMI_READ_NAME_REGEX.equals(readNameRegex);
        this.readNameRegex = readNameRegex;
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
    private boolean readLocationInformation(final String readName, final PhysicalLocation loc) {
        try {
            // Optimized version if using the default read name regex (== used on purpose):
            if (this.useOptimizedDefaultParsing) {
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
            } else if (this.useOptimizedUmiParsing) {
                if (!parseLastThreeNumericFieldsWithOptionalUmi(readName, tmpLocationFields)) {
                    if (null != log && !this.warnedAboutRegexNotMatching) {
                        this.log.warn(String.format("READ_NAME_FORMAT ILLUMINA_UMI did not match read name '%s'.  " +
                                        "Note that this message will not be emitted again even if other read names do not match.",
                                readName));
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
        catch (NumberFormatException nfe) {
            if (log != null && !this.warnedAboutRegexNotMatching) {
                this.log.warn("A field field parsed out of a read name was expected to contain an integer and did not. ",
                              "Read name: ", readName, ". Cause: ", nfe.getMessage());
                warnedAboutRegexNotMatching = true;
            }
            return false;
        }
    }

    public boolean addLocationInformation(final String readName, final PhysicalLocation loc){
        if (!readName.equals(readNameStored)) {
            if (readLocationInformation(readName, loc)) {
                readNameStored = readName;
                physicalLocationStored.setX(loc.getX());
                physicalLocationStored.setY(loc.getY());
                physicalLocationStored.setTile(loc.getTile());
                return true;
            }
            // return false if read name cannot be parsed
            return false;
        } else {
            loc.setTile(physicalLocationStored.getTile());
            loc.setX(physicalLocationStored.getX());
            loc.setY(physicalLocationStored.getY());
            return true;
        }
    }

    /**
     * Given a string, splits the string by the delimiter, and returns the the last three fields parsed as integers.  Parsing a field
     * considers only a sequence of digits up until the first non-digit character.  The three values are stored in the passed-in array.
     *
     * @throws NumberFormatException if any of the tokens that should contain numbers do not start with parsable numbers
     */
    public static int getLastThreeFields(final String readName, final char delim, final int[] tokens) throws NumberFormatException {
        int tokensIdx = 2; // start at the last token
        int numFields = 0;
        int i, endIdx;
        endIdx = readName.length();
        // find the last three tokens only
        for (i = readName.length() - 1; 0 <= i && 0 <= tokensIdx; i--) {
            if (readName.charAt(i) == delim || 0 == i) {
                numFields++;
                tokens[tokensIdx] = rapidParseInt(readName.substring((0 == i) ? 0 : (i+1), endIdx));
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
     * Parses {@code ...:tile:x:y} or {@code ...:tile:x:y:UMI} (the UMI optional, A/C/G/T/N/+/- only), storing tile, x and y in
     * {@code tokens}. Matches the regex {@code (?:.*:)?([0-9]+):([0-9]+):([0-9]+)(:[ACGTN+-]+)?$}.
     *
     * @return true if the name matched, false otherwise
     */
    private static boolean parseLastThreeNumericFieldsWithOptionalUmi(final String readName, final int[] tokens) {
        // Strip a trailing UMI field, if present. Its characters never include digits, so the two cases cannot overlap.
        final int lastColon = readName.lastIndexOf(':');
        final String head = (lastColon >= 0 && isUmi(readName, lastColon + 1))
                ? readName.substring(0, lastColon)
                : readName;

        // The last three fields of head are tile, x and y; their delimiters are the last colon (c2) and the two before it.
        final int c2 = head.lastIndexOf(':');
        if (c2 < 0) {
            return false;
        }
        final int c1 = head.lastIndexOf(':', c2 - 1);
        if (c1 < 0) {
            return false;
        }
        final int c0 = head.lastIndexOf(':', c1 - 1); // -1 when the tile field starts at the beginning of the name

        final String tileField = head.substring(c0 + 1, c1);
        final String xField = head.substring(c1 + 1, c2);
        final String yField = head.substring(c2 + 1);
        if (!isDigits(tileField) || !isDigits(xField) || !isDigits(yField)) {
            return false;
        }

        try {
            final int tile = Integer.parseInt(tileField);
            final int x = Integer.parseInt(xField);
            final int y = Integer.parseInt(yField);
            tokens[0] = tile;
            tokens[1] = x;
            tokens[2] = y;
            return true;
        } catch (final NumberFormatException nfe) {
            return false; // out of int range, as the regex path rejects it
        }
    }

    /** True if {@code readName} has at least one character from {@code start} on, and all of them are A, C, G, T, N, + or -. */
    private static boolean isUmi(final String readName, final int start) {
        if (start >= readName.length()) {
            return false;
        }
        for (int i = start; i < readName.length(); i++) {
            if ("ACGTN+-".indexOf(readName.charAt(i)) < 0) {
                return false;
            }
        }
        return true;
    }

    /** True if {@code s} is non-empty and consists only of the ASCII digits 0-9. */
    private static boolean isDigits(final String s) {
        if (s.isEmpty()) {
            return false;
        }
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) < '0' || s.charAt(i) > '9') {
                return false;
            }
        }
        return true;
    }

    /**
     * Very specialized method to rapidly parse a sequence of digits from a String up until the first
     * non-digit character.
     *
     * @throws NumberFormatException if the String does not start with an optional - followed by at least on digit
     */
    public static int rapidParseInt(final String input) throws NumberFormatException {
        final int len = input.length();
        int val = 0;
        int i = 0;
        boolean isNegative = false;

        if (0 < len && '-' == input.charAt(0)) {
            i = 1;
            isNegative = true;
        }

        boolean hasDigits = false;
        for (; i < len; ++i) {
            final char ch = input.charAt(i);
            if (Character.isDigit(ch)) {
                val = (val * 10) + (ch - 48);
                hasDigits = true;
            } else {
                break;
            }
        }

        if (!hasDigits) throw new NumberFormatException("String '" + input + "' did not start with a parsable number.");
        if (isNegative) val = -val;
        return val;
    }

}
