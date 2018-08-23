package picard.util;

/**
 * Created by farjoun on 8/22/18.
 */
public class SequenceUtil {
    /**
     * Gives the edit distance between this barcode and another of the same length.
     */
    public static byte calculateEditDistance(final String lhs, final String rhs) {
        if (lhs.length() != rhs.length()) {
            throw new IllegalArgumentException(String.format("Inputs must be of the same length: '%s'.length()!='%s'.length()", lhs, rhs));
        }
        byte tmp = 0;
        for (int i = 0; i < rhs.length(); ++i) {
            if (rhs.charAt(i) != lhs.charAt(i)) ++tmp;
        }
        return tmp;
    }
}

