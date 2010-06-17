package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;

/**
 * Misc utilities for working with Illuina specific files and data
 *
 * @author alecw@broadinstitute.org
 */
public class IlluminaUtil {
    /**
     * Standard Broad algorithm for creating a read name
     * @return read name in standard Broad format, without end suffix.
     */
    public static String makeReadName(final String runBarcode, final int lane, final int tile, final int xCoordinate,
                                      final int yCoordinate) {
        return runBarcode + ":" + lane + ":" + tile + ":" + xCoordinate + ":" + yCoordinate;
    }

    /**
     * Parse the tile # from the read name.
     * If we find that there are other elements needed from the read name, it might be a good idea to put
     * makeReadName() and various get..() methods into a new class.
     *
     * @param readName As produced by IlluminaUtil.makeReadName()
     * @return tile number, or null if read name is not in correct format.
     */
    public static Integer getTileFromReadName(final String readName) {
        final int first = readName.indexOf(':');
        if (first > 0) {
            final int second = readName.indexOf(':', first+1);
            if (second > 0) {
                final int third = readName.indexOf(':', second+1);
                if (third > 0) {
                    return Integer.parseInt(readName.substring(second+1, third));
                }
            }
        }

        return null;
    }

    /**
     * Convert from Solexa-scaled ASCII qualities to Phred-scaled binary.  The only difference is Solexa qualities have
     * 64 added to the phred binary to make them printable.
     *
     * @param solexaQualities Printable ASCII qualities.
     * @return binary Phred-scaled qualities.
     */
    public static byte[] makePhredBinaryFromSolexaQualityAscii_1_3(final String solexaQualities) {
        return makePhredBinaryFromSolexaQualityAscii_1_3(solexaQualities, 0, solexaQualities.length());
    }

    /**
     * Convert from Solexa-scaled ASCII qualities to Phred-scaled binary.  The only difference is Solexa qualities have
     * 64 added to the phred binary to make them printable.
     *
     * @param solexaQualities Printable ASCII qualities.
     * @param offset Character at which to start conversion.
     * @param length Number of characters to convert.
     * @return binary Phred-scaled qualities.
     */
    public static byte[] makePhredBinaryFromSolexaQualityAscii_1_3(final String solexaQualities, final int offset, final int length) {
        final byte[] quals = StringUtil.stringToBytes(solexaQualities, offset, length);
        SolexaQualityConverter.getSingleton().convertSolexa_1_3_QualityCharsToPhredBinary(quals);
        return quals;
    }

    /**
     * Converts from Solexa ASCII to Phred binary in place.  These are the older-style qualities
     * rather than Phred qualities with a different addend to make them printable.
     */
    public static void convertSolexaQualityAscii_1_1_ToPhredBinary(final byte[] solexaQualities) {
        SolexaQualityConverter.getSingleton().convertSolexaQualityCharsToPhredBinary(solexaQualities);
    }

    /**
     * Get a Solexa ASCII quality value from an array of strings that are integer qualities in this order:
     * [cycle-1-A, cycle-1-C, cycle-1-G, cycle-1-T, cycle-2-A, ...].  The best quality from the 4 qualities for
     * the cycle is found, and then it is ASCII-ized by adding 64.
     * @param qualities Array of integer quality strings.
     * @param cycleNumber Which cycle to get quality for.
     * @param formatter For converting decimal strings to ints.
     * @return best quality for the given cycle.
     * @throws net.sf.picard.PicardException if the best quality ASCII value is > 255.
     */
    public static byte getSolexaQualityCharFromFourQualities(final String[] qualities, final int cycleNumber, final FormatUtil formatter) {
        // It apparently is the case that all 4 qualities might be negative, but this appears to correspond to
        // an no-called base.
        int bestQuality = Integer.MIN_VALUE;
        final int startOffset = (cycleNumber - 1) * 4;
        for (int i = startOffset; i < startOffset + 4; ++i) {
            final int quality = formatter.parseInt(qualities[i]);
            if (quality > bestQuality) {
                bestQuality = quality;
            }
        }
        final int qualityAsCharacter = bestQuality + SolexaQualityConverter.SOLEXA_ADDEND;
        if (qualityAsCharacter > 255) {
            throw new PicardException("Quality too large: " + bestQuality);
        }
        return (byte)(qualityAsCharacter & 0xff);
    }

    /** Describes adapters used on each pair of strands */
    public static enum AdapterPair {

        PAIRED_END("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",  //58 bases)
                   "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"), // 61 bases

        INDEXED ("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                 "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"), // note  8 N's  // 67 bases

        SINGLE_END ("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                    "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"),

        ALTERNATIVE_SINGLE_END("AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGACGATC",
                          "TCGTATGCCGTCTTCTGCTTG");

        final String fivePrime, threePrime, fivePrimeReadOrder;
        final byte[]  fivePrimeBytes, threePrimeBytes, fivePrimeReadOrderBytes;

        private AdapterPair(final String fivePrime, final String threePrime) {  // constructor
            this.threePrime = threePrime;
            this.fivePrime = fivePrime;
            this.fivePrimeReadOrder = SequenceUtil.reverseComplement(fivePrime);
            this.threePrimeBytes = StringUtil.stringToBytes(threePrime);
            this.fivePrimeBytes = StringUtil.stringToBytes(fivePrime);
            this.fivePrimeReadOrderBytes = StringUtil.stringToBytes(fivePrimeReadOrder);

        }

        public String get3PrimeAdapter(){ return threePrime; }
        public String get5PrimeAdapter(){ return fivePrime; }
        public String get3PrimeAdapterInReadOrder(){ return threePrime; }
        public String get5PrimeAdapterInReadOrder() { return fivePrimeReadOrder; }
        public byte[] get3PrimeAdapterBytes() { return threePrimeBytes; }
        public byte[] get5PrimeAdapterBytes() { return fivePrimeBytes; }
        public byte[] get3PrimeAdapterBytesInReadOrder() { return threePrimeBytes; }
        public byte[] get5PrimeAdapterBytesInReadOrder()  { return fivePrimeReadOrderBytes; }
    }

}
