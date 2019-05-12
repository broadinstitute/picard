package picard.illumina;

import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.util.StringDistanceUtils;

/**
 * An enum controlling the type of distance metric used for barcode comparisons.
 */
public enum DistanceMetric implements CommandLineParser.ClpEnum {

    HAMMING("Simple edit distance without considering indels at all") {
        @Override
        int countMismatches(final byte[][] barcodeBytes, final byte[][] readSubsequence, final byte[][] qualities, final int minimumBaseQuality) {
            return StringDistanceUtils.countMismatches(barcodeBytes, readSubsequence, qualities, minimumBaseQuality);
        }
    },

    LEVENSHTEIN("Edit differences by global Levenshtein distance, allowing for gaps and aligning the entire strings.") {
        @Override
        int countMismatches(final byte[][] barcodeBytes, final byte[][] readSubsequence, final byte[][] qualities, final int minimumBaseQuality) {
            return StringDistanceUtils.countMismatchesLevenshtein(barcodeBytes, readSubsequence, qualities, minimumBaseQuality);
        }
    },

    /**
     * Based on https://www.pnas.org/content/115/27/E6217 FREE stands for "Filled/truncated Right End Edit"
     */
    FREE("Edit difference that doesn't include gaps at the end of the alignment that are forced by earlier gaps. Most appropriate for fixed length barcodes.") {
        @Override
        int countMismatches(final byte[][] barcodeBytes, final byte[][] readSubsequence, final byte[][] qualities, final int minimumBaseQuality) {
            return StringDistanceUtils.countMismatchesWithUnforcedIndels(barcodeBytes, readSubsequence, qualities, minimumBaseQuality);
        }
    };

    private final String helpdoc;

    DistanceMetric(String helpdoc){
        this.helpdoc = helpdoc;
    }

    abstract int countMismatches(byte[][] barcodeBytes, byte[][] readSubsequence, byte[][] qualities, int minimumBaseQuality);

    @Override
    public String getHelpDoc() {
        return this.helpdoc;
    }}
