package picard.illumina;

import htsjdk.samtools.util.SequenceUtil;
import picard.util.IlluminaUtil;

import java.util.Map;

public class BarcodeMatch {
    private boolean matched = false;
    private String barcode;
    private int mismatches;
    private int mismatchesToSecondBest;

    public boolean isMatched() {
        return matched;
    }

    public void setMatched(boolean matched) {
        this.matched = matched;
    }

    public String getBarcode() {
        return barcode;
    }

    public void setBarcode(String barcode) {
        this.barcode = barcode;
    }

    public int getMismatches() {
        return mismatches;
    }

    public void setMismatches(int mismatches) {
        this.mismatches = mismatches;
    }

    public int getMismatchesToSecondBest() {
        return mismatchesToSecondBest;
    }

    public void setMismatchesToSecondBest(int mismatchesToSecondBest) {
        this.mismatchesToSecondBest = mismatchesToSecondBest;
    }


    /**
     * Find the best barcode match for the given read sequence, and accumulate metrics
     *
     * @param readSubsequences   portion of read containing barcode
     * @param passingFilter      PF flag for the current read
     * @param maxNoCalls         The maximum number of no calls to allow in a match.
     * @param maxMismatches      The maximum number of mismatched calls to allow in a match.
     * @param minMismatchDelta   The minimum mismatch difference between the best and second best matches.
     * @param minimumBaseQuality The minimum base quality to allow for a matching base call.
     * @return perfect barcode string, if there was a match within tolerance, or null if not.
     */
    public static BarcodeMatch findBestBarcodeAndUpdateMetrics(final byte[][] readSubsequences,
                                                               final byte[][] qualityScores,
                                                               final boolean passingFilter,
                                                               final Map<String, BarcodeMetric> metrics,
                                                               final BarcodeMetric noMatchBarcodeMetric,
                                                               int maxNoCalls, int maxMismatches,
                                                               int minMismatchDelta,
                                                               int minimumBaseQuality) {
        BarcodeMetric bestBarcodeMetric = null;
        int totalBarcodeReadBases = 0;
        int numNoCalls = 0; // NoCalls are calculated for all the barcodes combined

        for (final byte[] bc : readSubsequences) {
            totalBarcodeReadBases += bc.length;
            for (final byte b : bc) if (SequenceUtil.isNoCall(b)) ++numNoCalls;
        }

        // PIC-506 When forcing all reads to match a single barcode, allow a read to match even if every
        // base is a mismatch.
        int numMismatchesInBestBarcode = totalBarcodeReadBases + 1;
        int numMismatchesInSecondBestBarcode = totalBarcodeReadBases + 1;

        for (final BarcodeMetric barcodeMetric : metrics.values()) {
            final int numMismatches = countMismatches(barcodeMetric.barcodeBytes, readSubsequences, qualityScores,
                    minimumBaseQuality);
            if (numMismatches < numMismatchesInBestBarcode) {
                if (bestBarcodeMetric != null) {
                    numMismatchesInSecondBestBarcode = numMismatchesInBestBarcode;
                }
                numMismatchesInBestBarcode = numMismatches;
                bestBarcodeMetric = barcodeMetric;
            } else if (numMismatches < numMismatchesInSecondBestBarcode) {
                numMismatchesInSecondBestBarcode = numMismatches;
            }
        }

        final boolean matched = bestBarcodeMetric != null &&
                numNoCalls <= maxNoCalls &&
                numMismatchesInBestBarcode <= maxMismatches &&
                numMismatchesInSecondBestBarcode - numMismatchesInBestBarcode >= minMismatchDelta;

        final BarcodeMatch match = new BarcodeMatch();

        // If we have something that's not a "match" but matches one barcode
        // slightly, we output that matching barcode in lower case
        if (numNoCalls + numMismatchesInBestBarcode < totalBarcodeReadBases) {
            match.setMismatches(numMismatchesInBestBarcode);
            match.setMismatchesToSecondBest(numMismatchesInSecondBestBarcode);
            match.setBarcode(bestBarcodeMetric != null ? bestBarcodeMetric.BARCODE.toLowerCase().replaceAll(IlluminaUtil.BARCODE_DELIMITER, "") : null);
        } else {
            match.setMismatches(totalBarcodeReadBases);
            match.setBarcode("");
        }

        if (matched) {
            ++bestBarcodeMetric.READS;
            if (passingFilter) {
                ++bestBarcodeMetric.PF_READS;
            }
            if (numMismatchesInBestBarcode == 0) {
                ++bestBarcodeMetric.PERFECT_MATCHES;
                if (passingFilter) {
                    ++bestBarcodeMetric.PF_PERFECT_MATCHES;
                }
            } else if (numMismatchesInBestBarcode == 1) {
                ++bestBarcodeMetric.ONE_MISMATCH_MATCHES;
                if (passingFilter) {
                    ++bestBarcodeMetric.PF_ONE_MISMATCH_MATCHES;
                }
            }

            match.setMatched(true);
            match.setBarcode(bestBarcodeMetric.BARCODE.replaceAll(IlluminaUtil.BARCODE_DELIMITER, ""));
        } else {
            ++noMatchBarcodeMetric.READS;
            if (passingFilter) {
                ++noMatchBarcodeMetric.PF_READS;
            }
        }

        return match;
    }

    /**
     * Compare barcode sequence to bases from read
     *
     * @return how many bases did not match
     */
    private static int countMismatches(final byte[][] barcodeBytes, final byte[][] readSubsequence,
                                       final byte[][] qualities, int minimumBaseQuality) {
        int numMismatches = 0;
        // Read sequence and barcode length may not be equal, so we just use the shorter of the two
        for (int j = 0; j < barcodeBytes.length; j++) {
            final int basesToCheck = Math.min(barcodeBytes[j].length, readSubsequence[j].length);
            for (int i = 0; i < basesToCheck; ++i) {
                if (!SequenceUtil.isNoCall(readSubsequence[j][i])) {
                    if (!SequenceUtil.basesEqual(barcodeBytes[j][i], readSubsequence[j][i])) ++numMismatches;
                    else if (qualities != null && qualities[j][i] < minimumBaseQuality) ++numMismatches;
                }
            }
        }
        return numMismatches;
    }
}