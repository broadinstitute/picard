package picard.illumina;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.util.BarcodeEditDistanceQuery;
import picard.util.IlluminaUtil;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * BarcodeExtractor is used to match barcodes and collect barcode match metrics.
 */
public class BarcodeExtractor {
    private final Map<String, BarcodeMetric> metrics = new HashMap<>();
    private final BarcodeMetric noMatch;
    private final Set<ByteString> barcodesBytes;
    private final int maxNoCalls, maxMismatches, minMismatchDelta, minimumBaseQuality;
    private final DistanceMetric distanceMode;
    private final static int INITIAL_LOOKUP_SIZE = 4096;
    private final ConcurrentHashMap<ByteString, BarcodeMatch> barcodeLookupMap = new ConcurrentHashMap<>(INITIAL_LOOKUP_SIZE);

    public BarcodeExtractor(final Map<String, BarcodeMetric> barcodeToMetrics,
                            final BarcodeMetric noMatchMetric,
                            final ReadStructure readStructure,
                            final int maxNoCalls,
                            final int maxMismatches,
                            final int minMismatchDelta,
                            final int minimumBaseQuality,
                            final DistanceMetric distanceMode) {
        this.maxNoCalls = maxNoCalls;
        this.maxMismatches = maxMismatches;
        this.minMismatchDelta = minMismatchDelta;
        this.minimumBaseQuality = minimumBaseQuality;
        this.distanceMode = distanceMode;

        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode =  generateNoMatchBarcode(readStructure);
        final byte[][] perfectScores = new byte[readStructure.sampleBarcodes.length()][];
        int index = 0;
        for (final ReadDescriptor d : readStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                perfectScores[index] = new byte[d.length];
                Arrays.fill(perfectScores[index], (byte) 60);
            }
        }

        this.noMatch = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

        Set<BarcodeExtractor.ByteString> barcodesBytes = new HashSet<>(barcodeToMetrics.size());
        for (final BarcodeMetric metric : barcodeToMetrics.values()) {
            barcodesBytes.add(new BarcodeExtractor.ByteString(metric.barcodeBytes));
            metrics.put(metric.BARCODE_WITHOUT_DELIMITER, metric);
        }
        this.barcodesBytes = barcodesBytes;

        // Prepopulate the lookup map with all perfect barcodes
        for(BarcodeMetric metric : barcodeToMetrics.values()){
            BarcodeExtractor.BarcodeMatch match = calculateBarcodeMatch(metric.barcodeBytes, perfectScores, true);
            barcodeLookupMap.put(new BarcodeExtractor.ByteString(metric.barcodeBytes), match);
        }

        // Prepopulate all no call barcode match
        BarcodeExtractor.BarcodeMatch noCallMatch = calculateBarcodeMatch(noMatchMetric.barcodeBytes,
                perfectScores, true);

        barcodeLookupMap.put(new BarcodeExtractor.ByteString(noMatchMetric.barcodeBytes), noCallMatch);
    }

    public Map<String, BarcodeMetric> getMetrics() {
        return this.metrics;
    }

    public BarcodeMetric getNoMatchMetric() {
        return this.noMatch;
    }

    /**
     * Find the best barcode match for the given read sequence, and accumulate metrics
     *
     * NOTE: the returned BarcodeMatch object will contain mismatches mismatchesToSecondBest values that may be
     * inaccurate as long as the conclusion match/no-match isn't affected. for example, mismatches and mismatchesToSecondBest
     * may be smaller than their true value if mismatches is truly larger than maxMismatches.
     * Also, mismatchesToSecondBest might be smaller than its true value if its true value is greater than
     * mismatches + minMismatchDelta. This is due to an optimization which allows the distance calculation to stop once
     * the conclusion (Match or no-Match) can be reached.
     *
     * @param readSubsequences The basecalls for the read
     * @param qualityScores The quality scores for the reads
     * @param isInlineMatching true if we are matching the barcodes inline during conversion, false if we are matching
     *                         barcodes separately (ExtractIlluminaBarcodes)
     * @return perfect barcode string, if there was a match within tolerance, or null if not.
     */
    BarcodeMatch findBestBarcode(final byte[][] readSubsequences,
                                 final byte[][] qualityScores,
                                 final boolean isInlineMatching) {
        final boolean canUseLookupTable = areAllQualitiesAboveMinimum(qualityScores, minimumBaseQuality);
        if (canUseLookupTable) {
            final ByteString barcodesAsString = new ByteString(readSubsequences);
            BarcodeMatch match = barcodeLookupMap.get(barcodesAsString);
            if (match == null) {
              match = calculateBarcodeMatch(readSubsequences, qualityScores, isInlineMatching);
            }
            if (match.isMatched()) {
                barcodeLookupMap.put(barcodesAsString, match);
            }
            return match;
        }
        else {
            return calculateBarcodeMatch(readSubsequences, qualityScores, isInlineMatching);
        }
    }

    BarcodeMatch calculateBarcodeMatch(final byte[][] readSubsequences,
                                       final byte[][] qualityScores,
                                       final boolean isInlineMatching) {
        String bestBarcode = null;
        final BarcodeMatch match = new BarcodeMatch();

        int totalBarcodeReadBases = 0;
        int numNoCalls = 0; // NoCalls are calculated for all the barcodes combined

        for (final byte[] bc : readSubsequences) {
            totalBarcodeReadBases += bc.length;
            for (final byte b : bc) {
                if (SequenceUtil.isNoCall(b)) {
                    ++numNoCalls;
                }
                if (isInlineMatching && numNoCalls > maxNoCalls) {
                    match.mismatches = totalBarcodeReadBases;
                    match.barcode = "";
                    match.matched = false;
                    return match;
                }
            }
        }

        // PIC-506 When forcing all reads to match a single barcode, allow a read to match even if every
        // base is a mismatch.
        int numMismatchesInBestBarcode = totalBarcodeReadBases + 1;
        int numMismatchesInSecondBestBarcode = totalBarcodeReadBases + 1;

        for (final ByteString barcodeBytes : barcodesBytes) {
            // need to add maxMismatches + minMismatchDelta together since the result might get used as numMismatchesInSecondBestBarcode
            final BarcodeEditDistanceQuery barcodeEditDistanceQuery = new BarcodeEditDistanceQuery(barcodeBytes.bytes, readSubsequences, qualityScores,
                    minimumBaseQuality, Math.min(maxMismatches, numMismatchesInBestBarcode) + minMismatchDelta);
            final int numMismatches = distanceMode.distance(barcodeEditDistanceQuery);

            if (numMismatches < numMismatchesInBestBarcode) {
                if (bestBarcode != null) {
                    numMismatchesInSecondBestBarcode = numMismatchesInBestBarcode;
                }
                numMismatchesInBestBarcode = numMismatches;
                bestBarcode = barcodeBytes.toString();
            } else if (numMismatches < numMismatchesInSecondBestBarcode) {
                numMismatchesInSecondBestBarcode = numMismatches;
            }
        }

        match.matched = bestBarcode != null &&
                numNoCalls <= maxNoCalls &&
                numMismatchesInBestBarcode <= maxMismatches &&
                numMismatchesInSecondBestBarcode - numMismatchesInBestBarcode >= minMismatchDelta;

        match.mismatches = numMismatchesInBestBarcode;
        match.mismatchesToSecondBest = numMismatchesInSecondBestBarcode;

        if (match.matched) {
            match.barcode = bestBarcode;
        } else {
            // If we have something that's not a "match" but matches one barcode
            // slightly, we output that matching barcode in lower case
            if (!isInlineMatching && numNoCalls + numMismatchesInBestBarcode < totalBarcodeReadBases && bestBarcode != null) {
                match.barcode = bestBarcode.toLowerCase();
            } else {
                match.mismatches = totalBarcodeReadBases;
                match.barcode = "";
            }
        }

        return match;
    }

    static void updateMetrics(final BarcodeMatch match,
                              final boolean passingFilter,
                              final Map<String, BarcodeMetric> metrics,
                              final BarcodeMetric noMatchBarcodeMetric) {
        if (match.matched) {
            final BarcodeMetric matchMetric = metrics.get(match.barcode);
            ++matchMetric.READS;
            if (passingFilter) {
                ++matchMetric.PF_READS;
            }
            if (match.mismatches == 0) {
                ++matchMetric.PERFECT_MATCHES;
                if (passingFilter) {
                    ++matchMetric.PF_PERFECT_MATCHES;
                }
            } else if (match.mismatches == 1) {
                ++matchMetric.ONE_MISMATCH_MATCHES;
                if (passingFilter) {
                    ++matchMetric.PF_ONE_MISMATCH_MATCHES;
                }
            }
        } else {
            ++noMatchBarcodeMetric.READS;
            if (passingFilter) {
                ++noMatchBarcodeMetric.PF_READS;
            }
        }
    }

    /**
     * Checks to ensure that all quality values are greater that the given cutoff such that comparisons can
     * be done just using the bases without further reference to quality scores.
     */
    private static boolean areAllQualitiesAboveMinimum(final byte[][] qualityScores, final int minimumBaseQuality) {
        if (qualityScores == null) return true;

        for (final byte[] qs : qualityScores) {
            for (final byte q : qs) {
                if (q < minimumBaseQuality) {
                    return false;
                }
            }
        }

        return true;
    }

    public int getMinimumBaseQuality() {
        return this.minimumBaseQuality;
    }

    public static String[] generateNoMatchBarcode(ReadStructure inputReadStructure) {
        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[inputReadStructure.sampleBarcodes.length()];
        int index = 0;
        for (final ReadDescriptor d : inputReadStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                noMatchBarcode[index++] = StringUtil.repeatCharNTimes('N', d.length);
            }
        }
        return noMatchBarcode;
    }

    /**
     * Utility class to hang onto data about the best match for a given barcode
     */
    public static class BarcodeMatch {
        private boolean matched;
        private String barcode;
        private int mismatches;
        private int mismatchesToSecondBest;

        public boolean isMatched() {
            return matched;
        }

        public String getBarcode() {
            return barcode;
        }

        public int getMismatches() { return mismatches; }

        public int getMismatchesToSecondBest() { return mismatchesToSecondBest; }
    }

    /**
     * Class to give a byte[][] a hashcode and equals without copying the whole contents into a String.
     */
    private static final class ByteString {
        private final byte[][] bytes;
        private final int hash;

        public ByteString(byte[][] bytes) {
            this.bytes = new byte[bytes.length][];
            System.arraycopy(bytes, 0, this.bytes, 0, bytes.length);

            // Pre-compute the hash-code
            int h = 0;
            for (final byte[] bs : this.bytes) {
                for (final byte b : bs) {
                    h = 31 * h + b;
                }
            }

            this.hash = h;

        }

        @Override
        public final int hashCode() {
            return this.hash;
        }

        @Override
        public boolean equals(final Object obj) {
            try {
                final ByteString that = (ByteString) obj;
                if (this.hash != that.hash) return false;
                if (this.bytes.length != that.bytes.length) return false;
                for (int i = 0; i < this.bytes.length; ++i) {
                    if (!Arrays.equals(this.bytes[i], that.bytes[i])) return false;
                }
                return true;
            }
            catch (final Exception e) {
                return false;
            }
        }

        @Override
        public String toString() {
            StringBuilder barcodeBuilder = new StringBuilder();
            for(byte[] barcode : bytes) {
                barcodeBuilder.append(new String(barcode, StandardCharsets.UTF_8));
            }
            return barcodeBuilder.toString();
        }
    }
}
