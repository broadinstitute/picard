/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.samtools;

/**
 * Utility methods for encoding and decoding the SQ tag value of SAMRecord.
 * 
 * @author alecw@broadinstitute.org
 */
public class SQTagUtil {
    /**
     * The ordinals of these are stored in the high-order 2 bits of each byte of the SQ tag.
     * Note that these have the convenient property that the binary complement of each ordinal, masked to
     * the two low-order bits, is the complementary base.
     */
    public enum SQBase {
        SQ_A('A'), SQ_C('C'), SQ_G('G'), SQ_T('T');
        private final Character base;

        SQBase(final Character base) {
            this.base = base;
        }

        public Character getBase() {
            return base;
        }
    }

    /**
     * For complementing SQBase ordinals.
     */
    private static final int COMPLEMENT_MASK = 3;

    private static final int QUALITY_MASK = 0x3f;
    public static final byte MAX_QUALITY = QUALITY_MASK;
    private static final int BASE_INDEX_SHIFT = 6;

    /**
     * Convert a pair of likelihoods into a value suitable for passing to baseAndProbDiffToSqValue.
     * @param secondBestLikelihood Probability of the 2nd-best base call.  1 > secondBestLikelihood > thirdBestLikelihood.
     * @param thirdBestLikelihood Probability of the 3rd-best base call.  thirdBestLikelihood > 0.
     * @return ratio of input probabilities for storing in SQ tag.
     */
    public static byte sqScaledProbabilityRatio(final double secondBestLikelihood, final double thirdBestLikelihood) {
        if (secondBestLikelihood >= 1.0 || thirdBestLikelihood <= 0 || thirdBestLikelihood > secondBestLikelihood) {
            throw new IllegalArgumentException("Likelihoods out of range.  second best: " + secondBestLikelihood +
            "; third best: " + thirdBestLikelihood);
        }
        // Cap value at QUALITY_MASK
        return (byte)(Math.min(Math.round(-10.0 * Math.log10(thirdBestLikelihood/secondBestLikelihood)), QUALITY_MASK));
    }

    /**
     * Compress a base and a log probabiliy difference (-10log10(p3/p2)) into
     * a single byte so that it can be output in a SAMRecord's SQ field.
     *
     * @param base  the 2nd-best base.
     * @param probRatio   the log probability difference between the secondary and tertiary bases (-10log10(p3/p2)),
     * rounded to an integer and capped so it fits in 6 bits.
     * @return a byte containing the index and the log probability difference.
     */
    public static byte baseAndProbDiffToSqValue(final SQBase base, final byte probRatio) {
        return baseAndProbDiffToSqValue(base.ordinal(), probRatio);
    }

    /**
     * Compress a base and a log probabiliy difference (-10log10(p3/p2)) into
     * a single byte so that it can be output in a SAMRecord's SQ field.
     *
     * @param base  the 2nd-best base (A=0, C=1, G=2, T=3).
     * @param probRatio   the log probability difference between the secondary and tertiary bases (-10log10(p3/p2)),
     * rounded to an integer and capped so it fits in 6 bits.  If this value is > MAX_QUALITY, it is truncated to that.
     * @return a byte containing the index and the log probability difference.
     */
    public static byte baseAndProbDiffToSqValue(final int base, final byte probRatio) {
        return (byte)((base << BASE_INDEX_SHIFT) | Math.min(probRatio, QUALITY_MASK));
    }

    /**
     * Retrieve SQ-scaled probability ratio from SQ value.
     * @param sqValue
     * @return the log probability difference between the secondary and tertiary bases (-10log10(p3/p2)).
     */
    public static byte sqValueToProbRatio(final byte sqValue) {
        return (byte)(sqValue & QUALITY_MASK);
    }

    /**
     * Retrieve the 2nd-best base call from SQ value.
     * @param sqValue
     * @return 2nd-best base call.
     */
    public static SQBase sqValueToBase(final byte sqValue) {
        return SQBase.values()[sqValueToBaseOrdinal(sqValue)];
    }

    /**
     * Retrieve the 2nd-best base call from SQ value.
     * @param sqValue
     * @return Ordinal of 2nd-best base call.
     */
    public static int sqValueToBaseOrdinal(final byte sqValue) {
        return (sqValue & 0xff) >>> BASE_INDEX_SHIFT;
    }


    /**
     * Reverses and complements the sqValues in place.
     * @param sqArray Array of SQ-values, with 2nd-best base in high-order 2 bits, and probability diff
     * in low-order 6 bits.
     */
    public static void reverseComplementSqArray(final byte[] sqArray) {
        final int lastIndex = sqArray.length - 1;

        int i, j;
        for (i=0, j=lastIndex; i<j; ++i, --j) {
            final byte tmp = complementSqValue(sqArray[i]);
            sqArray[i] = complementSqValue(sqArray[j]);
            sqArray[j] = tmp;
        }
        if (sqArray.length % 2 == 1) {
            sqArray[i] = complementSqValue(sqArray[i]);
        }
    }

    private static byte complementSqValue(final byte sqValue) {
        final byte probDiff = sqValueToProbRatio(sqValue);
        final int baseOrdinal = sqValueToBaseOrdinal(sqValue);
        final int complementOrdinal = COMPLEMENT_MASK & ~baseOrdinal;
        return baseAndProbDiffToSqValue(complementOrdinal, probDiff);
    }

}
