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
package net.sf.picard.util;

/**
 * Optimized method for converting Solexa ASCII qualities into Phred scores.
 * Pre-computes all values in order to eliminate repeated computation.
 */
public class SolexaQualityConverter {

    /**
     * This value is added to a Solexa quality score to make it printable ASCII
     */
    public static final int SOLEXA_ADDEND = 64;

    /**
     * This value is added to a Phred scord to make it printable ASCII
     */
    public static final int PHRED_ADDEND = 33;

    private static SolexaQualityConverter singleton = null;

    public static synchronized SolexaQualityConverter getSingleton()  {
        if (singleton == null) {
            singleton = new SolexaQualityConverter();
        }
        return singleton;
    }

    /**
     * Mapping from ASCII value in Gerald export file to phred score
     */
    private final byte[] phredScore = new byte[256];

    private SolexaQualityConverter() {
        for (int i = 0; i < SOLEXA_ADDEND; ++i) {
            phredScore[i] = 0;
        }
        for (int i = SOLEXA_ADDEND; i < phredScore.length; ++i) {
            phredScore[i] = convertSolexaQualityCharToPhredBinary(i);
        }
    }


    /** Converts a solexa character quality into a phred numeric quality. */
    private byte convertSolexaQualityCharToPhredBinary(final int solexaQuality) {
        return (byte) Math.round(10d * Math.log10(1d+Math.pow(10d, (solexaQuality - SOLEXA_ADDEND)/10d)));
    }

    /**
     * Convert a solexa quality ASCII character into a phred score.
     */
    public byte solexaCharToPhredBinary(final byte solexaQuality) {
        return phredScore[solexaQuality];
    }

    /**
     * @return a byte array that can be indexed by Solexa ASCII quality, with value
     * of corresponding Phred score.  Elements 0-63 are invalid because Solexa qualities
     * should all be >= 64.  Do not modify this array!
     */
    public byte[] getSolexaToPhredConversionTable() {
        return phredScore;
    }

    /**
     * Decodes an array of solexa quality ASCII chars into Phred numeric space.
     * Decode in place in order to avoid extra object allocation.
     */
    public void convertSolexaQualityCharsToPhredBinary(final byte[] solexaQuals) {
        for (int i=0; i<solexaQuals.length; ++i) {
            solexaQuals[i] = phredScore[solexaQuals[i]];
        }
    }

    /**
     * Decodes an array of solexa quality ASCII chars into Phred ASCII space.
     * Decode in place in order to avoid extra object allocation.
     */
    public void convertSolexaQualityCharsToPhredChars(final byte[] solexaQuals) {
        for (int i=0; i<solexaQuals.length; ++i) {
            solexaQuals[i] = (byte)((phredScore[solexaQuals[i]] + PHRED_ADDEND) & 0xff);
        }
    }

    /**
     * Casava 1.3 stores phred-scaled qualities, but non-standard because they have 64 added to them
     * rather than the standard 33.
     * @param solexaQuals qualities are converted in place.
     */
    public void convertSolexa_1_3_QualityCharsToPhredBinary(final byte[] solexaQuals) {
        for (int i=0; i<solexaQuals.length; ++i) {
            solexaQuals[i] -= SOLEXA_ADDEND;
        }
    }

    public void convertSolexa_1_3_QualityCharsToPhredBinary(int offset, int length, final byte[] solexaQuals) {
        final int limit = offset + length;
        for (int i=offset; i < limit; ++i) {
            solexaQuals[i] -= SOLEXA_ADDEND;
        }
    }
}
