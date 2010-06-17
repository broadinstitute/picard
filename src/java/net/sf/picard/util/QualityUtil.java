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

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;

/**
 * Utility class for working with quality scores and error probabilities.
 *
 * @author Tim Fennell
 */
public final class QualityUtil {
    private static final double[] errorProbabilityByPhredScore;

    static {
        errorProbabilityByPhredScore = new double[101];
        for (int i=0; i<errorProbabilityByPhredScore.length; ++i) {
            errorProbabilityByPhredScore[i] = 1d/Math.pow(10d, i/10d);
        }
    }

    /** Given a phred score between 0 and 100 returns the probability of error. */
    public static double getErrorProbabilityFromPhredScore(final int i) {
        return errorProbabilityByPhredScore[i];
    }

    /** Gets the phred score for any given probability of error. */
    public static int getPhredScoreFromErrorProbability(final double probability) {
        return (int) Math.round(-10 * Math.log10(probability));
    }

    /** Gets the phred score given the specified observations and errors. */
    public static int getPhredScoreFromObsAndErrors(final double observations, final double errors) {
        return getPhredScoreFromErrorProbability(errors / observations);
    }

    /**
     * Calculates the sum of error probabilities for all read bases in the SAM record. Takes
     * the SAM record as opposed to the qualities directly so that it can make sure to count
     * no-calls as 1 instead of what the quality score says.
     * */
    public static double sumOfErrorProbabilities(final SAMRecord rec) {
        final byte[] bases = rec.getReadBases();
        final byte[] quals = rec.getBaseQualities();

        double sum = 0;

        for (int i=0; i<bases.length; ++i) {
            if (SequenceUtil.isNoCall(bases[i])) ++sum;
            else sum += QualityUtil.getErrorProbabilityFromPhredScore(quals[i]);
        }

        return sum;
    }
}
