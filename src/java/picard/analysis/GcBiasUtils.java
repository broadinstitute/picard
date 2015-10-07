/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import java.io.File;

/** Utilities to calculate GC Bias
 * Created by kbergin on 9/23/15.
 */
public class GcBiasUtils {

    /////////////////////////////////////////////////////////////////////////////
    // Calculates GC as a number from 0 to 100 in the specified window.
    // If the window includes more than five no-calls then -1 is returned.
    /////////////////////////////////////////////////////////////////////////////
    public static int calculateGc(final byte[] bases, final int startIndex, final int endIndex, final CalculateGcState state) {
        if (state.init) {
            state.init = false;
            state.gcCount = 0;
            state.nCount = 0;
            for (int i = startIndex; i < endIndex; ++i) {
                final byte base = bases[i];
                if (SequenceUtil.basesEqual(base, (byte)'G') || SequenceUtil.basesEqual(base, (byte)'C')) ++state.gcCount;
                else if (SequenceUtil.basesEqual(base, (byte)'N')) ++state.nCount;
            }
        } else {
            final byte newBase = bases[endIndex - 1];
            if (SequenceUtil.basesEqual(newBase, (byte)'G') || SequenceUtil.basesEqual(newBase, (byte)'C')) ++state.gcCount;
            else if (newBase == 'N') ++state.nCount;

            if (SequenceUtil.basesEqual(state.priorBase, (byte)'G') || SequenceUtil.basesEqual(state.priorBase, (byte)'C')) --state.gcCount;
            else if (SequenceUtil.basesEqual(state.priorBase, (byte)'N')) --state.nCount;
        }
        state.priorBase = bases[startIndex];
        if (state.nCount > 4) return -1;
        else return (state.gcCount * 100) / (endIndex - startIndex);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculate number of 100bp windows in the refBases passed in that fall into
    // each gc content bin (0-100% gc)
    /////////////////////////////////////////////////////////////////////////////
    public static int[] calculateRefWindowsByGc(final int windows, final File referenceSequence, final int windowSize) {
        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceSequence);
        ReferenceSequence ref;

        final int [] windowsByGc = new int [windows];

        while ((ref = refFile.nextSequence()) != null) {
            final byte[] refBases = ref.getBases();
            StringUtil.toUpperCase(refBases);
            final int refLength = refBases.length;
            final int lastWindowStart = refLength - windowSize;

            final CalculateGcState state = new GcBiasUtils().new CalculateGcState();

            for (int i = 1; i < lastWindowStart; ++i) {
                final int windowEnd = i + windowSize;
                final int gcBin = calculateGc(refBases, i, windowEnd, state);
                if (gcBin != -1) windowsByGc[gcBin]++;
            }
        }

        return windowsByGc;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Calculate all the GC values for all windows
    /////////////////////////////////////////////////////////////////////////////
    public static byte [] calculateAllGcs(final byte[] refBases, final int lastWindowStart, final int windowSize) {

        final CalculateGcState state = new GcBiasUtils().new CalculateGcState();

        final int refLength = refBases.length;
        final byte[] gc = new byte[refLength + 1];

        for (int i = 1; i < lastWindowStart; ++i) {
            final int windowEnd = i + windowSize;
            final int windowGc = calculateGc(refBases, i, windowEnd, state);
            gc[i] = (byte) windowGc;
        }

        return gc;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Keeps track of current GC calculation state
    /////////////////////////////////////////////////////////////////////////////
    class CalculateGcState {
        boolean init = true;
        int nCount;
        int gcCount;
        byte priorBase;
    }
}
