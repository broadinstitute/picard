/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

import java.util.BitSet;

/**
 * Constants and methods used by BAM and Tribble indices
 */
public class GenomicIndexUtil {
    /**
     * Reports the total amount of genomic data that any bin can index.
     */
    public static final int BIN_GENOMIC_SPAN = 512*1024*1024;

    /**
     * What is the starting bin for each level?
     */
    public static final int[] LEVEL_STARTS = {0,1,9,73,585,4681};

    /**
     * Reports the maximum number of bins that can appear in a binning index.
     */
    public static final int MAX_BINS = 37450;   // =(8^6-1)/7+1

    public static final int MAX_LINEAR_INDEX_SIZE = MAX_BINS+1-LEVEL_STARTS[LEVEL_STARTS.length-1];


    /**
     * E.g. for a SAMRecord with no genomic coordinate.
     */
    public static final int UNSET_GENOMIC_LOCATION = 0;

    /**
     * calculate the bin given an alignment in [beg,end)
     * Copied from SAM spec.
     * @param beg 0-based start of read (inclusive)
     * @param end 0-based end of read (exclusive)
     */
    static int reg2bin(final int beg, int end)
    {
        --end;

        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return  ((1<<9)-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return  ((1<<6)-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return  ((1<<3)-1)/7 + (beg>>26);
        return 0;
    }

    // TODO: It is disturbing that reg2bin is 0-based, but regionToBins is 1-based.
    // TODO: It is also suspicious that regionToBins decrements endPos.  Test it!
    // TODO: However end is decremented in reg2bin so perhaps there is no conflict.
    /**
     * Get candidate bins for the specified region
     * @param startPos 1-based start of target region, inclusive.
     * @param endPos 1-based end of target region, inclusive.
     * @return bit set for each bin that may contain SAMRecords in the target region.
     */
    public static BitSet regionToBins(final int startPos, final int endPos) {
        final int maxPos = 0x1FFFFFFF;
        final int start = (startPos <= 0) ? 0 : (startPos-1) & maxPos;
        final int end = (endPos <= 0) ? maxPos : (endPos-1) & maxPos;
        if (start > end) {
            return null;
        }
        int k;
        final BitSet bitSet = new BitSet(GenomicIndexUtil.MAX_BINS);
        bitSet.set(0);
        for (k =    1 + (start>>26); k <=    1 + (end>>26); ++k) bitSet.set(k);
        for (k =    9 + (start>>23); k <=    9 + (end>>23); ++k) bitSet.set(k);
        for (k =   73 + (start>>20); k <=   73 + (end>>20); ++k) bitSet.set(k);
        for (k =  585 + (start>>17); k <=  585 + (end>>17); ++k) bitSet.set(k);
        for (k = 4681 + (start>>14); k <= 4681 + (end>>14); ++k) bitSet.set(k);
        return bitSet;
    }

}
