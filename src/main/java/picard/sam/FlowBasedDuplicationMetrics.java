/*
 * The MIT License
 *
 * Copyright (c) 2022 The Broad Institute
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

package picard.sam;

import htsjdk.samtools.SAMRecord;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;
import picard.util.MathUtil;

public class FlowBasedDuplicationMetrics extends DuplicationMetrics {

    /*
     * count of single end reads where the exact fragment length is known (i.e. clipped)
     */
    @MergeByAdding
    public long UNPAIRED_WITH_TLEN;

    /*
     * count of single end duplicates where the exact fragment length is
     * unknown (quality trimmed, not clipped)
     */
    @MergeByAdding
    public long UNPAIRED_DUPS_WITHOUT_TLEN;

    /*
     * count of duplicates where both ends are known
     */
    @MergeByAdding
    public long UNPAIRED_DUPS_WITH_TLEN;

    /**
     * The fraction of duplicated reads out of all reads with exact
     * fragment length unknown
     */
    @NoMergingIsDerived
    public Double UNPAIRED_DUP_RATE_WITHOUT_TLEN;

    /**
     * The fraction of duplicated reads out of all reads with exact fragment
     * length known
     */
    @NoMergingIsDerived
    public Double UNPAIRED_DUP_RATE_WITH_TLEN;


    @Override
    public void calculateDerivedFields() {
        super.calculateDerivedFields();

        UNPAIRED_DUP_RATE_WITHOUT_TLEN = MathUtil.divide(UNPAIRED_DUPS_WITHOUT_TLEN, UNPAIRED_READS_EXAMINED - UNPAIRED_WITH_TLEN);
        UNPAIRED_DUP_RATE_WITH_TLEN = MathUtil.divide(UNPAIRED_DUPS_WITH_TLEN, UNPAIRED_WITH_TLEN);
    }

    public void addDuplicateReadToMetrics(final SAMRecord rec) {
        super.addDuplicateReadToMetrics(rec);

        if (!rec.isSecondaryOrSupplementary() && !rec.getReadUnmappedFlag()) {
            if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                if ( ReadEndsForMarkDuplicates.isSingleEndReadKnownFragment(rec) ) {
                    ++UNPAIRED_DUPS_WITH_TLEN;
                } else {
                    ++UNPAIRED_DUPS_WITHOUT_TLEN;
                }
            }
        }
    }

    public void addReadToLibraryMetrics(final SAMRecord rec) {

        super.addReadToLibraryMetrics(rec);

        if (ReadEndsForMarkDuplicates.isSingleEndReadKnownFragment(rec)) {
            ++UNPAIRED_WITH_TLEN;
        }
    }

}
