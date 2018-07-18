/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class UmiUtilTest {

    @DataProvider(name = "topStrandDataProvider")
    private Object[][] testIsTopStrandDataProvider() {
        return new Object[][]{
                {100, 200, true, false, true, true},    // Read 1 in F1R2 pair
                {100, 200, false, false, true, false},  // Read 2 in F1R2 pair
                {100, 200, true, true, false, true},    // Read 1 in F2R1 pair
                {200, 100, false, false, true, true},   // Read 2 in F2R1 pair
                {100, 200, true, false, false, true},   // Read 1 in F1F2 pair
                {200, 100, false, false, false, true},  // Read 2 in F1F2 pair
                {100, 200, true, true, true, true},     // Read 1 in R1R2 pair
                {200, 100, false, true, true, true},    // Read 2 in R1R2 pair
        };
    }

    @Test(dataProvider = "topStrandDataProvider")
    public void testIsTopStrand(final int alignmentStart, final int mateAlignmentStart, final boolean firstOfPairFlag,
                              final boolean negativeStrandFlag, final boolean mateNegativeStrandFlag, final boolean topStrand) {

        SAMFileHeader header = new SAMFileHeader();
        SAMRecord rec = new SAMRecord(header);

        rec.setReadPairedFlag(true);

        rec.setCigarString("10M");
        rec.setAttribute("MC", "10M");
        rec.setAlignmentStart(alignmentStart);
        rec.setMateAlignmentStart(mateAlignmentStart);

        rec.setFirstOfPairFlag(firstOfPairFlag);
        rec.setReadNegativeStrandFlag(negativeStrandFlag);
        rec.setMateNegativeStrandFlag(mateNegativeStrandFlag);

        Assert.assertEquals(UmiUtil.isTopStrand(rec), topStrand);
    }
}
