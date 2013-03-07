/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package net.sf.picard.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class SamPairUtilTest {


    @Test(dataProvider = "testGetPairOrientation")
    public void testGetPairOrientation(final String testName,
                                       final int read1Start, final int read1Length, final boolean read1Reverse,
                                       final int read2Start, final int read2Length, final boolean read2Reverse,
                                       final SamPairUtil.PairOrientation expectedOrientation) {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 100000000));
        final SAMRecord rec1 = makeSamRecord(header, read1Start, read1Length, read1Reverse, true);
        final SAMRecord rec2 = makeSamRecord(header, read2Start, read2Length, read2Reverse, false);
        SamPairUtil.setMateInfo(rec1, rec2, header);
        Assert.assertEquals(SamPairUtil.getPairOrientation(rec1), expectedOrientation, testName + " first end");
        Assert.assertEquals(SamPairUtil.getPairOrientation(rec2), expectedOrientation, testName + " second end");
    }

    private SAMRecord makeSamRecord(final SAMFileHeader header, final int alignmentStart, final int readLength,
                                    final boolean reverse, final boolean firstOfPair) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReferenceIndex(0);
        final StringBuilder sb = new StringBuilder();
        final byte[] quals = new byte[readLength];
        for (int i = 0; i < readLength; ++i) {
            sb.append("A");
            quals[i] = 20;
        }
        rec.setReadString(sb.toString());
        rec.setBaseQualities(quals);
        rec.setAlignmentStart(alignmentStart);
        rec.setCigarString(readLength + "M");
        rec.setReadPairedFlag(true);
        rec.setReadNegativeStrandFlag(reverse);
        if (firstOfPair) rec.setFirstOfPairFlag(true);
        else rec.setSecondOfPairFlag(true);
        return rec;
    }

    @DataProvider(name = "testGetPairOrientation")
    public Object[][] testGetPairOrientationDataProvider() {
        /**
         * @param testName
         * @param read1Start
         * @param read1Length
         * @param read1Reverse
         * @param read2Start
         * @param read2Length
         * @param read2Reverse
         * @param expectedOrientation
         */
        return new Object[][] {
                {"normal innie", 1, 100, false, 500, 100, true, SamPairUtil.PairOrientation.FR},
                {"overlapping innie", 1, 100, false, 50, 100, true, SamPairUtil.PairOrientation.FR},
                {"second end enclosed innie", 1, 100, false, 50, 50, true, SamPairUtil.PairOrientation.FR},
                {"first end enclosed innie", 1, 50, false, 1, 100, true, SamPairUtil.PairOrientation.FR},
                {"completely overlapping innie", 1, 100, false, 1, 100, true, SamPairUtil.PairOrientation.FR},

                {"normal outie", 1, 100, true, 500, 100, false, SamPairUtil.PairOrientation.RF},
                {"nojump outie", 1, 100, true, 101, 100, false, SamPairUtil.PairOrientation.RF},

                {"forward tandem", 1, 100, true, 500, 100, true, SamPairUtil.PairOrientation.TANDEM},
                {"reverse tandem", 1, 100, false, 500, 100, false, SamPairUtil.PairOrientation.TANDEM},
                {"overlapping forward tandem", 1, 100, true, 50, 100, true, SamPairUtil.PairOrientation.TANDEM},
                {"overlapping reverse tandem", 1, 100, false, 50, 100, false, SamPairUtil.PairOrientation.TANDEM},
                {"second end enclosed forward tandem", 1, 100, true, 50, 50, true, SamPairUtil.PairOrientation.TANDEM},
                {"second end enclosed reverse tandem", 1, 100, false, 50, 50, false, SamPairUtil.PairOrientation.TANDEM},
                {"first end enclosed forward tandem", 1, 50, true, 1, 100, true, SamPairUtil.PairOrientation.TANDEM},
                {"first end enclosed reverse tandem", 1, 50, false, 1, 100, false, SamPairUtil.PairOrientation.TANDEM},
        };
    }
}
