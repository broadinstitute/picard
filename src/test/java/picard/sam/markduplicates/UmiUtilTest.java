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
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

public class UmiUtilTest {

    @DataProvider(name = "topStrandDataProvider")
    private Object[][] testIsTopStrandDataProvider() {
        final boolean read1 = true;
        final boolean read2 = false;
        return new Object[][]{
                {0, 100, 0, 200, read1, false, true, true, true, true},    // Read 1 in F1R2 pair, reads point inwards, top strand
                {0, 200, 0, 100, read2, true, false, true, true, true},    // Read 2 in F1R2 pair, reads point inwards, top strand
                {0, 200, 0, 100, read1, false, true, true, true, false},   // Read 1 in F1R2 pair, reads point outwards
                {0, 100, 0, 200, read2, true, false, true, true, false},   // Read 2 in F1R2 pair, reads point outwards
                {0, 100, 0, 200, read1, true, false, true, true, true},    // Read 1 in F2R1 pair
                {0, 200, 0, 100, read2, false, false, true, true, true},   // Read 2 in F2R1 pair
                {0, 100, 0, 200, read1, false, false, true, true, true},   // Read 1 in F1F2 pair
                {0, 200, 0, 100, read2, false, false, true, true, true},   // Read 2 in F1F2 pair
                {0, 100, 0, 200, read1, true, true, true, true, true},     // Read 1 in R1R2 pair
                {0, 200, 0, 100, read2, true, true, true, true, true},     // Read 2 in R1R2 pair
                {0, 100, 1, 200, read1, false, true, true, true, true},    // Read 1 in F1R2 chimera
                {0, 100, 1, 200, read2, true, false, true, true, false},   // Read 2 in F1R2 chimera
                {0, 100, 1, 200, read1, true, false, true, true, true},    // Read 1 in F2R1 chimera
                {0, 100, 1, 200, read2, false, false, true, true, false},  // Read 2 in F2R1 chimera
                {0, 100, 1, 200, read2, false, true, true, true, false},   // Read 2 in F2R1 chimera
                {0, 100, 0, 200, read1, false, true, false, true, true},   // Read 1 in F1R2 pair, top strand, read not mapped
                {0, 100, 0, 200, read1, false, true, true, false, false},  // Read 1 in F1R2 pair, top strand, mate not mapped
                {0, 100, 0, 200, read1, false, true, false, false, true},  // Read 1 in F1R2 pair, top strand, neither mapped
        };
    }

    @Test(dataProvider = "topStrandDataProvider")
    public void testIsTopStrand(final int referenceIndex, final int alignmentStart, final int mateReferenceIndex, final int mateAlignmentStart,
                                final boolean firstOfPairFlag, final boolean negativeStrandFlag, final boolean mateNegativeStrandFlag,
                                final boolean mapped, final boolean mateMapped,
                                final boolean topStrand) {

        final int readLength = 15;
        final int contigLength = 500;
        SAMFileHeader header = new SAMFileHeader();
        SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary();

        sequenceDictionary.addSequence(new SAMSequenceRecord("chr1", contigLength));
        sequenceDictionary.addSequence(new SAMSequenceRecord("chr2", contigLength));

        System.out.println(sequenceDictionary.getSequences());

        header.setSequenceDictionary(sequenceDictionary);

        SAMRecord rec = new SAMRecord(header);

        rec.setReadUnmappedFlag(!mapped);
        rec.setMateUnmappedFlag(!mateMapped);
        rec.setReadPairedFlag(true);

        rec.setCigarString(readLength + "M");
        rec.setAttribute("MC", readLength + "M");

        rec.setReferenceIndex(referenceIndex);
        rec.setAlignmentStart(alignmentStart);

        rec.setMateReferenceIndex(mateReferenceIndex);
        rec.setMateAlignmentStart(mateAlignmentStart);

        rec.setFirstOfPairFlag(firstOfPairFlag);
        rec.setReadNegativeStrandFlag(negativeStrandFlag);
        rec.setMateNegativeStrandFlag(mateNegativeStrandFlag);

        Assert.assertEquals(UmiUtil.isTopStrand(rec), topStrand);
    }

    @DataProvider(name = "brokenUmiProvider")
    private Object[][] testBrokenUmiDataProvider() {
        // The following are broken UMIs due to illegal characters that should result in a thrown exception when used.
        return new Object[][]{
                {"ATCxG-AGCG"},
                {"A1C"},
                {"@Agtc"},
                {"TA/AC"},
                {"AT:CG"}
        };
    }

    @Test(dataProvider = "brokenUmiProvider", expectedExceptions = PicardException.class)
    public void testBrokenUmi(final String brokenUmi) {
        SAMFileHeader header = new SAMFileHeader();
        SAMRecord rec = new SAMRecord(header);

        rec.setAttribute("RX", brokenUmi);

        // This should throw an exception due to a broken UMI in rec
        UmiUtil.getTopStrandNormalizedUmi(rec, "RX", true);
    }
}
