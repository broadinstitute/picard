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
package net.sf.samtools.util;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMTag;
import net.sf.samtools.SAMTextHeaderCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author alecw@broadinstitute.org
 */
public class SequenceUtilTest {
    private static final String HEADER = "@HD\tVN:1.0\tSO:unsorted\n";
    private static final String SEQUENCE_NAME=
        "@SQ\tSN:phix174.seq\tLN:5386\tUR:/seq/references/PhiX174/v0/PhiX174.fasta\tAS:PhiX174\tM5:3332ed720ac7eaa9b3655c06f6b9e196";

    @Test
    public void testExactMatch() {
        final SAMSequenceDictionary sd1 = makeSequenceDictionary(5386, "/seq/references/PhiX174/v0/PhiX174.fasta",
                "3332ed720ac7eaa9b3655c06f6b9e196");
        final SAMSequenceDictionary sd2 = makeSequenceDictionary(5386, "/seq/references/PhiX174/v0/PhiX174.fasta",
                "3332ed720ac7eaa9b3655c06f6b9e196");
        SequenceUtil.assertSequenceDictionariesEqual(sd1, sd2);
    }

    @Test(expectedExceptions = SequenceUtil.SequenceListsDifferException.class)
    public void testMismatch() {
        final SAMSequenceDictionary sd1 = makeSequenceDictionary(5386, "/seq/references/PhiX174/v0/PhiX174.fasta",
                "3332ed720ac7eaa9b3655c06f6b9e196");
        final SAMSequenceDictionary sd2 = makeSequenceDictionary(5386, "/seq/references/PhiX174/v0/PhiX174.fasta",
                "deadbeef");
        SequenceUtil.assertSequenceDictionariesEqual(sd1, sd2);
        Assert.fail();
    }

    @Test
    public void testFileColonDifference() {
        final SAMSequenceDictionary sd1 = makeSequenceDictionary(5386, "/seq/references/PhiX174/v0/PhiX174.fasta",
                "3332ed720ac7eaa9b3655c06f6b9e196");
        final SAMSequenceDictionary sd2 = makeSequenceDictionary(5386, "file:/seq/references/PhiX174/v0/PhiX174.fasta",
                "3332ed720ac7eaa9b3655c06f6b9e196");
        SequenceUtil.assertSequenceDictionariesEqual(sd1, sd2);
    }

    @Test
    public void testURDifferent() {
        final SAMSequenceDictionary sd1 = makeSequenceDictionary(5386, "/seq/references/PhiX174/v0/PhiX174.fasta",
                "3332ed720ac7eaa9b3655c06f6b9e196");
        final SAMSequenceDictionary sd2 = makeSequenceDictionary(5386, "file:/seq/references/PhiX174/v1/PhiX174.fasta",
                "3332ed720ac7eaa9b3655c06f6b9e196");
        SequenceUtil.assertSequenceDictionariesEqual(sd1, sd2);
    }

    private SAMSequenceDictionary makeSequenceDictionary(final int length, final String ur, final String m5) {
        final String s = HEADER +
                String.format("@SQ\tSN:phix174.seq\tLN:%d\tUR:%s\tAS:PhiX174\tM5:%s\n", length, ur, m5);
        return new SAMTextHeaderCodec().decode(new StringLineReader(s), null).getSequenceDictionary();
    }

    @Test(dataProvider = "makeReferenceFromAlignment")
    public void testMakeReferenceFromAlignment(final String seq, final String cigar, final String md,
                                               boolean includeReferenceBasesForDeletions,
                                               final String expectedReference) {
        final SAMRecord rec = new SAMRecord(null);
        rec.setReadName("test");
        rec.setReadString(seq);
        rec.setCigarString(cigar);
        rec.setAttribute(SAMTag.MD.name(), md);
        final byte[] refBases = SequenceUtil.makeReferenceFromAlignment(rec, includeReferenceBasesForDeletions);
        Assert.assertEquals(StringUtil.bytesToString(refBases), expectedReference);
    }

    @DataProvider(name = "makeReferenceFromAlignment")
    public Object[][] testMakeReferenceFromAlignmentDataProvider() {
        return new Object[][] {
                {"ACGTACGTACGT", "12M2H", "4GAAA4", true, "ACGTGAAAACGT"},
                {"ACGTACGTACGT", "2H12M", "12", false, "ACGTACGTACGT"},
                {"ACGTACGTACGT", "4M4I4M2H", "8", false, "ACGT----ACGT"},
                {"ACGTACGTACGT", "2S4M2I4M2S", "8", false, "00GTAC--ACGT00"},      
                {"ACGTACGTACGT", "6M2D6M2H", "4GA^TT0TG4", true, "ACGTGATTTGACGT"},
                {"ACGTACGTACGT", "6M2D6M2H", "4GA^TT0TG4", false, "ACGTGATGACGT"},
        };
    }
}
