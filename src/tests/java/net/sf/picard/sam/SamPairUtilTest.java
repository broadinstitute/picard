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

import net.sf.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class SamPairUtilTest {

    @Test(dataProvider = "testGetPairOrientation")
    public void testGetPairOrientation(final SAMRecord rec, final SamPairUtil.PairOrientation expectedOrientation)
    throws Exception {
        Assert.assertEquals(SamPairUtil.getPairOrientation(rec), expectedOrientation);
        final SAMRecord negatedInsertSize = (SAMRecord)rec.clone();
        negatedInsertSize.setInferredInsertSize(-rec.getInferredInsertSize());
        Assert.assertEquals(SamPairUtil.getPairOrientation(negatedInsertSize), expectedOrientation);

    }

    @DataProvider(name = "testGetPairOrientation")
    public Object[][] testGetPairOrientationDataProvider() throws Exception{
        final ArrayList<Object[]> ret = new ArrayList<Object[]>();

        SAMRecord rec = new SAMRecord(null);
        rec.setReadString("CTGTGCAGAGACAATACGGCTGGCCCGCACTGTGAGAAGTGCAGTGATGGGTACTATGGAGATTCAACTGCAGGCACCTCCTCCGATTGCCAACCCTGTC");
        rec.setBaseQualities(new byte[rec.getReadLength()]);
        rec.setReadPairedFlag(true);
        rec.setProperPairFlag(true);
        rec.setFirstOfPairFlag(true);
        rec.setReadNegativeStrandFlag(true);
        rec.setCigarString("3S97M");
        rec.setReferenceName("chr1");
        rec.setAlignmentStart(200);
        rec.setMateReferenceName("chr1");
        rec.setMateAlignmentStart(200);
        rec.setInferredInsertSize(-97);
        ret.add(new Object[]{rec, SamPairUtil.PairOrientation.FR});

        rec = new SAMRecord(null);
        rec.setReadString("GTGCAGAGACAATACGGCTGGGCCGCACTGTGAGAAGTGCAG");
        rec.setBaseQualities(new byte[rec.getReadLength()]);
        rec.setReadPairedFlag(true);
        rec.setProperPairFlag(true);
        rec.setSecondOfPairFlag(true);
        rec.setMateNegativeStrandFlag(true);
        rec.setCigarString("1S39M2S");
        rec.setReferenceName("chr1");
        rec.setAlignmentStart(200);
        rec.setMateReferenceName("chr1");
        rec.setMateAlignmentStart(200);
        rec.setInferredInsertSize(97);
        ret.add(new Object[]{rec, SamPairUtil.PairOrientation.FR});

        return ret.toArray(new Object[0][]);
    }
}
