package picard.sam.util;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Tuple;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.EnumSet;
import java.util.Set;

/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

public class PrimaryAlignmentKeyTest {

    final SAMFileHeader samHeader = new SAMFileHeader();

    @DataProvider(name="positiveTestData")
    public Object[][] positiveTestData() {
        return new Object[][]{
                // unpaired < first
                { makeRecordPair("r1", Collections.EMPTY_SET,
                        "r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.FIRST_OF_PAIR)), -1 },
                // unpaired < second
                { makeRecordPair("r1", Collections.EMPTY_SET,
                        "r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.SECOND_OF_PAIR)), -2 },
                // first < second
                { makeRecordPair("r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.FIRST_OF_PAIR),
                        "r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.SECOND_OF_PAIR)), -1 },

                // first > unpaired
                { makeRecordPair("r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.FIRST_OF_PAIR),
                        "r1", Collections.EMPTY_SET), 1 },
                // second > unpaired
                { makeRecordPair("r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.SECOND_OF_PAIR),
                        "r1", Collections.EMPTY_SET), 2 },
                // second > first
                { makeRecordPair("r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.SECOND_OF_PAIR),
                        "r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.FIRST_OF_PAIR)), 1 },

                // unpaired == unpaired
                { makeRecordPair("r1", Collections.EMPTY_SET, "r1", Collections.EMPTY_SET), 0 },
                // first == first
                { makeRecordPair("r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.FIRST_OF_PAIR),
                        "r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.FIRST_OF_PAIR)), 0 },
                // second == second
                { makeRecordPair("r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.SECOND_OF_PAIR),
                        "r1", EnumSet.of(SAMFlag.READ_PAIRED, SAMFlag.SECOND_OF_PAIR)), 0 },

                // different read names, "r1" < "r2"
                { makeRecordPair("r1", Collections.EMPTY_SET, "r2", Collections.EMPTY_SET), -1 }
        };
    }

    @Test(dataProvider="positiveTestData")
    public void testPrimaryAlignmentKeyCompare(
            Tuple<SAMRecord, SAMRecord> recPair,
            final int expectedCompare) {
        Assert.assertEquals(new PrimaryAlignmentKey(recPair.a).compareTo(new PrimaryAlignmentKey(recPair.b)), expectedCompare);
    }

    @DataProvider(name="negativeTestData")
    public Object[][] negativeTestData() {
        return new Object[][]{
                { makeSAMRecord("r1", EnumSet.of(SAMFlag.NOT_PRIMARY_ALIGNMENT)) },  // secondary!
                { makeSAMRecord("r1", EnumSet.of(SAMFlag.SUPPLEMENTARY_ALIGNMENT)) }, //supplementary
        };
    }

    // reject secondary and supplementary
    @Test(dataProvider="negativeTestData", expectedExceptions = IllegalArgumentException.class)
    public void testRejectNonPrimary(final SAMRecord rec) {
        new PrimaryAlignmentKey(rec);
    }

    private Tuple<SAMRecord, SAMRecord> makeRecordPair(
            final String firstReadName,
            final Set<SAMFlag> firstReadFlags,
            final String secondReadName,
            final Set<SAMFlag> secondReadFlags) {
        return new Tuple(
                makeSAMRecord(firstReadName, firstReadFlags),
                makeSAMRecord(secondReadName, secondReadFlags)
        );
    }

    private SAMRecord makeSAMRecord(final String readName, final Set<SAMFlag> readFlags)
    {
        final SAMRecord rec = new SAMRecord(samHeader);
        rec.setReadName(readName);
        rec.setFlags(readFlags.stream().map(SAMFlag::intValue).reduce(0, (a, b) -> a | b));
        return rec;
    }

}
