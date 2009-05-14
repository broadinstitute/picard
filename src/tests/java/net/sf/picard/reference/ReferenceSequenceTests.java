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

package net.sf.picard.reference;

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Random;

/**
 * Tests for the reading of reference sequences in various formats.
 *
 * @author Tim Fennell
 */
public class ReferenceSequenceTests {
    private static final byte[] BASES = "acgtACGTN".getBytes();
    private final Random random = new Random();

    @Test(dataProvider="fastaTestParameters")
    public void testSingleShortSequence(int chroms, int basesPerChrom) throws Exception {
        File f = makeRandomReference(chroms, basesPerChrom);
        ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(f);

        for (int i=1; i<=chroms; ++i) {
            ReferenceSequence seq = ref.nextSequence();
            Assert.assertNotNull(seq);
            Assert.assertEquals(seq.length(), basesPerChrom);
            Assert.assertEquals(seq.getName(), "chr" + i);
            Assert.assertEquals(seq.getContigIndex(), i-1);
        }

        Assert.assertNull(ref.nextSequence());
    }

    @DataProvider
    Object[][] fastaTestParameters() {
        return new Object[][] {
                new Object[] { 1,     60},
                new Object[] { 2,     60},
                new Object[] {10,     60},
                new Object[] { 1,   1000},
                new Object[] { 2,   1000},
                new Object[] {10,   1000},
                new Object[] { 1, 250000},
                new Object[] { 2, 250000},
                new Object[] {10, 250000}
        };
    }


    /** Utility method to write a random reference sequence of specified length. */
    private File makeRandomReference(int chroms, int basesPerChrom) throws Exception {
        File file = File.createTempFile("reference.", ".fasta");
        file.deleteOnExit();
        BufferedWriter out = new BufferedWriter(new FileWriter(file));

        for (int i=1; i<=chroms; ++i) {
            out.write("> chr" + i);
            out.newLine();

            for (int j=1; j<=basesPerChrom; ++j) {
                out.write(BASES[random.nextInt(BASES.length)]);

                if (j % 80 == 0 || j == basesPerChrom) out.newLine();
            }
        }

        out.flush();
        out.close();
        return file;
    }
}
