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
package net.sf.samtools;

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

import java.io.File;

import net.sf.samtools.util.CloseableIterator;

/**
 * @author alecw@broadinstitute.org
 */
public class BAMIteratorTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/samtools");

    @Test(dataProvider = "dataProvider")
    public void testIterateEmptyBam(final String bam) throws Exception {
        final SAMFileReader reader = new SAMFileReader(new File(TEST_DATA_DIR, bam));
        int numRecords = 0;
        for (final SAMRecord rec : reader) {
            ++numRecords;
        }
        Assert.assertEquals(numRecords, 0);
    }

    @Test(dataProvider = "dataProvider")
    public void testQueryUnmappedEmptyBam(final String bam) throws Exception {
        SAMFileReader reader = new SAMFileReader(new File(TEST_DATA_DIR, bam));
        CloseableIterator<SAMRecord> it = reader.queryUnmapped();
        int numRecords = 0;
        while (it.hasNext()) {
            it.next();
            ++numRecords;
        }
        Assert.assertEquals(numRecords, 0);
    }

    @DataProvider(name = "dataProvider")
    public Object[][] bams() {
        return new Object[][] {
                {"empty.bam"},
                {"empty_no_empty_gzip_block.bam"}
        };
    }
}
