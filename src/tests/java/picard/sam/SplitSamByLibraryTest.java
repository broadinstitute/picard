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
package picard.sam;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Iterator;

/**
 * Basic positive and negative tests for SplitSamByLibrary command-line program
 *
 * @author ktibbett@broadinstitute.org
 */
public class SplitSamByLibraryTest {

    @Test
    public void testNoLibrarySpecified() {
        SplitSamByLibrary splitter = new SplitSamByLibrary();
        splitter.INPUT = new File("testdata/picard/sam/invalid_coord_sort_order.sam");
        Assert.assertEquals(splitter.doWork(), SplitSamByLibrary.NO_LIBRARIES_SPECIFIED_IN_HEADER,
                "SAM file with no libraries should failed but didn't.");

    }

    @Test
    public void basicPositiveTest() {
        SplitSamByLibrary splitter = new SplitSamByLibrary();
        splitter.INPUT = new File("testdata/picard/sam/split_test.sam");
        Assert.assertEquals(splitter.doWork(), 0, "SAM file split should have succeeded but didn't.");

        File f = new File("unknown.sam");
        Assert.assertTrue(f.exists(), "uknown.sam should exist but doesn't");
        Assert.assertEquals(countReads(f), 2, "unknown.sam has the wrong number of reads");
        f.delete();

        f = new File("lib-1.sam");
        Assert.assertTrue(f.exists(), "lib-1.sam should exist but doesn't");
        Assert.assertEquals(countReads(f), 6, "lib-1.sam has the wrong number of reads");
        f.delete();

        f = new File("lib-2.sam");
        Assert.assertFalse(f.exists(), "lib-2.sam should not exist but does");
        if (f.exists()) f.delete();

        f = new File("lib-3.sam");
        Assert.assertTrue(f.exists(), "lib-3.sam should exist but doesn't");
        Assert.assertEquals(countReads(f), 2, "lib-3.sam has the wrong number of reads");
        f.delete();

    }

    @Test
    public void testNoUnknownFile() {
        SplitSamByLibrary splitter = new SplitSamByLibrary();
        splitter.INPUT = new File("testdata/picard/sam/split_test2.sam");
        Assert.assertEquals(splitter.doWork(), 0, "SAM file split should have succeeded but didn't.");

        // The unknown file should exist and have two reads
        File f = new File("unknown.sam");
        Assert.assertFalse(f.exists(), "uknown.sam should not exist but does");
        if (f.exists()) f.delete();

        f = new File("lib-1.sam");
        Assert.assertTrue(f.exists(), "lib-1.sam should exist but doesn't");
        Assert.assertEquals(countReads(f), 4, "lib-1.sam has the wrong number of reads");
        f.delete();

        f = new File("lib-3.sam");
        Assert.assertTrue(f.exists(), "lib-3.sam should exist but doesn't");
        Assert.assertEquals(countReads(f), 2, "lib-3.sam has the wrong number of reads");
        f.delete();

    }

    private int countReads(File samFile) {
        SamReader reader = SamReaderFactory.makeDefault().open(samFile);
        int count = 0;
        for (Iterator it = reader.iterator(); it.hasNext(); ) {
            it.next();
            count++;
        }
        CloserUtil.close(reader);
        return count;

    }
}
