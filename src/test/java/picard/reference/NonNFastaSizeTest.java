/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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
package picard.reference;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.lang.Exception;
import java.lang.String;

/**
 * @author ebanks
 */

public class NonNFastaSizeTest {

	private static final String REFERENCE = "testdata/picard/reference/test.fasta";

	@Test
	public void noIntervals() throws IOException {
        final File input = new File(REFERENCE);
        final File outfile   = File.createTempFile("nonNcount", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath()
        };
        Assert.assertEquals(new NonNFastaSize().instanceMain(args), 0);

        final BufferedReader reader = IOUtil.openFileForBufferedReading(outfile);
        final String count = reader.readLine();

        try {
            Assert.assertEquals(Long.parseLong(count), 1008);
        } catch (Exception e) {
            System.err.println("Failed to read in count because of error: " + e.getMessage());
        }
    }

    @Test
    public void withIntervals() throws IOException {
        final File input = new File(REFERENCE);
        final File outfile = File.createTempFile("nonNcount", ".txt");
        final File intervals = new File("testdata/picard/reference/test.intervals");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "INTERVALS=" + intervals.getAbsolutePath()
        };
        Assert.assertEquals(new NonNFastaSize().instanceMain(args), 0);

        final BufferedReader reader = IOUtil.openFileForBufferedReading(outfile);
        final String count = reader.readLine();

        try {
            Assert.assertEquals(Long.parseLong(count), 53);
        } catch (Exception e) {
            System.err.println("Failed to read in count because of error: " + e.getMessage());
        }
    }
}
