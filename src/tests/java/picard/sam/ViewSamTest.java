/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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

import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.LineReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;

public class ViewSamTest {
    /**
     * Confirm that ViewSam retains whatever version number was in the input header.
     */
    @Test
    public void testHeaderVersion() throws Exception {
        final String oldVersionHeader = "@HD\tVN:1.3\tSO:unsorted";
        final File inputSam = File.createTempFile("ViewSamTest.input.", ".sam");
        inputSam.deleteOnExit();
        final AsciiWriter writer = new AsciiWriter(new FileOutputStream(inputSam));
        writer.write(oldVersionHeader);
        writer.write("\n");
        writer.close();
        final File viewSamOutputFile = File.createTempFile("ViewSamTest.output.", ".sam");
        viewSamOutputFile.deleteOnExit();
        final ViewSam viewSam = new ViewSam();
        viewSam.INPUT = inputSam.getAbsolutePath();
        final PrintStream viewSamPrintStream = new PrintStream(viewSamOutputFile);
        // TODO - Should switch over to using invocation via new PicardCommandLine() - BUT the test here is accessing class members directly.
        Assert.assertEquals(viewSam.writeSamText(viewSamPrintStream), 0);
        viewSamPrintStream.close();
        final LineReader viewSamInputReader = new BufferedLineReader(new FileInputStream(viewSamOutputFile));
        Assert.assertEquals(viewSamInputReader.readLine(), oldVersionHeader);
    }
}
