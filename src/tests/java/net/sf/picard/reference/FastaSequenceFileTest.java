/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

import net.sf.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * @author alecw@broadinstitute.org
 */
public class FastaSequenceFileTest {
    @Test
    public void testTrailingWhitespace() throws Exception {
        final File fasta = File.createTempFile("test", ".fasta");
        fasta.deleteOnExit();
        final PrintWriter writer = new PrintWriter(fasta);
        final String chr1 = "chr1";
        writer.println(">" + chr1);
        final String sequence = "ACGTACGT";
        writer.println(sequence);
        writer.println(sequence + " \t");
        writer.close();
        final FastaSequenceFile fastaReader = new FastaSequenceFile(fasta, true);
        final ReferenceSequence referenceSequence = fastaReader.nextSequence();
        Assert.assertEquals(referenceSequence.getName(), chr1);
        Assert.assertEquals(StringUtil.bytesToString(referenceSequence.getBases()), sequence + sequence);
    }

    @Test
    public void testIntermediateWhitespace() throws Exception {
        final File fasta = File.createTempFile("test", ".fasta");
        fasta.deleteOnExit();
        final PrintWriter writer = new PrintWriter(fasta);
        final String chr1 = "chr1";
        writer.println(">" + chr1 + " extra stuff after sequence name");
        final String sequence = "ACGTACGT";
        writer.println(sequence + "  ");
        writer.println(sequence + " \t");
        writer.println(sequence);
        writer.close();
        final FastaSequenceFile fastaReader = new FastaSequenceFile(fasta, true);
        final ReferenceSequence referenceSequence = fastaReader.nextSequence();
        Assert.assertEquals(referenceSequence.getName(), chr1);
        Assert.assertEquals(StringUtil.bytesToString(referenceSequence.getBases()), sequence + sequence + sequence);
    }
}
