/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class SamFileConverterTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/SamFileConverterTest");
    private static final File unmappedSam = new File(TEST_DATA_DIR, "unmapped.sam");
    private static final File unmappedBam = new File(TEST_DATA_DIR, "unmapped.bam");
    private static final File unmappedCram = new File(TEST_DATA_DIR, "unmapped.cram");

    @Test
    public void testSAMToBAM() {
        convertFile(unmappedSam, unmappedBam, ".bam");
    }

    @Test
    public void testSAMToCRAM() {
        convertFile(unmappedSam, unmappedCram, ".cram");
    }

    @Test
    public void testBAMToCRAM() {
        convertFile(unmappedBam, unmappedCram, ".cram");
    }

    @Test
    public void testBAMToSAM() {
        convertFile(unmappedBam, unmappedSam, ".sam");
    }

    @Test
    public void testCRAMToBAM() {
        convertFile(unmappedCram, unmappedBam, ".bam");
    }

    @Test
    public void testCRAMToSAM() {
        convertFile(unmappedCram, unmappedSam, ".sam");
    }


    private void convertFile(final File inputFile, final File fileToCompare, final String extension) {
        final SamFormatConverter samFormatConverter = new SamFormatConverter();
        final List<File> samFiles = new ArrayList<File>();
        final ValidateSamFile validateSamFile = new ValidateSamFile();
        final CompareSAMs compareSAMs = new CompareSAMs();

        samFormatConverter.INPUT = inputFile;
        try {
            samFormatConverter.OUTPUT = File.createTempFile("SamFileConverterTest." + inputFile.getName(), extension);
        } catch (final IOException e) {
            e.printStackTrace();
        }
        samFormatConverter.doWork();

        validateSamFile.INPUT = samFormatConverter.OUTPUT;
        assertEquals(validateSamFile.doWork(), 0);

        samFiles.add(samFormatConverter.OUTPUT);
        samFiles.add(fileToCompare);

        compareSAMs.samFiles = samFiles;
        compareSAMs.doWork();

        assertTrue(compareSAMs.areEqual());
    }
}
