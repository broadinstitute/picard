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
package net.sf.samtools.example;

import net.sf.samtools.*;

import java.io.File;

public class ExampleSamUsage {
    /**
     * Read a SAM or BAM file, convert each read name to upper case, and write a new
     * SAM or BAM file.
     */
    public void convertReadNamesToUpperCase(final File inputSamOrBamFile, final File outputSamOrBamFile) {

        // Open the input file.  Automatically detects whether input is SAM or BAM
        // and delegates to a reader implementation for the appropriate format.
        final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);

        // makeSAMorBAMWriter() writes a file in SAM text or BAM binary format depending
        // on the file extension, which must be either .sam or .bam.

        // Since the SAMRecords will be written in the same order as they appear in the input file,
        // and the output file is specified as having the same sort order (as specified in
        // SAMFileHeader.getSortOrder(), presorted == true.  This is much more efficient than
        // presorted == false, if coordinate or queryname sorting is specified, because the SAMRecords
        // can be written to the output file directly rather than being written to a temporary file
        // and sorted after all records have been sent to outputSam.

        final SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),
                true, outputSamOrBamFile);

        for (final SAMRecord samRecord : inputSam) {
            // Convert read name to upper case.
            samRecord.setReadName(samRecord.getReadName().toUpperCase());
            outputSam.addAlignment(samRecord);
        }

        outputSam.close();
        inputSam.close();
    }
}
