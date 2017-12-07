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

package picard.sam.markduplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.util.RuntimeIOException;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.io.IOException;

public class ParallelMarkDuplicatesTester extends AbstractMarkDuplicatesCommandLineProgramTester {

    @Override
    protected CommandLineProgram getProgram() {
        return new ParallelMarkDuplicates();
    }

    @Override
    protected File createInputFile() {
        // Create the input file
        final File input = new File(outputDir, "input.bam");
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(samRecordSetBuilder.getHeader(), true, input);
        samRecordSetBuilder.getRecords().forEach(writer::addAlignment);
        writer.close();

        File index = new File(outputDir, "input.bam.bai");
        try (SamReader reader = SamReaderFactory.makeDefault().referenceSequence(Defaults.REFERENCE_FASTA)
                .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                .open(input)) {
            BAMIndexer.createIndex(
                    reader,
                    index);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
        return input;
    }
}
