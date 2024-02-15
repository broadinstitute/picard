/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

import htsjdk.samtools.DuplicateScoringStrategy;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.io.IOException;

public class CollectDuplicateMetricsTester extends AbstractMarkDuplicatesCommandLineProgramTester {

    @Override
    public File getOutput() {
        final File output;
        try {
            output = File.createTempFile("CollectDuplicateMetrics", ".sam");
        } catch (IOException e) {
            throw new PicardException("problems creating output file", e);
        }
        output.deleteOnExit();

        try(SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(samRecordSetBuilder.getHeader(), true, output.toPath(), this.fastaFiles.get(samRecordSetBuilder.getHeader()))) {
            samRecordSetBuilder.getRecords().forEach(a -> {
                final SAMRecord b = a.deepCopy();
                b.setDuplicateReadFlag(duplicateFlags.get(samRecordToDuplicatesFlagsKey(a)));
                writer.addAlignment(b);
            });
        }

        return output;
    }

    @Override
    protected void customPretestCallback() {
        super.customPretestCallback();

        getArgs().removeIf(s->s.startsWith("OUTPUT="));
        getArgs().removeIf(s->s.startsWith("READ_NAME_REGEX="));
        getArgs().removeIf(s->s.startsWith("DUPLICATE_SCORING_STRATEGY="));
        getArgs().removeIf(s->s.startsWith("TAGGING_POLICY="));
        getArgs().removeIf(s->s.startsWith("ASSUME_SORT_ORDER="));

        getArgs().removeIf(s->s.startsWith("INPUT="));
        getArgs().add("INPUT=" + getOutput().getAbsoluteFile());
    }

    @Override
    protected CommandLineProgram getProgram() {
        return new CollectDuplicateMetrics();
    }

    public CollectDuplicateMetricsTester( final SAMFileHeader.SortOrder sortOrder) {
        super(DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES, sortOrder, false);
    }

    @Override
    public void recordOpticalDuplicatesMarked() {}

    @Override
    void updateExpectationsHook() {
        // CollectDuplicateMetrics cannot identify Optical duplicates....
        setExpectedOpticalDuplicate(0);
    }
}
