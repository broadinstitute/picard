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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.argumentcollections.EmptyOutputArgumentCollection;
import picard.cmdline.argumentcollections.OutputArgumentCollection;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.AbstractMarkDuplicatesCommandLineProgram;
import picard.sam.markduplicates.util.LibraryIdGenerator;

import java.io.File;

/**
 *
 * Collect DuplicateMark'ing metrics from an input file that was already Duplicate-Marked.
 *
 */
@CommandLineProgramProperties(
        summary = CollectDuplicateMetrics.USAGE_SUMMARY + CollectDuplicateMetrics.USAGE_DETAILS,
        oneLineSummary = CollectDuplicateMetrics.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature

public class CollectDuplicateMetrics extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "Collect Duplicate metrics from marked file.";
    static final String USAGE_DETAILS = "\n" +
            "This tool only collects the duplicate metrics from a file that has already been duplicate-marked. " +
            "The resulting metrics file will always have a READ_PAIR_OPTICAL_DUPLICATES=0 and as a result the ESTIMATED_LIBRARY_SIZE will be slightly incorrect. ";

    @Argument(shortName = "M",
            doc = "File to write duplication metrics to.")
    public File METRICS_FILE;

    @Override
    protected OutputArgumentCollection getOutputArgumentCollection() {
        return new EmptyOutputArgumentCollection();
    }

    private LibraryIdGenerator libraryIdGenerator;
    private SAMFileHeader header;

    @Override
    protected boolean requiresReference() {
        return false;
    }

    @Override
    protected boolean usesNoRefReads() {
        return true;
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {

        libraryIdGenerator = new LibraryIdGenerator(header);
        this.header = header;
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {

        final DuplicationMetrics metrics = AbstractMarkDuplicatesCommandLineProgram.addReadToLibraryMetrics(rec, header, libraryIdGenerator, false);
        final boolean isDuplicate = rec.getDuplicateReadFlag();

        if (isDuplicate) {
            metrics.addDuplicateReadToMetrics(rec);
        }
    }

    @Override
    protected void finish() {
        // Write out the metrics
        AbstractMarkDuplicatesCommandLineProgram.finalizeAndWriteMetrics(libraryIdGenerator, getMetricsFile(), METRICS_FILE);
    }
}
