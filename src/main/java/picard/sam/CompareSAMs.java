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
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.PositionalArguments;
import picard.sam.util.SamComparison;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Rudimentary SAM comparer.  Compares headers, and if headers are compatible enough, compares SAMRecords,
 * looking only at basic alignment info.  Summarizes the number of alignments that match, mismatch, are missing, etc.
 *
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = CompareSAMs.USAGE_SUMMARY + CompareSAMs.USAGE_DETAILS,
        oneLineSummary = CompareSAMs.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class)
@DocumentedFeature
public class CompareSAMs extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Compare two input \".sam\" or \".bam\" files.  ";
    static final String USAGE_DETAILS = "This tool initially compares the headers of SAM or BAM files. " +
            " If the file headers are comparable, the tool will examine and compare readUnmapped flag, reference name, " +
            "start position and strand between the SAMRecords. The tool summarizes information on the number of read " +
            "pairs that match or mismatch, and of reads that are missing or unmapped (stratified by direction: " +
            "forward or reverse)." +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CompareSAMs \\<br />" +
            "      file_1.bam \\<br />" +
            "      file_2.bam" +
            "</pre>" +
            "<hr />";
    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> samFiles;

    /**
     * Do the work after command line has been parsed. RuntimeException may be
     * thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    @Override
    protected int doWork() {
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);

        try (final SamReader samReader1 = samReaderFactory.open(samFiles.get(0));
             final SamReader samReader2 = samReaderFactory.open(samFiles.get(1)))
        {
            final SamComparison comparison = new SamComparison(samReader1, samReader2);
            comparison.printReport();
            if (comparison.areEqual()) {
                System.out.println("SAM files match.");
            } else {
                System.out.println("SAM files differ.");
            }
            return comparison.areEqual() ? 0 : 1;
        } catch (IOException e) {
            throw new RuntimeIOException("Error opening file", e);
        }
    }

}
