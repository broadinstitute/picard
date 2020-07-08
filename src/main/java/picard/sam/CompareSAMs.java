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
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.PositionalArguments;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.sam.util.SAMComparisonArgumentCollection;
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
    static final String USAGE_SUMMARY = "Compare two input SAM/BAM/CRAM files.  ";
    static final String USAGE_DETAILS = "This tool initially compares the headers input files. " +
            " If the file headers are comparable, the tool can perform either strict comparisons for which " +
            "each alignment and the header must be identical, or a more lenient check of \"equivalence\", where reads with mapping quality < LOW_MQ_THRESHOLD " +
            "are allowed to have different alignments, duplicate marks are allowed to differ to account for " +
            "ambiguities in selecting the representative read of a duplicate set, and some differences in headers is allowed.  By default, alignment comparisons, " +
            "duplicate marking comparisons, and header comparisons are performed in the strict mode.  Results of comparison are summarised in " +
            " an output metrics file." +
            "<h3>Usage example:</h3>" +
            "<h4>CompareSAMs for exact matching:</h4>" +
            "<pre>" +
            "java -jar picard.jar CompareSAMs \\<br />" +
            "      file_1.bam \\<br />" +
            "      file_2.bam \\<br />" +
            "      O=comparison.tsv" +
            "</pre>" +
            "\n" +
            "<h4>CompareSAMs for \"equivalence\":</h4>" +
            "java -jar picard.jar CompareSAMs \\<br />" +
            "      file_1.bam \\<br />" +
            "      file_2.bam \\<br />" +
            "      LENIENT_LOW_MQ_ALIGNMENT=true \\<br />" +
            "      LENIENT_DUP=true \\<br />" +
            "      O=comparison.tsv" +
            "</pre>" +
            "\n" +
            "<h4>CompareSAMs for \"equivalence\", allow certain differences in header:</h4>" +
            "java -jar picard.jar CompareSAMs \\<br />" +
            "      file_1.bam \\<br />" +
            "      file_2.bam \\<br />" +
            "      LENIENT_LOW_MQ_ALIGNMENT=true \\<br />" +
            "      LENIENT_DUP=true \\<br />" +
            "      LENIENT_HEADER=true \\<br />" +
            "      O=comparison.tsv" +
            "<hr />";
    @PositionalArguments(doc = "Exactly two input SAM/BAM/CRAM files to compare to one another.", minElements = 2, maxElements = 2)
    public List<File> SAM_FILES;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file to write comparison results to.", optional = true)
    public File OUTPUT;

    @ArgumentCollection
    public SAMComparisonArgumentCollection samComparisonArgumentCollection = new SAMComparisonArgumentCollection();



    /**
     * Do the work after command line has been parsed. RuntimeException may be
     * thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    @Override
    protected int doWork() {
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);

        try (final SamReader samReader1 = samReaderFactory.open(SAM_FILES.get(0));
             final SamReader samReader2 = samReaderFactory.open(SAM_FILES.get(1)))
        {
            final SamComparison comparison = new SamComparison(samReader1, samReader2,
                    SAM_FILES.get(0).getAbsolutePath(), SAM_FILES.get(1).getAbsolutePath(), samComparisonArgumentCollection);
            if (OUTPUT != null) {
                comparison.writeReport(OUTPUT, getDefaultHeaders());
            }
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
