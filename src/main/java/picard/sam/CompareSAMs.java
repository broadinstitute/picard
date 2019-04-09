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
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
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
            " If the file headers are comparable, the tool Can perform either a naive comparison for which " +
            "each alignment must be identical, or a more sophisticated check of \"equivalence\", where mapping quality" +
            "0 reads are allowed to have different alignments, and duplicate marks are allowed to differ to account for " +
            "ambiguities in selecting the representative read of a duplicate set.  Results of comparison are summarised in " +
            " an output metrics file." +
            "<h3>Usage example:</h3>" +
            "<h4>CompareSAMs for exact matching:</h4>" +
            "<pre>" +
            "java -jar picard.jar CompareSAMs \\<br />" +
            "      I=file_1.bam \\<br />" +
            "      I=file_2.bam \\<br />" +
            "      O=comparison.tsv" +
            "</pre>" +
            "\n" +
            "<h4>CompareSAMs for \"equivalence\":</h4>" +
            "java -jar picard.jar CompareSAMs \\<br />" +
            "      I=file_1.bam \\<br />" +
            "      I=file_2.bam \\<br />" +
            "      MQ0_MATCH=true \\<br />" +
            "      LENIENT_DUP=true \\<br />" +
            "      O=comparison.tsv" +
            "</pre>" +
            "\n" +
            "<h4>CompareSAMs for \"equivalence\", allow name of species in header to differ:</h4>" +
            "java -jar picard.jar CompareSAMs \\<br />" +
            "      I=file_1.bam \\<br />" +
            "      I=file_2.bam \\<br />" +
            "      MQ0_MATCH=true \\<br />" +
            "      LENIENT_DUP=true \\<br />" +
            "      IGNORE_SPECIES=true \\<br />" +
            "      O=comparison.tsv" +
            "<hr />";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Exactly two input .sam or .bam files to compare to one another.", minElements = 2, maxElements = 2)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file to write comparison results to.")
    public File OUTPUT;

    @Argument(doc = "Allow species name field in header lines to differ.")
    boolean IGNORE_SPECIES;

    @Argument(doc = "Allow assembly name field in header lines to differ.")
    boolean IGNORE_ASSEMBLY;

    @Argument(doc = "Allow UR field in header lines to differ.")
    boolean IGNORE_UR;

    @Argument(doc = "Allow M5 field in header lines to differ.")
    boolean IGNORE_M5;

    @Argument(doc = "Allow Program Record lines in header to differ.")
    boolean IGNORE_PG;

    @Argument(doc = "Perform lenient checking of duplicate marks.  In this mode , will reduce the number of mismatches by allowing the choice of the representative read in each duplicate set" +
            " to differ bewteen the input files, as long as the duplicate sets agree.")
    boolean LENIENT_DUP;

    @Argument(doc = "Count reads which have mapping quality below 3 in both files but are mapped to different locations as matches.  By default we count such reads as mismatching.")
    boolean MQ0_MATCH;

    /**
     * Do the work after command line has been parsed. RuntimeException may be
     * thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    @Override
    protected int doWork() {
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);

        try (final SamReader samReader1 = samReaderFactory.open(INPUT.get(0));
             final SamReader samReader2 = samReaderFactory.open(INPUT.get(1)))
        {
            final SamComparison comparison = new SamComparison(samReader1, samReader2, IGNORE_SPECIES, IGNORE_ASSEMBLY,
                    IGNORE_UR, IGNORE_M5, IGNORE_PG, LENIENT_DUP, MQ0_MATCH, INPUT.get(0).getAbsolutePath(), INPUT.get(1).getAbsolutePath());
            comparison.writeReport(OUTPUT, getDefaultHeaders());
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
