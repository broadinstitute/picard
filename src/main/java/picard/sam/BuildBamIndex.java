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

package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.nio.file.Path;

/**
 * Command line program to generate a BAM index (.bai) file from a BAM (.bam) file
 *
 * @author Martha Borkan
 */
@CommandLineProgramProperties(
        summary = BuildBamIndex.USAGE_SUMMARY + BuildBamIndex.USAGE_DETAILS,
        oneLineSummary = BuildBamIndex.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class BuildBamIndex extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Generates a BAM index \".bai\" file.  ";
    static final String USAGE_DETAILS = "This tool creates an index file for the input BAM that allows fast look-up of data in a " +
            "BAM file, lke an index on a database. Note that this tool cannot be run on SAM files, and that the input BAM file must be " +
            "sorted in coordinate order." +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar BuildBamIndex \\<br />" +
            "      I=input.bam" +
            "</pre>" +
            "<hr />";
    private static final Log log = Log.getInstance(BuildBamIndex.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "A BAM file or GA4GH URL to process. Must be sorted in coordinate order.")
    public PicardHtsPath INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai.\n" +
                    "If INPUT is a URL and OUTPUT is unspecified, defaults to a file in the current directory.", optional = true)
    public PicardHtsPath OUTPUT;

    // The legacy behavior is to use the ".bam.bai"; change it to ".bai"
    public boolean USE_SANE_EXTENSION = true;

    /**
     * Main method for the program.  Checks that all input files are present and
     * readable and that the output file can be written to.  Then iterates through
     * all the records generating a BAM Index, then writes the bai file.
     */
    protected int doWork() {
        final Path inputPath = INPUT.toPath();

        // set default output file - input-file.bai
        if (OUTPUT == null) {
            final String baseFileName = inputPath.toUri().toString();

            // only BAI indices can be created for now, although CSI indices can be read as well
            // tsato: would be good to update htsjdk BamFileIoUtils.isBamFile() to support Path
            // tsato: since we check whether the input is bam below, this check is unnecessary.
            // csi file...is that inherent to a bam file? Why would anyone input a csi index to buildbamindex....
            if (baseFileName.endsWith(FileExtensions.BAM)) {
                final int index = baseFileName.lastIndexOf('.');
                OUTPUT = new PicardHtsPath(baseFileName.substring(0, index) + FileExtensions.BAI_INDEX);
            } else {
                // tsato: why give two different kinds of output?
                // Given them a benefit of the doubt...https://github.com/broadinstitute/picard/pull/998/files
                OUTPUT = new PicardHtsPath(baseFileName + FileExtensions.BAI_INDEX);
            } // tsato: I guess the possible scenario is that there is already an index in the same directory, and we want to use
        } // the same kind (csi or bai). Or, this index type information is somehow already encoded in a bam file.

        final SamReader bam;
        bam = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE)
                .disable(SamReaderFactory.Option.EAGERLY_DECODE)
                .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                .open(SamInputResource.of(inputPath));

        if (bam.type() != SamReader.Type.BAM_TYPE && bam.type() != SamReader.Type.BAM_CSI_TYPE) { // PR: ? what did this PR do exactly?
            throw new SAMException("Input file must be bam file, not sam file.");
        } // tsato: see BamFileREader.getIndexType...

        if (!bam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
            throw new SAMException("Input bam file must be sorted by coordinate");
        }

        BAMIndexer.createIndex(bam, OUTPUT.toPath());

        log.info("Successfully wrote bam index file " + OUTPUT.toPath().toUri());
        CloserUtil.close(bam);
        return 0;
    }
}
