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

package picard.fastq;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;

/**
 * Converts a BAM file into a BFQ (binary fastq formatted) file.
 * <p>
 * The BFQ format is the input format to some tools like Maq aligner.
 * </p>
 * <h3>Input</h3>
 * <p>A single BAM file to convert</p>.
 * <h3>Output</h3>
 * <p>One or two FASTQ files depending on whether the BAM file contains single- or
 * paired-end sequencing data. You must indicate the output directory that will contain these files (<code>ANALYSIS_DIR</code>)
 * and the output file name prefix (<code>OUTPUT_FILE_PREFIX</code>).</p>
 * <h3>Usage example:</h3>
 * <pre>
 *     java -jar picard.jar BamToBfq \
 *             I=input.bam \
 *             ANALYSIS_DIR=output_dir \
 *             OUTPUT_FILE_PREFIX=output_name \
 *             PAIRED_RUN=false
 * </pre>
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = "<p>" + BamToBfq.USAGE_SUMMARY + ".</p>" + BamToBfq.USAGE_DETAILS,
        oneLineSummary = BamToBfq.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class BamToBfq extends CommandLineProgram {
    static final String USAGE_SUMMARY =
            "Converts a BAM file into a BFQ (binary fastq formatted) file";
    static final String USAGE_DETAILS =
            "<p>The BFQ format is the input format to some tools like Maq aligner.</p>" +
            "<h3>Input</h3>" +
            "<p>A single BAM file to convert</p>" +
            "<h3>Output</h3>" +
            "<p>One or two FASTQ files depending on whether the BAM file contains single- or " +
            "paired-end sequencing data. You must indicate the output directory that will contain these files (<code>ANALYSIS_DIR</code>) " +
            "and the output file name prefix (<code>OUTPUT_FILE_PREFIX</code>).</p>" +
            "<h3>Usage example:</h3>" +
            "<pre>" +
            "java -jar picard.jar BamToBfq \\\n" +
            "     I=input.bam \\\n" +
            "     ANALYSIS_DIR=output_dir \\\n" +
            "     OUTPUT_FILE_PREFIX=output_name \\\n" +
            "     PAIRED_RUN=false" +
            "</pre><hr />";
    // The following attributes define the command-line arguments

    @Argument(doc="The BAM file to parse.", shortName=StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc="The analysis directory for the binary output file. ")
    public File ANALYSIS_DIR;

    @Argument(doc="Flowcell barcode (e.g. 30PYMAAXX).  ", shortName="F", mutex="OUTPUT_FILE_PREFIX")
    public String FLOWCELL_BARCODE;

    @Argument(doc="Lane number. ", shortName= StandardOptionDefinitions.LANE_SHORT_NAME, optional=true,mutex="OUTPUT_FILE_PREFIX")
    public Integer LANE;

    @Argument(doc="Prefix for all output files", mutex={"FLOWCELL_BARCODE","LANE"})
    public String OUTPUT_FILE_PREFIX;

    @Argument(doc="Number of reads to align (null = all).", shortName="NUM", optional=true)
    public Integer READS_TO_ALIGN;

    @Argument(doc="Number of reads to break into individual groups for alignment", shortName="CHUNK")
    public Integer READ_CHUNK_SIZE = 2000000;

    @Argument(doc="Whether this is a paired-end run. ", shortName="PE")
    public Boolean PAIRED_RUN;

    @Argument(doc="Deprecated option; use READ_NAME_PREFIX instead", mutex="READ_NAME_PREFIX", shortName="RB", optional=true)
    public String RUN_BARCODE;

    @Argument(doc="Prefix to be stripped off the beginning of all read names  (to make them short enough to run in Maq)", optional=true)
    public String READ_NAME_PREFIX;

    @Argument(doc="Whether to include non-PF reads", shortName="NONPF", optional=true)
    public Boolean INCLUDE_NON_PF_READS = false;

    @Argument(doc="Whether to clip adapters from the reads")
    public boolean CLIP_ADAPTERS = true;

    @Argument(doc="The number of bases from each read to write to the bfq file.  If this is non-null, then " +
            "only the first BASES_TO_WRITE bases from each read will be written.", optional=true)
    public Integer BASES_TO_WRITE = null;

    protected int doWork() {

        String outputPrefix = ANALYSIS_DIR.getAbsolutePath();
        if (!outputPrefix.endsWith("/")) {
            outputPrefix += "/";
        }
        outputPrefix += OUTPUT_FILE_PREFIX + ".";

        SamToBfqWriter writer = new SamToBfqWriter(INPUT, REFERENCE_SEQUENCE, outputPrefix, READS_TO_ALIGN,
                READ_CHUNK_SIZE, PAIRED_RUN, READ_NAME_PREFIX,
                INCLUDE_NON_PF_READS, CLIP_ADAPTERS, BASES_TO_WRITE);
        writer.writeBfqFiles();
        return 0;
    }

    protected String[] customCommandLineValidation() {

        if (OUTPUT_FILE_PREFIX == null) {
            OUTPUT_FILE_PREFIX = FLOWCELL_BARCODE + "." + LANE;
        }
        if (READ_NAME_PREFIX == null) {
            READ_NAME_PREFIX = RUN_BARCODE + ":";
        }
        return null;
    }
}
