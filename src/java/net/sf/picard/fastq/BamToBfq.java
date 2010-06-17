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

package net.sf.picard.fastq;

import net.sf.picard.fastq.BamToBfqWriter;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;

import java.io.File;

/**
 * CommandLineProgram to generate to bfq files for use by the Maq aligner
 *
 * @author ktibbett@broadinstitute.org
 */
public class BamToBfq extends CommandLineProgram {

    private static final String PROGRAM_VERSION = "1.0";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Create BFQ files for use by the Maq aligner.";

    @Option(doc="The BAM file to parse.", shortName=StandardOptionDefinitions.INPUT_SHORT_NAME) public File INPUT;
    @Option(doc="The analysis directory for the binary output file. ") public File ANALYSIS_DIR;
    @Option(doc="Flowcell barcode (e.g. 30PYMAAXX).  ", shortName="F") public String FLOWCELL_BARCODE;
    @Option(doc="Lane number. ", shortName= StandardOptionDefinitions.LANE_SHORT_NAME, optional=true) public Integer LANE;
    @Option(doc="Number of reads to align (null = all).", shortName="NUM", optional=true) public Integer READS_TO_ALIGN;
    @Option(doc="Number of reads to break into individual groups for alignment", shortName="CHUNK") public Integer READ_CHUNK_SIZE = 2000000;
    @Option(doc="Whether this is a paired-end run. ", shortName="PE") public Boolean PAIRED_RUN;
    @Option(doc="Deprecated option; use READ_NAME_PREFIX instead", mutex="READ_NAME_PREFIX", shortName="RB", optional=true) public String RUN_BARCODE;
    @Option(doc="Prefix to be stripped off the beginning of all read names  (to make them short enough to run in Maq)", optional=true) public String READ_NAME_PREFIX;
    @Option(doc="Whether to include non-PF reads", shortName="NONPF", optional=true) public Boolean INCLUDE_NON_PF_READS = false;
    @Option(doc="Whether to clip adapters from the reads") public boolean CLIP_ADAPTERS = true;

    protected int doWork() {

        String outputPrefix = ANALYSIS_DIR.getAbsolutePath();
        if (!outputPrefix.endsWith("/")) {
            outputPrefix += "/";
        }
        outputPrefix += FLOWCELL_BARCODE + "." + LANE + ".";

        BamToBfqWriter writer = new BamToBfqWriter(INPUT, outputPrefix, READS_TO_ALIGN,
                READ_CHUNK_SIZE, PAIRED_RUN, RUN_BARCODE != null ? RUN_BARCODE + ":" : READ_NAME_PREFIX,
                INCLUDE_NON_PF_READS, CLIP_ADAPTERS);
        writer.writeBfqFiles();
        return 0;
    }

    public static void main(String[] argv) {
        System.exit(new BamToBfq().instanceMain(argv));
    }

}