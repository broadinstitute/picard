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
package net.sf.picard.sam;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * A command-line tool to merge BAM/SAM alignment info from a third-party aligner with the data in an
 * unmapped BAM file, producing a third BAM file that has alignment data and all the additional data
 * from the unmapped BAM
 *
 * @author ktibbett@broadinstitute.org
 */
public class MergeBamAlignment extends CommandLineProgram {
    // Usage and parameters

    @Usage
    public String USAGE = getStandardUsagePreamble() +  "Merges alignment data from a SAM or BAM " +
            "file with additional data stored in an unmapped BAM file and produces a third SAM " +
            "or BAM file of aligned and unaligned reads.  NOTE that this program expects to " +
            "find a sequence dictionary in the same directory as REFERENCE_SEQUENCE and expects it " +
            "to have the same base name as the reference fasta except with the extension '.dict'";
    @Option(shortName="UNMAPPED", doc="Original SAM or BAM file of unmapped reads, which must " +
            "be in queryname order.")
        public File UNMAPPED_BAM;
    @Option(shortName="ALIGNED", doc="SAM or BAM file with alignment data.")
        public File ALIGNED_BAM;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Merged SAM or BAM file " +
        "to write to.") public File OUTPUT;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Path to the fasta " +
        "file for the reference sequence.") public File REFERENCE_SEQUENCE;
    @Option(shortName=StandardOptionDefinitions.PROGRAM_RECORD_ID_SHORT_NAME,
            doc="The program group ID of the aligner (if not supplied by the aligned file).",
            optional=true) public String PROGRAM_RECORD_ID;
    @Option(shortName="PG_VERSION", doc="The version of the program group (if not supplied by " +
            "the aligned file).", optional=true) public String PROGRAM_GROUP_VERSION;
    @Option(shortName="PG_COMMAND", doc="The command line of the program group (if not supplied " +
            "by the aligned file).", optional=true) public String PROGRAM_GROUP_COMMAND_LINE;
    @Option(shortName="PG_NAME", doc="The name of the program group (if not supplied by " +
            "the aligned file).", optional=true)public String PROGRAM_GROUP_NAME;
    @Option(doc="Whether this is a paired-end run. ", shortName="PE") public Boolean PAIRED_RUN;
    @Option(doc="The expected jump size (required if this is a jumping library).", shortName="JUMP",
            optional=true) public Integer JUMP_SIZE;
    @Option(doc="Whether to clip adapters where identified.") public boolean CLIP_ADAPTERS = true;
    @Option(doc="Whether the lane is bisulfite sequence (used when caculating the NM tag).")
        public boolean IS_BISULFITE_SEQUENCE = false;
    @Option(doc="Whether to output only aligned reads.  ") public boolean ALIGNED_READS_ONLY = false;
    @Option(doc="The maximum number of insertions or deletions permitted for an alignment to be " +
            "included.  Alignments with more than this many insertions or deletions will be ignored.",
            shortName="MAX_GAPS") public int MAX_INSERTIONS_OR_DELETIONS = 1;
    @Option(doc="Reserved alignment attributes (tags starting with X, Y, or Z) that should be " +
            "brought over from the alignment data when merging.")
    public List<String> ATTRIBUTES_TO_RETAIN = new ArrayList<String>();


    private static final Log log = Log.getInstance(MergeBamAlignment.class);

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new MergeBamAlignment().instanceMain(argv));
    }

    @Override
    protected int doWork() {

        // Check the files are readable/writable
        SAMProgramRecord prod = null;
        if (PROGRAM_RECORD_ID != null) {
            prod = new SAMProgramRecord(PROGRAM_RECORD_ID);
            prod.setProgramVersion(PROGRAM_GROUP_VERSION);
            prod.setCommandLine(PROGRAM_GROUP_COMMAND_LINE);
            prod.setProgramName(PROGRAM_GROUP_NAME);
        }
        SamAlignmentMerger merger = new SamAlignmentMerger (UNMAPPED_BAM, OUTPUT,
            REFERENCE_SEQUENCE, prod, CLIP_ADAPTERS, IS_BISULFITE_SEQUENCE, PAIRED_RUN,
            JUMP_SIZE != null, ALIGNED_READS_ONLY, ALIGNED_BAM,
            MAX_INSERTIONS_OR_DELETIONS, ATTRIBUTES_TO_RETAIN);
        merger.mergeAlignment();
        return 0;
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     * @return null if command line is valid.  If command line is invalid, returns
     *         an array of error messages to be written to the appropriate place.
     */
    protected String[] customCommandLineValidation() {
        if ((PROGRAM_RECORD_ID != null || PROGRAM_GROUP_VERSION != null ||
             PROGRAM_GROUP_COMMAND_LINE != null) &&
            (PROGRAM_RECORD_ID == null || PROGRAM_GROUP_VERSION == null ||
             PROGRAM_GROUP_COMMAND_LINE == null)) {

            return new String[] {"PROGRAM_RECORD_ID, PROGRAM_GROUP_VERSION, and " +
                    "PROGRAM_GROUP_COMMAND_LINE must all be supplied or none should " +
                    "be included."};
        }
        return null;
    }

}
