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
import net.sf.samtools.SAMFileHeader.SortOrder;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
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
    @Option(shortName="ALIGNED", doc="SAM or BAM file(s) with alignment data.", mutex={"READ1_ALIGNED_BAM","READ2_ALIGNED_BAM"}, optional=true)
        public List<File> ALIGNED_BAM;
    @Option(shortName="R1_ALIGNED", doc="SAM or BAM file(s) with alignment data from the first read of a pair.", mutex={"ALIGNED_BAM"}, optional=true)
        public List<File> READ1_ALIGNED_BAM;
    @Option(shortName="R2_ALIGNED", doc="SAM or BAM file(s) with alignment data from the second read of a pair.", mutex={"ALIGNED_BAM"}, optional=true)
        public List<File> READ2_ALIGNED_BAM;
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
    @Option(doc="The expected jump size (required if this is a jumping library).  Deprecated.  " +
            "Use EXPECTED_ORIENTATIONS instead", shortName="JUMP", mutex="EXPECTED_ORIENTATIONS",
            optional=true) public Integer JUMP_SIZE;
    @Option(doc="Whether to clip adapters where identified.") public boolean CLIP_ADAPTERS = true;
    @Option(doc="Whether the lane is bisulfite sequence (used when caculating the NM tag).")
        public boolean IS_BISULFITE_SEQUENCE = false;
    @Option(doc="Whether to output only aligned reads.  ") public boolean ALIGNED_READS_ONLY = false;
    @Option(doc="The maximum number of insertions or deletions permitted for an alignment to be " +
            "included.  Alignments with more than this many insertions or deletions will be ignored.  " +
            "Set to -1 to allow any number of insertions or deletions.",
            shortName="MAX_GAPS") public int MAX_INSERTIONS_OR_DELETIONS = 1;
    @Option(doc="Reserved alignment attributes (tags starting with X, Y, or Z) that should be " +
            "brought over from the alignment data when merging.")
    public List<String> ATTRIBUTES_TO_RETAIN = new ArrayList<String>();
    @Option(shortName="R1_TRIM", doc="The number of bases trimmed from the beginning of read 1 prior to alignment")
    public int READ1_TRIM = 0;
    @Option(shortName="R2_TRIM", doc="The number of bases trimmed from the beginning of read 2 prior to alignment")
    public int READ2_TRIM = 0;
    @Option(shortName="ORIENTATIONS", doc="The expected orientation of proper read pairs.  Replaces JUMP_SIZE",
            mutex = "JUMP_SIZE", optional=true)
    public List<SamPairUtil.PairOrientation> EXPECTED_ORIENTATIONS;
    @Option(shortName=StandardOptionDefinitions.SORT_ORDER_SHORT_NAME,
            doc="The order in which the merged reads should be output.")
    public SortOrder SORT_ORDER = SortOrder.coordinate;

    @Option(doc="For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.")
    public boolean CLIP_OVERLAPPING_READS = true;

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
        // TEMPORARY FIX until internal programs all specify EXPECTED_ORIENTATIONS
        if (JUMP_SIZE != null) {
            EXPECTED_ORIENTATIONS = Arrays.asList(new SamPairUtil.PairOrientation[]{SamPairUtil.PairOrientation.RF});
        }
        else if (EXPECTED_ORIENTATIONS == null || EXPECTED_ORIENTATIONS.isEmpty()) {
            EXPECTED_ORIENTATIONS = Arrays.asList(new SamPairUtil.PairOrientation[]{SamPairUtil.PairOrientation.FR});
        }

        SamAlignmentMerger merger = new SamAlignmentMerger (UNMAPPED_BAM, OUTPUT,
            REFERENCE_SEQUENCE, prod, CLIP_ADAPTERS, IS_BISULFITE_SEQUENCE, PAIRED_RUN,
            ALIGNED_READS_ONLY, ALIGNED_BAM, MAX_INSERTIONS_OR_DELETIONS,
            ATTRIBUTES_TO_RETAIN, READ1_TRIM, READ2_TRIM,
            READ1_ALIGNED_BAM, READ2_ALIGNED_BAM, EXPECTED_ORIENTATIONS, SORT_ORDER);
        merger.setClipOverlappingReads(CLIP_OVERLAPPING_READS);
        merger.setMaxRecordsInRam(MAX_RECORDS_IN_RAM);
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

        boolean r1sExist = READ1_ALIGNED_BAM != null && READ1_ALIGNED_BAM.size() > 0;
        boolean r2sExist = READ2_ALIGNED_BAM != null && READ2_ALIGNED_BAM.size() > 0;
        if ((r1sExist && !r2sExist) || (r2sExist && !r1sExist)) {
            return new String[] {"READ1_ALIGNED_BAM and READ2_ALIGNED_BAM " +
                    "must both be supplied or neither should be included.  For " +
                    "single-end read use ALIGNED_BAM."};
        }
        if (ALIGNED_BAM == null || ALIGNED_BAM.size() == 0 && !(r1sExist && r2sExist)) {
            return new String[] {"Either ALIGNED_BAM or the combination of " +
                    "READ1_ALIGNED_BAM and READ2_ALIGNED_BAM must be supplied."};

        }

        return null;
    }

}
