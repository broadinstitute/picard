/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecordDuplicateComparator;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Alpha;

/**
 * This is a simple tool to mark duplicates using the UmiAwareDuplicateSetIterator, DuplicateSet, and SAMRecordDuplicateComparator.
 *
 * It makes use of the UmiAwareDuplicateSetIterator which is a wrapper around the DuplicateSetIterator.  It makes use
 * of the fact that duplicate sets with UMIs can be broken up into subsets based on information contained in the UMI.
 *
 * Users should continue to use MarkDuplicates in general, the main motivation for this tool is to provide a way to
 * mark duplicates using information from UMIs.
 *
 * @author fleharty
 */
@CommandLineProgramProperties(
        usage = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. " +
                "All records are then written to the output file with the duplicate records flagged.",
        usageShort = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules.",
        programGroup = Alpha.class
)

public class UmiAwareMarkDuplicatesWithMateCigar extends SimpleMarkDuplicatesWithMateCigar {

    @Option(shortName = "EDIT_DISTANCE_TO_JOIN",
            doc = "This option specifies the edit distance of UMIs to join", optional = true)
    public int EDIT_DISTANCE_TO_JOIN = 1;


    @Option(shortName = "ADD_INFERRED_UMI",
            doc = "This option adds the inferred UMI to the bam in the inferred UMI tag (by default this tag is RI).", optional = true)
    public boolean ADD_INFERRED_UMI = false;

    @Option(shortName = "UMI_TAG_NAME",
            doc = "Tag name to use for UMI (default is RX)", optional = true)
    public String UMI_TAG_NAME = "RX";

    @Option(shortName = "INFERRED_UMI_TAG_NAME",
            doc = "Tag name to use for inferred UMI (default is RI)", optional = true)
    public String INFERRED_UMI_TAG_NAME = "RI";

    private final Log log = Log.getInstance(UmiAwareMarkDuplicatesWithMateCigar.class);

    /** Stock main method. */
    public static void main(final String[] args) {
        new UmiAwareMarkDuplicatesWithMateCigar().instanceMainWithExit(args);
    }

    @Override
    protected CloseableIterator<DuplicateSet> getDuplicateSetIterator(final SamHeaderAndIterator headerAndIterator, final SAMRecordDuplicateComparator comparator) {
        return new UmiAwareDuplicateSetIterator(
                    new DuplicateSetIterator(headerAndIterator.iterator,
                    headerAndIterator.header,
                    false,
                    comparator), EDIT_DISTANCE_TO_JOIN, ADD_INFERRED_UMI, UMI_TAG_NAME, INFERRED_UMI_TAG_NAME);
    }
}