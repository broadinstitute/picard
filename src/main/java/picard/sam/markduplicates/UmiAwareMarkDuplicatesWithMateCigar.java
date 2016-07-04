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
 * This is a simple tool to mark duplicates making use of UMIs in the reads.
 *
 * It makes use of the fact that duplicate sets with UMIs can be broken up into subsets based on
 * information contained in the UMI.  Since UMIs may contain sequencing errors, this tool allows
 * for UMIs that are different but within a given edit distance to be considered to be part of the
 * same duplicate set.
 *
 * Users should continue to use MarkDuplicates in general, the main motivation for this tool is to provide a way to
 * mark duplicates using information from UMIs.
 *
 * @author fleharty
 */
@CommandLineProgramProperties(
        usage = UmiAwareMarkDuplicatesWithMateCigar.USAGE_SUMMARY + UmiAwareMarkDuplicatesWithMateCigar.USAGE_DETAILS,
        usageShort = UmiAwareMarkDuplicatesWithMateCigar.USAGE_SUMMARY,
        programGroup = Alpha.class
)
public class UmiAwareMarkDuplicatesWithMateCigar extends SimpleMarkDuplicatesWithMateCigar {
    static final String USAGE_SUMMARY = "Identifies duplicate reads using information from read positions and UMIs." +
            "All records are then written to the output file with the duplicate records flagged.";
    static final String USAGE_DETAILS = "<p>UmiAwareMarkDuplicatesWithMateCigar locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are " +
            "defined as originating from a single fragment of DNA. </p>" +
            "<p>This tool identifies a duplicate set by assuming that all members of a duplicate set must have the same start and end position," +
            "and must also have a sufficiently similar UMIs.  Sufficiently similar is parameterized by MAX_EDIT_DISTANCE_TO_JOIN which indicates" +
            "the edit distance between UMIs that shall be considered to be part of the same original molecule.</p>" +
            "<p>This tool is not intended to be used on data without UMIs, see MarkDuplicates for marking duplicates that" +
            "do not have UMIs.</p>";

    @Option(shortName = "MAX_EDIT_DISTANCE_TO_JOIN", doc = "Largest edit distance that UMIs must have in order to be considered as coming from distinct source molecules.", optional = true)
    public int MAX_EDIT_DISTANCE_TO_JOIN = 1;

    @Option(shortName = "UMI_TAG_NAME", doc = "Tag name to use for UMI", optional = true)
    public String UMI_TAG_NAME = "RX";

    @Option(shortName = "ASSIGNED_UMI_TAG", doc = "Tag name to use for assigned UMI", optional = true)
    public String ASSIGNED_UMI_TAG = "MI";

    // Since we inherit from SimpleMarkDuplicatesWithMateCigar, it is useful for us to also inherit the tests
    // which do not contain UMIs.  By default, we don't allow for missing UMIs, but for the inherited tests
    // we allow for missing UMIs.
    @Option(doc = "Allow for missing UMIs if data doesn't have UMIs.  This option is intended to be used only for testing the code.  Use SimpleMarkDuplicatesWithMateCigar if data has missing UMIs.", optional = true)
    public boolean ALLOW_MISSING_UMIS = false;

    private final Log log = Log.getInstance(UmiAwareMarkDuplicatesWithMateCigar.class);

    @Override
    protected CloseableIterator<DuplicateSet> getDuplicateSetIterator(final SamHeaderAndIterator headerAndIterator, final SAMRecordDuplicateComparator comparator) {
        return new UmiAwareDuplicateSetIterator(
                    new DuplicateSetIterator(headerAndIterator.iterator,
                    headerAndIterator.header,
                    false,
                    comparator), MAX_EDIT_DISTANCE_TO_JOIN, UMI_TAG_NAME, ASSIGNED_UMI_TAG, ALLOW_MISSING_UMIS);
    }
}
