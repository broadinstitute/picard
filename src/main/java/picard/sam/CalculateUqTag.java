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
package picard.sam;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.programgroups.SamOrBam;

import java.util.stream.StreamSupport;

/**
 * @author Jose Soto
 */
@CommandLineProgramProperties(
        usage = CalculateUqTag.USAGE_SUMMARY + CalculateUqTag.USAGE_DETAILS,
        usageShort = CalculateUqTag.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class CalculateUqTag extends SetNmMdAndUqTags {
    static final String USAGE_SUMMARY = "Calculate the UQ tag in a SAM file.  ";
    static final String USAGE_DETAILS = "This tool takes in a SAM or BAM file (sorted by coordinate) and calculates the UQ tag by comparing with the reference."+
            "<br />" +
            "This shouold be used when wanting to calculate the UQ tag for every record in your SAM or BAM file"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CalculateUqTag \\<br />" +
            "      I=sorted.bam \\<br />" +
            "      O=UQtagged.bam \\<br />"+
            "</pre>" +
            "<hr />";

    private final Log log = Log.getInstance(CalculateUqTag.class);

    public static void main(final String[] argv) {
        new CalculateUqTag().instanceMainWithExit(argv);
    }

    @Override
    protected int doWork() {
        validateInputs();

        final SamReader reader = createSamReader();

        final SAMFileWriter writer = createSamFileWriter(reader, log);

        final ReferenceSequenceFileWalker refSeq = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);

        StreamSupport.stream(reader.spliterator(), false)
                .peek(rec -> {if (!rec.getReadUnmappedFlag()) AbstractAlignmentMerger.calculateUq(rec, refSeq, IS_BISULFITE_SEQUENCE);})
                .forEach(writer::addAlignment);

        CloserUtil.close(reader);
        writer.close();
        return 0;
    }
}
