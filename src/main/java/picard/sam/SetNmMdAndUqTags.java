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

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.stream.StreamSupport;

/**
 * Fixes the NM, MD, and UQ tags in a SAM or BAM file.
 *
 * <p>This tool takes in a coordinate-sorted SAM or BAM file and calculates the NM, MD, and UQ
 * tags by comparing with the reference. </p>
 *
 * <p>This may be needed when MergeBamAlignment was run with SORT_ORDER other than 'coordinate'
 * and thus could not fix these tags then. The input must be coordinate sorted in order to run.
 * If specified, the MD and NM tags can be ignored and only the UQ tag be set.</p>
 *
 * <h3>Inputs</h3>
 * <p>
 *     <li> The BAM or SAM file to fix </li>
 *     <li> A reference sequence </li>
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A BAM or SAM output file with recalculated NM, MD, and UQ tags
 * </p>
 *
 * <p>
 * <h4>Usage example: </h4>
 * <h3>Fix the tags in a BAM file:</h3>
 * <pre>
 *     java -jar picard.jar SetNmMdAndUqTags \
 *          R=reference_sequence.fasta \
 *          I=sorted.bam \
 *          O=fixed.bam
 * </pre>
 * </p>
 *
 * @author Yossi Farjoun
 */

@CommandLineProgramProperties(
        summary = SetNmMdAndUqTags.USAGE_DETAILS,
        oneLineSummary = SetNmMdAndUqTags.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class SetNmMdAndUqTags extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Fixes the NM, MD, and UQ tags in a SAM file ";
    static final String USAGE_DETAILS = "This tool takes in a coordinate-sorted SAM or BAM and calculates"+
            "the NM, MD, and UQ tags by comparing with the reference."+
            "<br />" +
            "This may be needed when MergeBamAlignment was run with SORT_ORDER other than 'coordinate' and thus"+
            "could not fix these tags then. The input must be coordinate sorted in order to run. If specified,"+
            "the MD and NM tags can be ignored and only the UQ tag be set."+
            "<br />"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SetNmMdAndUqTags <br />" +
            "      R=reference_sequence.fasta <br />" +
            "      I=sorted.bam <br />" +
            "      O=fixed.bam <br />"+
            "</pre>";

    @Argument(doc = "The BAM or SAM file to fix. ", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "The fixed BAM or SAM output file. ", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "Whether the file contains bisulfite sequence (used when calculating the NM tag).")
    public boolean IS_BISULFITE_SEQUENCE = false;

    @Argument(doc = "Only set the UQ tag, ignore MD and NM.")
    public boolean SET_ONLY_UQ = false;

    @Override
    protected boolean requiresReference() {
        return true;

    }

    private final Log log = Log.getInstance(SetNmMdAndUqTags.class);

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        if (reader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new SAMException("Input must be coordinate-sorted for this program to run. Found: " + reader.getFileHeader().getSortOrder());
        }

        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);
        writer.setProgressLogger(
                new ProgressLogger(log, (int) 1e7, "Wrote", "records"));

        final ReferenceSequenceFileWalker refSeqWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);

        StreamSupport.stream(reader.spliterator(), false)
                .peek(rec -> fixRecord(rec, refSeqWalker))
                .forEach(writer::addAlignment);
        CloserUtil.close(reader);
        writer.close();
        return 0;
    }

    private void fixRecord(SAMRecord record, ReferenceSequenceFileWalker refSeqWalker){
        if (!record.getReadUnmappedFlag()) {
            if (SET_ONLY_UQ) {
                AbstractAlignmentMerger.fixUq(record, refSeqWalker, IS_BISULFITE_SEQUENCE);
            } else {
                AbstractAlignmentMerger.fixNmMdAndUq(record, refSeqWalker, IS_BISULFITE_SEQUENCE);
            }
        }
    }
}
