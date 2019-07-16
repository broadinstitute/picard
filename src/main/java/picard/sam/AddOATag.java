/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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
import htsjdk.samtools.util.*;
import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

@CommandLineProgramProperties(
        summary = AddOATag.USAGE_DETAILS,
        oneLineSummary = AddOATag.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class AddOATag extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Record current alignment information to OA tag.";
    static final String USAGE_DETAILS = "This tool takes in an aligned SAM or BAM and adds the " +
            "OA tag to every aligned read unless an interval list is specified, where it only adds the tag to reads " +
            "that fall within the intervals in the interval list. This can be useful if you are about to realign but want " +
            "to keep the original alignment information as a separate tag." +
            "<br />"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar AddOATag \\<br />" +
            "      L=some_picard.interval_list \\<br />" +
            "      I=sorted.bam \\<br />" +
            "      O=fixed.bam <br />"+
            "</pre>";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "SAM or BAM input file")
    public String INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "SAM or BAM file to write merged result to")
    public String OUTPUT;

    @Argument(shortName = "L", doc = "If provided, only records that overlap given interval list will have the OA tag added.", optional = true)
    public File INTERVAL_LIST;

    private static final Log log = Log.getInstance(AddOATag.class);

    @Override
    protected int doWork() {
        try (final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(IOUtil.getPath(INPUT));
             final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(reader.getFileHeader(), true, IOUtil.getPath(OUTPUT), REFERENCE_SEQUENCE)) {
                writer.setProgressLogger(
                        new ProgressLogger(log, (int) 1e7, "Wrote", "records"));

                final OverlapDetector overlapDetector = getOverlapDetectorFromIntervalListFile(INTERVAL_LIST, 0, 0);
                for (final SAMRecord rec : reader) {
                    if (overlapDetector == null || overlapDetector.overlapsAny(rec)) {
                        setOATag(rec);
                    }
                    writer.addAlignment(rec);
                }
        } catch (IOException e) {
            log.error(e);
            return 1;
        }

        return 0;
    }

    // Take an interval list file and convert it to an overlap detector, can add left and right padding
    static OverlapDetector<Interval> getOverlapDetectorFromIntervalListFile(final File intervalList, final int lhsBuffer, final int rhsBuffer) {
        if (intervalList == null) {
            return null;
        }
        List<Interval> intervals = IntervalList.fromFile(intervalList).uniqued().getIntervals();
        OverlapDetector<Interval> detector = new OverlapDetector<>(lhsBuffer, rhsBuffer);
        detector.addAll(intervals, intervals);
        return detector;
    }

    // format OA tag string according to the spec
    //TODO: Move this to htsjdk once https://github.com/samtools/hts-specs/pull/193 is merged
    private void setOATag(SAMRecord rec) {
        if (rec.getReferenceName().contains(",")) {
            throw new PicardException(String.format("Reference name for record %s contains a comma character.", rec.getReadName()));
        }
        final String oaValue;
        if (rec.getReadUnmappedFlag()) {
            oaValue = String.format("*,0,%s,*,255,;", rec.getReadNegativeStrandFlag() ? Strand.NEGATIVE : Strand.POSITIVE);
        } else {
            oaValue = String.format("%s,%s,%s,%s,%s,%s;",
                    rec.getReferenceName(),
                    rec.getAlignmentStart(),
                    rec.getReadNegativeStrandFlag() ? Strand.NEGATIVE : Strand.POSITIVE,
                    rec.getCigarString(),
                    rec.getMappingQuality(),
                    Optional.ofNullable(rec.getAttribute(SAMTag.NM.name())).orElse(""));
        }
        rec.setAttribute(SAMTag.OA.name(), Optional.ofNullable(rec.getAttribute(SAMTag.OA.name())).orElse("") +  oaValue);
    }
}
