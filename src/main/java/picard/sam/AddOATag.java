package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.util.List;
import java.util.Optional;

@CommandLineProgramProperties(
        summary = AddOATag.USAGE_DETAILS,
        oneLineSummary = AddOATag.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class AddOATag extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Add the OA tag to reads";
    static final String USAGE_DETAILS = "This tool takes in an aligned SAM or BAM and adds the " +
            "OA tag to every aligned read unless an interval list is specified, where it only adds the tag to reads " +
            "that fall within the intervals in the interval list" +
            "<br />"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar AddOATag \\<br />" +
            "      L=some_picard.interval_list \\<br />" +
            "      I=sorted.bam \\<br />" +
            "      O=fixed.bam <br />"+
            "</pre>";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "SAM or BAM input file")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "SAM or BAM file to write merged result to")
    public File OUTPUT;

    @Argument(shortName = "OV", doc = "Whether or not to override any existing OA tag values", optional = true)
    public Boolean OVERWRITE_TAG = true;

    @Argument(shortName = "L", doc = "An interval list file that contains the locations of reads to add the OA tag to", optional = true)
    public File INTERVAL_LIST;

    private static final Log log = Log.getInstance(AddOATag.class);

    public static void main(final String[] argv) {
        System.exit(new AddOATag().instanceMain(argv));
    }

    @Override
    protected int doWork() {
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);
        writer.setProgressLogger(
                new ProgressLogger(log, (int) 1e7, "Wrote", "records"));

        OverlapDetector overlapDetector = getOverlapDetectorFromIntervalListFile(INTERVAL_LIST, 0, 0);
        for (SAMRecord rec : reader) {
            if (!rec.getReadUnmappedFlag() && (overlapDetector == null || overlapDetector.overlapsAny(rec))) {
                setOATag(rec, OVERWRITE_TAG);
            }
            writer.addAlignment(rec);
        }

        CloserUtil.close(reader);
        writer.close();

        return 0;
    }

    // Take an interval list file and convert it to an overlap detector, can add left and right padding
    static OverlapDetector<Interval> getOverlapDetectorFromIntervalListFile(File intervalList, int lhsBuffer, int rhsBuffer) {
        if (intervalList == null) {
            return null;
        } else {
            List<Interval> intervals = IntervalList.fromFile(intervalList).uniqued().getIntervals();
            OverlapDetector<Interval> detector = new OverlapDetector<>(lhsBuffer, rhsBuffer);
            detector.addAll(intervals, intervals);
            return detector;
        }
    }

    // format OA tag string according to the spec
    private void setOATag(SAMRecord rec, Boolean overrideTag) {
        String OAValue = String.format("%s,%s,%s,%s,%s,%s;",
                (rec.getReferenceName().replace(",", "_")),
                (rec.getAlignmentStart()),
                ((rec.getReadNegativeStrandFlag() ? "-" : "+")),
                (rec.getCigarString()),
                (rec.getMappingQuality()),
                (Optional.ofNullable(rec.getAttribute(SAMTag.NM.name())).orElse("").toString()));
        if (overrideTag) {
            rec.setAttribute("OA", OAValue);
        } else {
            rec.setAttribute("OA",Optional.ofNullable(rec.getAttribute("OA")).orElse("") +  OAValue);
        }
    }
}
