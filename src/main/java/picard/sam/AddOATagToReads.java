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
import java.util.StringJoiner;

@CommandLineProgramProperties(
        summary = AddOATagToReads.USAGE_DETAILS,
        oneLineSummary = AddOATagToReads.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class AddOATagToReads extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Add the OA tag to reads";
    static final String USAGE_DETAILS = "This tool takes in an aligned SAM or BAM and adds the " +
            "OA tag to every aligned read unless an interval list is specified, where it only adds the tag to reads " +
            "that fall within the intervals in the interval list"+
            "<br />"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar AddOATagToReads \\<br />" +
            "      L=some_picard.interval_list \\<br />" +
            "      I=sorted.bam \\<br />" +
            "      O=fixed.bam <br />"+
            "</pre>";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "SAM or BAM input file")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "SAM or BAM file to write merged result to")
    public File OUTPUT;

    @Argument(shortName = "OV", doc = "Whether or not to override any existing OA tags", optional = true)
    public Boolean OVERRIDE_TAGS = true;

    @Argument(shortName = "L", doc = "An interval list file that contains the locations of reads to add the OA tag to", optional = true)
    public File INTERVAL_LIST;

    private static final Log log = Log.getInstance(AddOATagToReads.class);

    public static void main(final String[] argv) {
        System.exit(new AddOATagToReads().instanceMain(argv));
    }

    @Override
    protected int doWork() {
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);
        writer.setProgressLogger(
                new ProgressLogger(log, (int) 1e7, "Wrote", "records"));

        OverlapDetector overlapDetector = getOverlapDetector();
        for (SAMRecord rec : reader) {
            if (!rec.getReadUnmappedFlag() && (overlapDetector == null || overlapDetector.overlapsAny(rec))) {
                setOATag(rec, OVERRIDE_TAGS);
            }
            writer.addAlignment(rec);
        }

        CloserUtil.close(reader);
        writer.close();

        return 0;
    }

    private OverlapDetector<Interval> getOverlapDetector() {
        if (INTERVAL_LIST == null) {
            return null;
        } else {
            List<Interval> intervals = IntervalList.fromFile(INTERVAL_LIST).uniqued().getIntervals();
            OverlapDetector<Interval> detector = new OverlapDetector<>(0, 0);
            detector.addAll(intervals, intervals);
            return detector;
        }
    }

    private void setOATag(SAMRecord rec, Boolean overrideTag) {
        StringJoiner oa_value = new StringJoiner(",");
        oa_value.add(rec.getReferenceName().replace(",", "_"));
        oa_value.add(String.valueOf(rec.getAlignmentStart()));
        oa_value.add((rec.getReadNegativeStrandFlag() ? "-" : "+"));
        oa_value.add(rec.getCigarString());
        oa_value.add(String.valueOf(rec.getMappingQuality()));
        oa_value.add(Optional.ofNullable(rec.getAttribute(SAMTag.NM.name())).orElse("").toString());
        if (overrideTag) {
            rec.setAttribute("OA", oa_value.toString());
        } else {
            rec.setAttribute("OA",Optional.ofNullable(rec.getAttribute("OA")).orElse("") + ";" +  oa_value.toString());
        }
    }
}
