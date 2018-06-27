package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;

@CommandLineProgramProperties(
        summary = MoveReadNameToTag.USAGE_SUMMARY + MoveReadNameToTag.USAGE_DETAILS,
        oneLineSummary = MoveReadNameToTag.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
public class MoveReadNameToTag extends CommandLineProgram {
    static final String USAGE_SUMMARY = "NONE";
    static final String USAGE_DETAILS = "NONE";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM file")
    public File INPUT;

    @Argument(shortName = "TILE_IN_READ_NAME", doc = "Put only the tile in the readname")
    public Boolean TILE_IN_RN = false;

    @Argument(shortName = "IN_READ_NAME", doc = "Alter read name, not tags")
    public Boolean IN_RN = false;

    @Argument(shortName = "TILE", doc = "Output TILE as tag")
    public Boolean TILE_TAG = false;

    @Argument(shortName = "XY_FULL")
    public Boolean XY_FULL_TAG = false;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output BAM file")
    public File OUTPUT;

    private final Log log = Log.getInstance(MoveReadNameToTag.class);
    private int warnings = 0;
    private int uniqueReadNumber = 0;

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        final ProgressLogger progress = new ProgressLogger(log, 1000000);
        final SamReaderFactory readerFactory = SamReaderFactory.makeDefault();
        final SamReader reader = readerFactory.referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileHeader header = reader.getFileHeader();

        final SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
        final SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, OUTPUT);

        for(SAMRecord rec : reader) {
            String readName = rec.getReadName();
            String[] splitReadName = readName.split(":");
            if(splitReadName.length != 5) {
                log.warn("Suspicious read name:" + readName);
                warnings++;
                if(warnings > 50) {
                    return 1;
                }
                continue;
            }
            if(IN_RN) {
                String newReadName = splitReadName[2] + ":" + splitReadName[3] + ":" + splitReadName[4];
                rec.setReadName(newReadName);
            }
            if(TILE_IN_RN) {
                if(progress.getCount() % 2 == 0) {
                    uniqueReadNumber++;
                }
                String newReadName = String.valueOf(uniqueReadNumber) + ":" + splitReadName[2];
                rec.setReadName(newReadName);
            }
            if(TILE_TAG) {
                rec.setAttribute("XT", Integer.parseInt(splitReadName[2]));
            }
            if(XY_FULL_TAG) {
                rec.setAttribute("XX", Integer.parseInt(splitReadName[3]));
                rec.setAttribute("XY", Integer.parseInt(splitReadName[4]));
            }
            writer.addAlignment(rec);
            progress.record(rec);
        }

        writer.close();
        return 0;
    }
}
