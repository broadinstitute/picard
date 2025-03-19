package picard.sam;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import htsjdk.samtools.util.Log;
import picard.sam.util.ReadNameParser;
import picard.sam.util.PhysicalLocation;
import java.io.File;
import java.io.IOException;

@CommandLineProgramProperties(
        summary = "Filters out reads with specific flowcell coordinates",
        oneLineSummary = "Removes reads from specific flowcell positions from BAM/CRAM files",
        programGroup = ReadDataManipulationProgramGroup.class
)
public class FilterReadsByFlowCellLocation extends CommandLineProgram {
    private static final Log logger = Log.getInstance(FilterReadsByFlowCellLocation.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Input BAM/CRAM file")
    public String INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output BAM/CRAM file")
    public String OUTPUT;

    @Argument(shortName = "X",
            doc = "X coordinate to filter (default: 1000)",
            optional = true)
    public int X_COORD = 1000;

    @Argument(shortName = "Y",
            doc = "Y coordinate to filter (default: 1000)",
            optional = true)
    public int Y_COORD = 1000;

    private final ReadNameParser readNameParser = new ReadNameParser(ReadNameParser.DEFAULT_READ_NAME_REGEX);

    private boolean hasFlowcellCoordinates(String readName) {
        class ReadLocation implements PhysicalLocation {
            private short libraryId;
            private int x = -1, y = -1; // Default to invalid values
            private short tile;

            @Override
            public void setLibraryId(short libraryId) {
                this.libraryId = libraryId;
            }

            @Override
            public short getLibraryId() {
                return libraryId;
            }

            @Override
            public void setX(int x) {
                this.x = x;
            }

            @Override
            public int getX() {
                return x;
            }

            @Override
            public void setY(int y) {
                this.y = y;
            }

            @Override
            public int getY() {
                return y;
            }

            @Override
            public void setReadGroup(short readGroup) {}

            @Override
            public short getReadGroup() {
                return 0;
            }

            @Override
            public void setTile(short tile) {
                this.tile = tile;
            }

            @Override
            public short getTile() {
                return tile;
            }
        }

        ReadLocation location = new ReadLocation();
        try {
            readNameParser.addLocationInformation(readName, location);
        } catch (Exception e) {
            logger.warn("Failed to parse read name: " + readName, e);
            return false;  // Keep the read if parsing fails
        }

        if (location.getX() == -1 || location.getY() == -1) {
            return false;  // Keep the read if coordinates are invalid
        }

        return location.getX() == X_COORD && location.getY() == Y_COORD;
    }

    @Override
    protected int doWork() {
        final SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(REFERENCE_SEQUENCE)
                .open(new File(INPUT));

        final SAMFileHeader header = reader.getFileHeader();
        final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeWriter(header, true, new File(OUTPUT), REFERENCE_SEQUENCE);

        int totalReads = 0;
        int filteredReads = 0;

        try {
            for (final SAMRecord read : reader) {
                totalReads++;
                if (hasFlowcellCoordinates(read.getReadName())) {
                    filteredReads++;
                    continue;
                }
                writer.addAlignment(read);
            }
        } finally {
            try {
                reader.close();
            } catch (IOException e) {
                logger.error("Error closing input file", e);
            }
            writer.close();
        }

        logger.info("Processed " + totalReads + " total reads");
        logger.info("Filtered " + filteredReads + " reads at flowcell position " + X_COORD + ":" + Y_COORD);
        logger.info("Wrote " + (totalReads - filteredReads) + " reads to output");

        return 0;
    }

    public static void main(String[] args) {
        new FilterReadsByFlowCellLocation().instanceMain(args);
    }
}
