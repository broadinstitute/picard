package picard.sam;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import htsjdk.samtools.util.Log;  // Added this import
import java.io.File;
import java.io.IOException;

@CommandLineProgramProperties(
        summary = "Filters out reads with specific flowcell coordinates",
        oneLineSummary = "Removes reads from specific flowcell positions from BAM/CRAM files",
        programGroup = ReadDataManipulationProgramGroup.class
)
public class FilterFlowCellEdgeReads extends CommandLineProgram {
    // Initialize logger
    private static final Log logger = Log.getInstance(FilterFlowCellEdgeReads.class);

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

    private boolean hasFlowcellCoordinates(String readName) {
        // Parse Illumina read name format
        // Example format: @HWUSI-EAS100R:6:73:941:1973#0/1
        // or: @EAS139:136:FC706VJ:2:2104:15343:197393
        try {
            String[] parts = readName.split(":");
            if (parts.length >= 6) {  // Ensure we have enough parts
                // The last two numbers are typically X and Y coordinates
                int x = Integer.parseInt(parts[parts.length-2]);
                int y = Integer.parseInt(parts[parts.length-1].split("[#/]")[0]);  // Remove any trailing /1 or #0

                return x == X_COORD && y == Y_COORD;
            }
        } catch (NumberFormatException | ArrayIndexOutOfBoundsException e) {
            // If we can't parse the coordinates, assume it doesn't match
            return false;
        }
        return false;
    }

    @Override
    protected int doWork() {
        final SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(REFERENCE_SEQUENCE)
                .open(new File(INPUT));

        final SAMFileHeader header = reader.getFileHeader();
        final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeWriter(header, true, new File(OUTPUT), REFERENCE_SEQUENCE);

        // Process reads
        int totalReads = 0;
        int filteredReads = 0;

        try {
            for (final SAMRecord read : reader) {
                totalReads++;

                // Check if read has the specified flowcell coordinates
                if (hasFlowcellCoordinates(read.getReadName())) {
                    filteredReads++;
                    continue; // Skip this read
                }

                // Write read to output if it doesn't match filter criteria
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
        new FilterFlowCellEdgeReads().instanceMain(args);
    }
}
