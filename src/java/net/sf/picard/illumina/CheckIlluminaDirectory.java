package net.sf.picard.illumina;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.illumina.parser.*;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.*;

/**
 * Program to check a lane of an Illumina output directory.  This program checks that files exist, are non-zero in length, for every tile/cycle and
 * specified data type.  If NO data type is specified then the default data types used by IlluminaBasecallsToSam are used.
 */
public class CheckIlluminaDirectory extends CommandLineProgram {
    private static final Log log = Log.getInstance(CheckIlluminaDirectory.class);

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() +
            "Check that the files to provide the data specified by DATA_TYPES are available, exist, and are reasonably sized for every tile/cycle.  " +
            "Reasonably sized means non-zero sized for files that exist per tile and equal size for binary files that exist per cycle/per tile. " +
            "CheckIlluminaDirectory  DOES NOT check that the individual records in a file are well-formed.\n";

    @Option(doc="The basecalls output directory. ", shortName="B")
    public File BASECALLS_DIR;

    @Option(doc="The data types that should be available for each tile/cycle.  If this value remains null then the data types that are used in" +
            "IlluminaBaseCallsToSam which is a superset of those used in ExtractIlluminaBarcodes.  These data types vary slightly depending on" +
            "whether or not the run is barcoded so READ_STRUCTURE should be the same as that passed to IlluminaBaseCallsToSam.  Therefore, if you omit this option " +
            "and IlluminaDirIntegrityChecker passes then both those programs should complete UNLESS the individual records of the files themselves are spurious. ",
            shortName="DT",
            optional=true)
    public final Set<IlluminaDataType> DATA_TYPES = new TreeSet<IlluminaDataType>();

    @Option(doc=IlluminaBasecallsToSam.READ_STRUCTURE_DOC + "  Note:  If you want to check whether or not a future IlluminaBasecallsToSam or ExtractIlluminaBarcodes " +
            "run will fail then be sure to use the exact same READ_STRUCTURE that you would pass to these programs for this run.", shortName="RS")
    public String READ_STRUCTURE;

    @Option(doc="Lane number. ", shortName= StandardOptionDefinitions.LANE_SHORT_NAME, minElements = 1)
    public List<Integer> LANES;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CheckIlluminaDirectory().instanceMainWithExit(argv);
    }


    @Override
    protected int doWork() {
        final ReadStructure readStructure = new ReadStructure(READ_STRUCTURE);
        if(DATA_TYPES.size() == 0) {
            DATA_TYPES.addAll(Arrays.asList(IlluminaBasecallsToSam.DATA_TYPES_NO_BARCODE));
        }

        final List<Integer> failingLanes = new ArrayList<Integer>();
        int totalFailures = 0;

        final int [] expectedCycles = new OutputMapping(readStructure).getOutputCycles();
        log.info("Checking lanes(" + StringUtil.join(",", LANES) + " in basecalls directory (" + BASECALLS_DIR.getAbsolutePath() + ")\n");
        log.info("Expected cycles: " + StringUtil.intValuesToString(expectedCycles));

        for(final Integer lane : LANES) {
            final IlluminaFileUtil fileUtil   = new IlluminaFileUtil(BASECALLS_DIR, lane);
            final List<Integer> expectedTiles  = fileUtil.getExpectedTiles();

            log.info("Checking lane " + lane);
            log.info("Expected tiles: "  + StringUtil.join(", ", expectedTiles));

            final int numFailures = verifyLane(fileUtil, expectedTiles, expectedCycles, DATA_TYPES);

            if(numFailures > 0) {
                log.info("Lane " + lane + " FAILED " + " Total Errors: " + numFailures);
                failingLanes.add(lane);
                totalFailures += numFailures;
            } else {
                log.info("Lane " + lane + " SUCCEEDED ");
            }
        }

        int status = 0;
        if(totalFailures == 0) {
            log.info("SUCCEEDED!  All required files are present and non-empty.");
        } else {
            status = totalFailures;
            log.info("FAILED! There were " + totalFailures + " in the following lanes: " + StringUtil.join(", ", failingLanes));
        }

        return status;
    }

    /**
     * Use fileUtil to find the data types that would be used by IlluminaDataProvider.  Verify that for the expected tiles/cycles/data types that all
     * the files needed to provide their data is present.  This method logs every error that is found and returns the number of errors found
     * @param fileUtil A file util paramterized with the directory/lane to check
     * @param expectedTiles The tiles we expect to be available/well-formed
     * @param cycles The cycles we expect to be available/well-formed
     * @param dataTypes The data types we expect to be available/well-formed
     * @return The number of errors found/logged for this directory/lane
     */
    private static final int verifyLane(final IlluminaFileUtil fileUtil, final List<Integer> expectedTiles, final int[] cycles, final Set<IlluminaDataType> dataTypes) {
        if(expectedTiles.size() == 0) {
            throw new PicardException("0 input tiles were specified!  Check to make sure this lane is in the InterOp file!");
        }

        if(cycles.length == 0) {
            throw new PicardException("0 output cycles were specified!");
        }

        int numFailures = 0;

        //find what request IlluminaDataTypes we have files for and select the most preferred file format available for that type
        final Map<IlluminaFileUtil.SupportedIlluminaFormat, Set<IlluminaDataType>> formatToDataTypes = IlluminaDataProviderFactory.determineFormats(dataTypes, fileUtil);

        //find if we have any IlluminaDataType with NO available file formats and, if any exist, throw an exception
        final Set<IlluminaDataType> unmatchedDataTypes = IlluminaDataProviderFactory.findUnmatchedTypes(dataTypes, formatToDataTypes);
        if(unmatchedDataTypes.size() > 0) {
            log.info("Could not find a format with available files for the following data types: " + StringUtil.join(", ", new ArrayList<IlluminaDataType>(unmatchedDataTypes)));
            numFailures += unmatchedDataTypes.size();
        }

        for(final IlluminaFileUtil.SupportedIlluminaFormat format : formatToDataTypes.keySet()) {
            final List<String> failures = fileUtil.getUtil(format).verify(expectedTiles, cycles);
            numFailures += failures.size();
            for(final String failure : failures) {
                log.info(failure);
            }
        }

        return numFailures;
    }


    protected String[] customCommandLineValidation() {
        IoUtil.assertDirectoryIsReadable(BASECALLS_DIR);
        final List<String> errors = new ArrayList<String>();

        for(final Integer lane : LANES) {
            if(lane < 1) {
                errors.add("LANES must be greater than or equal to 1.  LANES passed in " + StringUtil.join(", ", LANES));
                break;
            }
        }

        if(errors.size() == 0) {
            return null;
        } else {
            return errors.toArray(new String[errors.size()]);
        }
    }
}
