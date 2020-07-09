package picard.cmdline;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Utility class for CommandLine Program testing.
 */
public abstract class CommandLineProgramTest {

    public static final File REFERENCE_TEST_DIR = new File("testdata/picard/reference");
    public static final File CHR_M_REFERENCE = new File(REFERENCE_TEST_DIR,"chrM.reference.fasta");
    public static final File CHR_M_DICT = new File(REFERENCE_TEST_DIR,"chrM.reference.dict");

    public abstract String getCommandLineProgramName();

    /**
     * For testing support.  Given a name of a Picard CommandLineProgram and it's arguments, builds the arguments appropriate for calling the
     * program through PicardCommandLine
     *
     * @param programName Picard command line program to run
     * @param args Arguments to pass to command line program
     * @return String[] of command line arguments
     */
    public String[] makePicardCommandLineArgs(final String programName, final List<String> args) {
        final String[] picardCommandLineArgs = new String[args.size() + 1];
        picardCommandLineArgs[0] = programName;
        int i = 1;
        for (final String arg : args) {
            picardCommandLineArgs[i++] = arg;
        }
        return picardCommandLineArgs;
    }

    public String[] makePicardCommandLineArgs(final List<String> args) {
        return makePicardCommandLineArgs(getCommandLineProgramName(), args);
    }

    public String[] makePicardCommandLineArgs(final Map<String, String> kwargs) {
        final List<String> args = new ArrayList<>();
        for (final String key : kwargs.keySet()) {
            args.add(key + "=" + kwargs.get(key));
        }
        return makePicardCommandLineArgs(args);
    }

    public String[] makePicardCommandLineArgs(final String[] args) {
        return makePicardCommandLineArgs(Arrays.asList(args));
    }

    public int runPicardCommandLine(final List<String> args) {
        return new PicardCommandLine().instanceMain(makePicardCommandLineArgs(args));
    }

    public int runPicardCommandLine(final String[] args) {
        return new PicardCommandLine().instanceMain(makePicardCommandLineArgs(args));
    }

    public int runPicardCommandLine(final String programName, final String[] args) {
        return new PicardCommandLine().instanceMain(makePicardCommandLineArgs(programName, Arrays.asList(args)));
    }

    public int runPicardCommandLine(final Map<String, String> kwargs) {
        return new PicardCommandLine().instanceMain(makePicardCommandLineArgs(kwargs));
    }
}
