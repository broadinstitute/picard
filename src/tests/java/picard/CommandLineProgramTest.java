package picard;

import picard.cmdline.PicardCommandLine;

import java.util.Arrays;
import java.util.List;

/**
 * Utility class for CommandLine Program testing.
 */
public abstract class CommandLineProgramTest {
    public abstract String getCommandLineProgramName();

    /**
     * For testing support.  Given a name of a Picard CommandLineProgram and it's arguments, builds the arguments appropriate for calling the
     * program through PicardCommandLine
     *
     * @param args
     * @return String[] of command line arguments
     */
    public String[] makePicardCommandLineArgs(final List<String> args) {
        final String[] picardCommandLineArgs = new String[args.size() + 2];
        picardCommandLineArgs[0] = "-T";
        picardCommandLineArgs[1] = getCommandLineProgramName();
        int i = 2;
        for (final String arg : args) {
            picardCommandLineArgs[i++] = arg;
        }
        return picardCommandLineArgs;
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
}
