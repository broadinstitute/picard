package picard.cmdline;

import org.apache.commons.io.FileUtils;
import org.testng.annotations.AfterClass;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
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


    // A per-test-class directory that will be deleted after the tests are complete.
    private File tempOutputDir;


    /**
     * returns an directory designated for output which will be deleted after the test class is tested
     */
    public File getTempOutputDir() {
        if (tempOutputDir == null) {
            try {
                tempOutputDir = Files.createTempDirectory(this.getClass().getName()).toFile();
            } catch (IOException e) {
                throw new PicardException("Couldn't create temp directory");
            }
        }
        return tempOutputDir;
    }

    /**
     * returns an file designated for output which will be deleted (together with the entire subdirectory)
     * after the test class is tested
     * @throws IOException when there's a problem creating the file
     */

    public File getTempOutputFile(final String prefix,final String extension) throws IOException {
        return File.createTempFile(prefix, extension, getTempOutputDir());
    }
    @AfterClass
    final void cleanup_temp_dir() throws IOException {
        if (tempOutputDir != null) {
            FileUtils.deleteDirectory(tempOutputDir);
        }
    }

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
