package picard.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.illumina.CheckIlluminaDirectory;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.*;

/**
 * Created by farjoun on 9/10/15.
 */
public class PicardCommandLineTest {

    @Test
    public void TestPicardPublic() { // this tests fails if any CLP in picard is missing its @CommandLineProgramProperties annotation
        PicardCommandLine picardCommandLine = new PicardCommandLine();
        picardCommandLine.instanceMain(new String[]{""});
    }

    @Test
    public void testProcessCommandLinePrograms() {
        final List<Class<CommandLineProgram>> allCLPs = new ArrayList<>();

        PicardCommandLine.processAllCommandLinePrograms(
                Collections.singletonList("picard"),
                (Class<CommandLineProgram> clazz, CommandLineProgramProperties clProperties) -> {
                    allCLPs.add(clazz);
                });
        Assert.assertTrue(allCLPs.size() > 0);
        allCLPs.forEach(clazz -> Assert.assertTrue(CommandLineProgram.class.isAssignableFrom(clazz)));
    }

    // Find every CommandLineProgram in picard, instantiate it, and initialize a Barclay parser using the
    // resulting object. This catches any CommandLineProgram that depends on using duplicate (of overrideable) arguments,
    // which are not accepted by the Barclay parser, and have no other tests where this would surface.
    //
    // NOTE that it does NOT test that those tools actually work correctly, only that they aren't immediately
    // rejected with a CommandLineParserInternalException from the Barclay parser.
    @Test
    public void testLaunchAllCommandLineProgramsWithBarclayParser() {
        PicardCommandLine.processAllCommandLinePrograms(
                Collections.singletonList("picard"),
                (Class<CommandLineProgram> clazz, CommandLineProgramProperties clProperties) -> {
                    // Check for missing annotations
                    if (null != clProperties) {
                        try {
                            final Object commandLineProgram = clazz.newInstance();
                            try {
                                new CommandLineArgumentParser(commandLineProgram);
                            } catch (CommandLineException.CommandLineParserInternalException e) {
                                throw new RuntimeException("Barclay command line parser internal exception parsing class: " + clazz.getName(), e);
                            }
                        } catch (IllegalAccessException | InstantiationException e) {
                            throw new RuntimeException("Failure instantiating command line program: " + clazz.getName(), e);
                        }
                    }
                }
        );
    }


    @Test
    public void testPrintUsage() throws Exception {
        PrintStream stdOut = System.out;
        ByteArrayOutputStream stdOutData = new ByteArrayOutputStream();
        PrintStream stdOutProxy = new PrintStream(stdOutData);
        System.setOut(stdOutProxy);

        String preDescription =
                "\u001B[1m\u001B[31mUSAGE:  \u001B[32m<program name>\u001B[1m\u001B[31m [-h]\n" +
                "\n" +
                "\u001B[0m\u001B[1m\u001B[31mAvailable Programs:\n" +
                "\u001B[0m\u001B[37m--------------------------------------------------------------------------------------\n" +
                "\n" +
                "\u001B[0m";
        String checkIlluminaWithDesc = "\u001B[1m\u001B[31mUSAGE:  \u001B[32m<program name>\u001B[1m\u001B[31m [-h]\n" +
                "\n" +
                "\u001B[0m\u001B[1m\u001B[31mAvailable Programs:\n" +
                "\u001B[0m\u001B[37m--------------------------------------------------------------------------------------\n" +
                "\u001B[0m\u001B[31mBase Calling:                                    Tools that process sequencing machine data, e.g. Illumina base calls, and detect sequencing level attributes, e.g. adapters\u001B[0m\n" +
                "\u001B[32m    CheckIlluminaDirectory                       \u001B[33m\u001B[36mAsserts the validity for specified Illumina basecalling data.  \u001B[0m\n" +
                "\n" +
                "\u001B[37m--------------------------------------------------------------------------------------\n" +
                "\n" +
                "\u001B[0m";
        String checkIlluminaOnlyName = "CheckIlluminaDirectory\n";


        Method method = PicardCommandLine.class.getDeclaredMethod("printUsage", Set.class, String.class, Boolean.TYPE, Boolean.TYPE);
        method.setAccessible(true);

        Set<Class<?>> classes = new HashSet<>();
        method.invoke(null, classes, "", false, true);
        Assert.assertEquals(stdOutData.toString(), preDescription);
        stdOutData.reset();
        method.invoke(null, classes, "", true, true);
        Assert.assertEquals(stdOutData.toString(), "");

        classes.add(CheckIlluminaDirectory.class);
        method.invoke(null, classes, "", false, true);
        Assert.assertEquals(stdOutData.toString(), checkIlluminaWithDesc);
        stdOutData.reset();
        method.invoke(null, classes, "", true, true);
        Assert.assertEquals(stdOutData.toString(), checkIlluminaOnlyName);

        System.setOut(stdOut);
    }
}